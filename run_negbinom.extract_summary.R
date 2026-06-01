options(future.globals.maxSize = 32000 * 1024^2)
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(MASS))
suppressPackageStartupMessages(library(data.table))

parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-j", "--job", help="Integer used as index in the list of gRNAs")
parser$add_argument("-s", "--sgrnas", nargs='+',
                    help="List of folders containing sgRNAs data for Read10X function. Each folder should have 3 files: MarketMatrix file (matrix.mtx.gz); Cell barcodes file (barcodes.tsv.gz); sgRNA identifiers file (features.tsv.gz).")
parser$add_argument("-t", "--transcript-data", nargs='+',
                    help="List of folders containing transcriptomic data for Read10X function. Each folder should have 3 files: MarketMatrix file (matrix.mtx.gz); Cell barcodes file (barcodes.tsv.gz); Gene/transcript identifiers file (features.tsv.gz).")
parser$add_argument("--de-method", default="negbinom", help="DE method label (used in output naming)")
parser$add_argument("--target-region-chrom", help="Chromosome of the target region use for DE analysis")
parser$add_argument("--target-region-start", help="Start coordinate of the target region use for DE analysis")
parser$add_argument("--target-region-end", help="End coordinate of the target region use for DE analysis")
parser$add_argument("--target-region-bed", help="BED file with target region(s) to use for DE analysis")
parser$add_argument("--gene-annotations-file", help="Bed6 file with gene definitions. Do not include headers. Expected columns: c('chrom', 'start', 'end', 'gene_symbol', 'score', 'strand')")
parser$add_argument("--gene-symbol-whitelist", help="Gene symbols (HUGO names) to be used for DE analysis")
parser$add_argument("--output-basename", help="Basename used to create output file (args$output-basename + args$job + .txt)")
parser$add_argument("--cell-barcodes-whitelist", help="When specified, only restrict analysis to these cell barcodes", required=F)
parser$add_argument("--do-not-rename-barcodes", action="store_true", default=FALSE,
                    help="Do not rename cell barcodes (see code)")
parser$add_argument("--pseudocount", help="Pseudocount used to compute log2 fold-changes (Default: 1)", default=1)
parser$add_argument("--latentvar", help="Latent variable to use as covariate in the NB model (e.g. nCount_RNA)")
parser$add_argument("--sctransform", action="store_true", default=FALSE,
                    help="Use SCTransform-ed UMI counts instead of raw UMI counts")
parser$add_argument("-sm", "--sctransform-max-ncell", default=10000, type="double",
                    help="Maximum number of cells used in SCTransform")
parser$add_argument("-smo", "--sctransform-mu-offset", action="store_true", default=FALSE,
                    help="Use SCTransform Mu's as an offset in negative binomial")
parser$add_argument("--rdata",
                    help="RData file to bypass loading GEX and gRNA data and upstream processing before DE analysis")
parser$add_argument("--save-rdata",
                    help="Path to save RData checkpoint after data loading and before DE analysis")
parser$add_argument("--min-cells-per-gene-and-group",
                    help="Minimum number of cells per group expressing a gene to perform DE analysis")
parser$add_argument("--min-percent-cells-across-groups",
                    help="Minimum percentage of cells across groups expressing a gene to perform DE analysis")
parser$add_argument("-f", "--min-nfeature-rna", default=20, type="double", help="Min. nFeature_RNA used to filter cells")
parser$add_argument("-c", "--max-ncount-rna",   default=1500, type="double", help="Max. nCount_RNA used to filter cells (to avoid doublets)")
parser$add_argument("-m", "--max-mt-percent",   default=20, type="double", help="Max. percentage of MT- reads per cell")
parser$add_argument("--sgrnas-binary-mixture-model", action="store_true", default=FALSE, help="Used if the gRNA MarketMatrix file(s) contains binary assignment of guides to cells")
parser$add_argument("--apply-bias-reduction", action="store_true", default=FALSE,
                    help="Apply neg. binomial model with bias reduction (brnb)")

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()
args
job <- as.integer(args$job)

# Quick check of incompatible options
if (!is.null(args$min_cells_per_gene_and_group) && !is.null(args$min_percent_cells_across_groups)){
    stop("--min-cells-per-gene-and-group and --min-percent-cells-across-groups options are mutually exclusive")
}

if (args$apply_bias_reduction){
    suppressPackageStartupMessages(library(brglm2))
}

if (!is.null(args$rdata)){
    args2 <- args
    load(args$rdata)
    args <- args2
    job <- as.integer(args$job)

    # Handle backward compatibility with older Seurat object versions
    cells <- tryCatch({
        print(cells)
        cells
    }, error = function(e) {
        UpdateSeuratObject(cells)
    })

}else{
    gRNAs.counts <- purrr::reduce(purrr::map(args$sgrnas, Read10X, gene.column=1), cbind)

    if (!is.null(args$cell_barcodes_whitelist)){
        cell_barcodes_whitelist <- read.table(args$cell_barcodes_whitelist, sep="\t")$V1
        gRNAs.counts <- gRNAs.counts[, colnames(gRNAs.counts) %in% cell_barcodes_whitelist]
    }

    # The combine function does not assign the '-1' suffix to the first batch of barcodes;
    # append '-1' to all barcodes that lack a batch suffix so they match across assays.
    if (is.null(args$do_not_rename_barcodes)){
        no_suffix <- !grepl("-", colnames(gRNAs.counts))
        colnames(gRNAs.counts)[no_suffix] <- paste0(colnames(gRNAs.counts)[no_suffix],
               rep("-1", sum(no_suffix)))
    }

    # preserve original gRNA names (saved in row names)
    original_rownames <- row.names(gRNAs.counts)

    if (is.null(args$sgrnas_binary_mixture_model) | !args$sgrnas_binary_mixture_model){
        # Create a filter to identify gRNA presence in cells with at least 5 UMIs
        filter_5tags <- gRNAs.counts>=5

        # Divide each UMI count by the total number of UMI counts per cell (analogous to library size)
        gRNAs.counts@x <- gRNAs.counts@x / rep.int(Matrix::colSums(gRNAs.counts), diff(gRNAs.counts@p))

        # Keep only gRNAs: 1) representing >=0.5% of cell UMI counts, and 2) at least 5 UMI counts per gRNA
        gRNAs.counts <- Matrix::drop0((gRNAs.counts>0.005)*filter_5tags)
    }

    # Load 10X transcriptomic reads
    cells.counts <- purrr::reduce(purrr::map(args$transcript_data, Read10X), cbind)

    if (!is.null(args$cell_barcodes_whitelist)){
        cells.counts <- cells.counts[, colnames(cells.counts) %in% colnames(gRNAs.counts)]
    }

    # get genes within target region for DE testing
    if (!is.null(args$target_region_bed)){
        target_region <- GRanges(read.csv(args$target_region_bed, sep='\t', header=F,
                                          col.names=c('chrom', 'start', 'end')))
    } else {
        target_region <- GRanges(data.frame(chrom=c(args$target_region_chrom),
                                            start=c(as.numeric(args$target_region_start)),
                                            end=c(as.numeric(args$target_region_end))))
    }

    # Create Seurat objects
    if (!is.null(names(cells.counts))){
       cells <- CreateSeuratObject(counts=cells.counts$'Gene Expression')
    } else {
       cells <- CreateSeuratObject(counts=cells.counts)
    }
    cat("done converting GEX counts to CreateSeuratObject\n")
    cells.counts <- NULL

    if (!is.null(args$gene_symbol_whitelist)){
        genes <- read.delim(args$gene_symbol_whitelist, header = FALSE)
    } else {
        genes <- GRanges(read.csv(args$gene_annotations_file, header=F,
                                  sep='\t',
                                  col.names=c('chrom', 'start', 'end',
                                              'gene_symbol', 'score', 'strand')))
        genes <- data.frame(V1=unique(genes[queryHits(findOverlaps(genes, target_region, ignore.strand=TRUE)), ]$gene_symbol))
    }

    # Add MT percentage and filter cells:
    #    - fewer than min-nfeature-rna genes detected
    #    - more than max-ncount-rna total UMIs (likely doublets)
    #    - more than max-mt-percent mitochondrial reads
    cells[["percent.mt"]] <- PercentageFeatureSet(cells, pattern = "^MT-")

    expr_nFeatures <- FetchData(cells, vars = 'nFeature_RNA')
    expr_nCounts   <- FetchData(cells, vars = 'nCount_RNA')
    expr_percentMT <- FetchData(cells, vars = 'percent.mt')
    keep_cells <- which(expr_nFeatures > args$min_nfeature_rna & expr_nCounts < args$max_ncount_rna & expr_percentMT < args$max_mt_percent)
    cells <- cells[, keep_cells]

    if (!is.null(args$sctransform) & args$sctransform){
        # Add SCTransformed data layer (use maximum sctransform-max-ncell cells)
        cells <- SCTransform(object = cells,
                             ncells = min(args$sctransform_max_ncell,
                                          dim(cells)[2]))
        # SCTransform will filter out genes not sufficiently expressed across the dataset
        cells <- cells[intersect(row.names(cells[['RNA']]), row.names(cells[['SCT']])),]
    }

    if (!is.null(args$save_rdata)){
        save.image(paste0(args$save_rdata, ".tmp"))
        cat("RData checkpoint saved:", paste0(args$save_rdata, ".tmp"), "\n")
    }

    # Subset object to only interested genes
    cells <- cells[row.names(cells[['RNA']]) %in% genes$V1, ]

    # Add gRNAs data
    cells[["gRNAs"]] <- CreateAssayObject(counts = gRNAs.counts[, keep_cells])
    cat("done adding gRNA counts to CreateSeuratObject\n")
    gRNAs.counts <- NULL

    # LogNormalize data
    cells <- NormalizeData(object = cells, assay = "RNA")

    if (!is.null(args$save_rdata)){
        save.image(args$save_rdata)
        cat("RData saved:", args$save_rdata, "\n")
    }
}
grna_ix <- as.integer(args$job)
cat(grna_ix)
cat("\nNumber of cells with gRNA: ")
cat(length(which(as.vector(GetAssayData(cells[['gRNAs']])[grna_ix,]>0))))
cat("\n")
Idents(object = cells) <- 0
Idents(object = cells, cells = which(as.vector(GetAssayData(cells[['gRNAs']])[grna_ix,]>0))) <- 1
Idents(object = cells, cells = which(as.vector(colSums(GetAssayData(cells[['gRNAs']]))==0))) <- 2

mast_results_selected <- NA

    # Subset data
    if (!is.null(args$sctransform) & args$sctransform){
        data.use <-  GetAssayData(object = cells, assay='SCT', layer = 'counts')
    } else {
        data.use <-  GetAssayData(object = cells, assay='RNA', layer = 'counts')
    }
    cellnames.use <-  colnames(x = data.use)

    # Paste IdentsToCells
    cells_tmp <- IdentsToCells(
        object = cells,
        ident.1 = 1,
        ident.2 = 0,
        cellnames.use = cellnames.use
      )

    # Subset data to use
    data.use <- data.use[, c(cells_tmp$cells.1, cells_tmp$cells.2)]

    # add latent variables
    if (!is.null(args$latentvar)){
        latent.vars <- FetchData(
              object = cells,
              vars = args$latentvar,
              cells = c(cells_tmp$cells.1, cells_tmp$cells.2)
            )
    }else{
        latent.vars <- NULL
    }

    # Construct contrast based on group + latent_vars
    group.info <- data.frame(row.names = c(cells_tmp$cells.1, cells_tmp$cells.2))
    latent.vars <- latent.vars %||% group.info
    group.info[cells_tmp$cells.1, "group"] <- "Group1"
    group.info[cells_tmp$cells.2, "group"] <- "Group2"
    group.info[, "group"] <- factor(x = group.info[, "group"])
    latent.vars.names <- c("condition", colnames(x = latent.vars))
    latent.vars <- cbind(latent.vars, group.info)


    genes_ix_to_test <- 1:nrow(data.use)
    if (!is.null(args$min_cells_per_gene_and_group)){
        genes_min_filter <- lapply(1:nrow(data.use), function(x){
                tmp <- latent.vars
                tmp[, "GENE"] <- as.numeric(x = data.use[x, ])
                sum(tmp[tmp$group == "Group1", "GENE"] >0)>as.numeric(args$min_cells_per_gene_and_group) &
                sum(tmp[tmp$group == "Group2", "GENE"] >0)>as.numeric(args$min_cells_per_gene_and_group)
        })
        genes_ix_to_test <- which(unlist(genes_min_filter))
    }
    if (!is.null(args$min_percent_cells_across_groups)){
        min_percent_cells <- as.numeric(args$min_percent_cells_across_groups)
        if (min_percent_cells > 1){
            min_percent_cells <- min_percent_cells / 100
        }
        ncells_thres <- round(ncol(cells) * min_percent_cells)
        genes_min_filter <- lapply(1:nrow(data.use), function(x){
            tmp <- latent.vars
            tmp[, "GENE"] <- as.numeric(x = data.use[x, ])
            sum(tmp[, "GENE"] >0)>ncells_thres
        })
        genes_ix_to_test <- which(unlist(genes_min_filter))
    }

    if(!is.null(args$sctransform_mu_offset) & args$sctransform_mu_offset){
        regressor_data_orig <- model.matrix(
            as.formula(gsub('^y', '', cells[['SCT']]@misc$vst.out$model_str)),
            cells[['SCT']]@misc$vst.out$cell_attr
        )
        coefs <- cells[['SCT']]@misc$vst.out$model_pars_fit[,-1]
        mu <- exp(tcrossprod(coefs, regressor_data_orig)) %>% t
        thetas <- cells[['SCT']]@misc$vst.out$model_pars_fit %>% as.data.frame
    }

    # Run negbinom across all genes
    start.time <- Sys.time()
    foo_res <- lapply(genes_ix_to_test, function(x){
        tryCatch({
            fmla2 <- paste(c("GENE ~ group", args$latentvar), collapse = " + ")
            tmp <- latent.vars
            tmp[, "GENE"] <- as.numeric(x = data.use[x, ])
            if (!is.null(args$sctransform_mu_offset) && args$sctransform_mu_offset){
                tmp[rownames(mu), "mu"] = mu[,row.names(data.use)[x]]
                fmla2 <- paste(c("GENE ~ group", args$latentvar, "mu"), collapse = " + ")
            }

            if(args$apply_bias_reduction){
                negbinom.res <- brnb(formula = fmla2, data = tmp, link="log",
                                     transformation="inverse",type="MPL_Jeffreys")
                coef_rows <- if (!is.null(args$sctransform_mu_offset) && args$sctransform_mu_offset) c(1,2,3) else c(1,2)
                foo2 <- as.data.frame(summary(negbinom.res)$coefficients[coef_rows,])
            } else {
                negbinom.res <- glm.nb(formula = fmla2, data = tmp)
                coef_rows <- if (!is.null(args$sctransform_mu_offset) && args$sctransform_mu_offset) c(1,2,3) else c(1,2)
                foo2 <- as.data.frame(summary(negbinom.res)$coef[coef_rows,])
            }
            foo2$gene <- row.names(data.use)[x]
            foo2$grna <- row.names(GetAssayData(cells[['gRNAs']]))[grna_ix]
            foo2$negbinom_component <- row.names(foo2)
            foo2
        }, error = function(e) {
            NA
        })
    })
    foo_res <- data.table::rbindlist(foo_res[!is.na(foo_res)])
    end.time <- Sys.time()
    time.taken <- end.time - start.time
    time.taken

write.table(foo_res,
            paste(args$output_basename, args$job, "txt", sep="."),
            sep='\t', quote=F,  row.names=F)
