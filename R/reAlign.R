#' Local re-align reads
#' 
#' To align reads to reference and alternative sequences based on
#' given indels/variants.
#'
#' @param var A GRanges object with REF and ALT in metadata columns.
#' @param bam a character vector of BAM files.
#' @param genome The path of the reference genome file, or a BSgenome
#'     reference object.
#' @param ... More options for verIndel.
#' @import GenomicFiles
#' @importFrom grDevices col2rgb rainbow
#' @import SummarizedExperiment
#' @importFrom S4Vectors DataFrame
#' @importFrom methods as is
#' @importFrom stats na.action na.omit
#' @export
reAlign <- function(var, bam, genome = BSgenome.Hsapiens.1000genomes.hs37d5,
                    ...){
    stopifnot(is(var, "GRanges"))
    vl <- splitAsList(var, seq_along(var))
    r1 <- reduceByFile(vl, bam, MAP=MAP, genome = genome)

    ## counts
    adp <- lapply(r1, function(y){
        lapply(y, function(x){
            c(colSums(as.matrix(mcols(x)[,-(1:5)]), na.rm = TRUE), DP=length(x))
        })
    })
    adp <- do.call(cbind, adp)

    ## label with YC
    r1_yc <- list()
    for(i in seq(r1)){
        r1_y <- list()
        for(j in seq(r1[[i]])){
            y <- r1[[i]][[j]]
            if(length(y)==0){
                elementMetadata(y)$YC <- character()
            }else{
                elementMetadata(y)$YC <- "128,128,128"
                ref <- as.character(var$REF)[j]
                alts <- as.character(as.list(var$ALT)[[j]])
                seqs <- c(ref, alts)
                cols <- apply(col2rgb(c("grey20", rainbow(length(alts)))),
                              2, paste, collapse=",")
                for(s in seq(seqs)){
                    a1 <- which(elementMetadata(y)[[seqs[s]]])
                    elementMetadata(y)$YC[a1] <- cols[s]
                }
            }
            r1_y[[j]] <- y
        }
        r1_yc[[i]] <- r1_y
    }
    alns <- do.call(cbind, r1_yc)
    ##viewBam(aligns = aln, bamParam=list(tag="YC"))
    rownames(adp) <- rownames(alns) <- names(var)
    colnames(adp) <- colnames(alns) <- basename(bam)
    coldat <- DataFrame(bams=bam)

    SummarizedExperiment(assays = SimpleList(ADP=adp, Align=alns),
                         rowRanges = var,
                         colData = coldat)
}
