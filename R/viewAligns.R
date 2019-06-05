#' viewAligns
#'
#' view BAM files in a given range using igvR.
#'
#' @param grange A GRanges object.
#' @param bams PATHs of the BAM files.
#' @param aligns A `GAlignments` or `GAlignmentsList` object.
#' @param genome A character string, one of "hg38", "hg19", "mm10",
#'     "tair10".
#' @param window The windown to extend for the GRanges.
#' @param bamParam The option list for `ScanBamParam`.
#' @param ... Other options for `igvR`.
#' @export
#' @import igvR
viewAligns <- function(grange, bams = NULL, aligns = NULL,
                    genome = "hg19", window = 50,
                    bamParam = list(), ...){
    igv <- igvR(...)
    setGenome(igv, genome)
    Sys.sleep(5)

    if(is.null(aligns)){
        stopifnot(!missing(grange) | !missing(bams))
        stopifnot(is(grange, "GenomicRanges"))
        stopifnot(length(grange) == 1)
        showGenomicRegion(igv, as.character(grange + window))
        ## param <- ScanBamParam(which=grange, what = scanBamWhat(), ...)
        bamParam$which <- grange
        bamParam$what <- scanBamWhat()
        param <- do.call(ScanBamParam, bamParam)
        aligns <- list()
        for(i in seq_along(bams)){
            aligns[[i]] <- readGAlignments(bams[i], use.names=TRUE, param=param)
        }
        names(aligns) <- basename(bams)
    }else{
        if(is(aligns, "GAlignments")){
            aligns <- GAlignmentsList(list(aligns))
        }
        grl <- GRangesList(lapply(aligns, granges))
        grange <- GRanges(unique(unlist(seqnames(grl))),
                          IRanges(min(start(grl)), max(start(grl))))
        showGenomicRegion(igv, as.character(grange + window))
    }

    for(i in seq_along(aligns)){
        ## remove NA in GAlignments
        ga <- aligns[[i]]
        taglist <- scanBamWhat()
        if(!is.null(bamParam$tag)){
            taglist <- c(taglist, bamParam$tag)
        }
        elementMetadata(ga) <- elementMetadata(ga)[,colnames(elementMetadata(ga))
                            %in% taglist]
        ga.na <- na.omit(elementMetadata(ga))
        if(!is.null(na.action(ga.na))){
            ga <- ga[-as.integer(na.action(ga.na))]
        }
        if(is.null(names(aligns)[i])){
            tname <- as.character(i)
        }else{
            tname <- names(aligns)[i]
        }
        track <- GenomicAlignmentTrack(tname, ga)
        displayTrack(igv, track)
        Sys.sleep(5)
    }
}
