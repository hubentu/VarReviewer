
readsOut <- function(reads, bam.out){
    adp <- lapply(reads, function(x){
        c(colSums(as.matrix(mcols(x)[,-(1:5)]), na.rm = TRUE), DP=length(x))
    })

    if(!missing(bam.out)){
        reads <- lapply(reads, function(x){
            F <- as.list(mcols(x)[,-c(1:5)])
            for(i in 1:length(F)){
                names(x)[which(F[[i]])] <- paste(names(x)[which(F[[i]])],
                                                 names(F)[i], sep=":")
            }
            mcols(x) <- mcols(x)[,1:3]
            x
        })
        reads <- unlist(GAlignmentsList(reads), use.names = FALSE)
        rtracklayer::export(reads, BamFile(bam.out), format="bam")
    }
    return(adp)
}
## reads <- unlist(Res[[1]])
## readsOut(reads)

## MAP <- function(range, file, ...){
##     lapply(seq(range),
##            function(x){
##                verIndel(range=range[x], bam=file, ...)
##            })
## }

## MAP <- function(range, file, ...){
##     verIndel(range, bam = file, ...)
## }

VerINDEL <- function(vranges, bams, outDir, genome=BSgenome.Hsapiens.1000genomes.hs37d5, reduceBy=c("file", "range"),  ...){
    reduceBy <- match.arg(reduceBy, c("file", "range"))
    if(reduceBy=="file"){
        readsList <- reduceByFile(vranges, files=bams, MAP=MAP,
                                  iterate = FALSE, genome=genome, ...)
        readsList <- lapply(readsList, unlist)
    }else{
        readsList <- reduceByRange(vranges, files=bams, MAP=MAP,
                                   iterate = FALSE, genome=genome, ...)
        readsList <- lapply(readsList, unlist)
        readsDF <- do.call(cbind, readsList)
        readsList <- split(readsDF, 1:nrow(readsDF))
    }
    if(!missing(outDir)){
        if(!exists(outDir)){
            dir.create(outDir, showWarnings = FALSE)
        }
        outBams <- file.path(outDir, basename(bams))
        adpList <- lapply(seq(readsList), function(i){
            readsOut(readsList[[i]], bam.out=outBams[i])
        })
        coldat <- S4Vectors::DataFrame(bams=bams, outBams=outBams)
    }else{
        adpList <- lapply(seq(readsList), function(i){
            readsOut(readsList[[i]])
        })
        coldat <- S4Vectors::DataFrame(bams=bams)
    }
    ADP <- do.call(cbind, adpList)
    reads <- do.call(cbind, readsList)
    rownames(ADP) <- rownames(reads) <- names(vranges)
    colnames(ADP) <- colnames(reads) <- basename(bams)
    SummarizedExperiment(assays = SimpleList(ADP=ADP, reads=reads),
                         rowRanges = vranges,
                         colData = coldat)
}

## Summarize
