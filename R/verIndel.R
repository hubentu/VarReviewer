#' Verify Indel alignments
#' @importFrom S4Vectors elementMetadata
#' @import VariantAnnotation
#' @export

verIndel <- function(range, Ref = NULL, Alt = NULL, bam, genome=BSgenome.Hsapiens.1000genomes.hs37d5, mapq=20, baseQ13=NULL, trim.softclip=FALSE){
    if(is(range, "CollapsedVCF") && (missing(Ref) | missing(Alt))){
        Ref <- NULL; Alt <- NULL
        vrange <- range[1,1]
    }else if(is(range, "GRanges") && (missing(Ref) | missing(Alt))){
        if(all(c("REF", "ALT") %in% colnames(elementMetadata(range)))){
            ##vrange <- VCF(vrange, fixed=mcols(vrange))
            vrange <- VCF(range, fixed=S4Vectors::DataFrame(REF=DNAStringSet(elementMetadata(range)$REF),
                                                  ALT=DNAStringSetList(elementMetadata(range)$ALT)))
        }else{
            ##return(range)
            stop("ref and alt missing for GRange input")
        }
    }else if(is(range, "GRanges") & !missing(Ref) & !missing(Alt)){
        vrange <- VCF(range, fixed=S4Vectors::DataFrame(REF=DNAStringSet(Ref), ALT=DNAStringSetList(Alt)))
    }else{
        stop("Ref and Alt missing for GRange input")
    }

    v1e <- VariantAnnotation::expand(vrange)
    ## check repeat regions
    v1rg <- lapply(v1e, repReg, genome)
    v1rg <- GenomicRanges::reduce(unlist(GRangesList(v1rg)))

    ## load reads
    param <- ScanBamParam(which=v1rg,
                          flag=scanBamFlag(isDuplicate=FALSE),
                          what=c("seq", "mapq", "qual"))
    reads <- readGAlignments(bam, param=param, use.names=TRUE)

    ## filter reads with qualities
    reads <- filterReads(reads, baseQ13=baseQ13, mapq=mapq)

    ## trimsoftclip reads
    reads <- trimSoftclip(reads, trim=trim.softclip)

    if(length(reads)>0){
        ## ref/alt alignment
        aln <- RefAltAln(vrange, reads, v1rg, genome)

        ## mismatch tables
        f.ref <- AlnCheck(aln$Ref, aln$rRange)
        f.alt <- lapply(1:length(aln$Alt), function(x){
            AlnCheck(aln$Alt[[x]], range=aln$aRange[x])
        })
        F <- cbind(f.ref, do.call(cbind, f.alt))
        colnames(F) <- as.character(c(VariantAnnotation::ref(vrange),
                                      unlist(VariantAnnotation::alt(vrange))))
        ## update cigar
        reads <- mCigar(aln, F, reads)
    }else{
        F <- data.frame(matrix(logical(), 0,
                               lengths(VariantAnnotation::alt(vrange)) + 1))
        colnames(F) <- as.character(c(VariantAnnotation::ref(vrange),
                                      unlist(VariantAnnotation::alt(vrange))))
    }

    ## ## export BAM
    ## if(!is.null(bam.out)){
    ##     names(reads)[which(f.ref)] <- paste(names(reads)[which(f.ref)], colnames(F)[1], sep=":")
    ##     for(i in 1:length(f.alt)){
    ##         names(reads)[which(f.alt[[i]])] <- paste(names(reads)[which(f.alt[[i]])],
    ##                                                  colnames(F)[i+1], sep=":")
    ##     }
    ##     rtracklayer::export(reads, BamFile(bam.out), format="bam")
    ## }

    mcols(reads) <- S4Vectors::DataFrame(mcols(reads), F)
    return(reads)
}

#' MAP function for reduceByFile
#' 
MAP <- function(range, file, ...){
    verIndel(range, bam = file, ...)
}
