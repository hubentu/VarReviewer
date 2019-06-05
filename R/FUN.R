#' Class to add ref/alt to GRanges
#' 
#' @importFrom Biostrings DNAStringSet DNAStringSetList
#' @importClassesFrom Biostrings DNAStringSet DNAStringSetList
## setClass("CollapsedVRanges", contains="GRanges",
##          prototype=prototype(
##              S4Vectors::DataFrame(REF=DNAStringSet(), ALT=DNAStringSetList())
##          ))
## setMethod("ref", "CollapsedVRanges",
##           function(x){
##               mcols(x)$REF
##           })
## setMethod("alt", "CollapsedVRanges",
##           function(x){
##               mcols(x)$ALT
##           })

.checkExRanges <- function(vrange){
    c1 <- is(vrange, "ExpandedVRanges")
    c2 <- is(vrange, "ExpandedVCF")
    c3 <- length(vrange)==1
    ## stopifnot(c1 | c2 && c3)
    if(!(c1 | c2 | c3)){
        stop("vrange must be ExpandedVRanges or ExpandedVCF")
    }else if(!c3){
        stop("the length of vrange must be 1")
    }
}

#' Extract reference sequence
#'
#' Get the segments of reference sequence from BSgenome library or a fastq file
#'
#' @param genome a character string for fasta file location, or a BSgenome reference object
#' @param gr a GRange object indicate the range to retrieve
getReferenceSeq <- function(gr, genome=BSgenome.Hsapiens.UCSC.hg38){
    if(is.character(genome)){ # assume to be the reference fasta file location
    # read from .fa file, return A DNAStringSet instance of length 1
        scanFa(genome, param=gr)
    }else if(is(genome, "BSgenome")){
        suppressWarnings(seqlevelsStyle(gr) <- seqlevelsStyle(genome)[1])
        seq <- getSeq(genome, seqnames(gr), start=start(gr), width=width(gr), strand=strand(gr), as.character=FALSE) # a "DNAString" instance
      # the following convertion may not be absolute necessary for our program
        ##seq <- DNAStringSet(seq)
        ##names(seq) <- seqnames(gr)
        return(seq)
    }else{
        stop("Not supported yet")
    }
}

#' check repeat region
repReg <- function(vrange, genome, check.range=NULL){
    .checkExRanges(vrange)
    chr <- as.character(seqnames(vrange))
    Ref <- VariantAnnotation::ref(vrange)
    Alt <- VariantAnnotation::alt(vrange)
    pos <- start(vrange)

    if(width(Ref) > width(Alt)){ # if ref bigger than alt, deletion, get characters that remain after taking alt from ref
        alta <- sub(Alt, "", Ref) #del
    }else{ # if alt bigger than ref, insertion, get characters that remain after taking ref from alt
        alta <- sub(Ref, "", Alt) #in
    }
    
    ##extract all overlapped reads
    if(is.null(check.range)){
        ##mw <- max(nchar(ref), nchar(alt))
        ##which <- GRanges(chr, IRanges(pos+1, width=mw-1))
        
        ##extract Ref and Alt: A-repeat*-B
        i <- nchar(alta) # size of indel

        repeat{
            which <- GenomicRanges::GRanges(chr, IRanges::IRanges(pos, pos+i))
            ## refseq <- scanFa(genome, param=which)
            refseq <- getReferenceSeq(which, genome=genome)
            subseq <- sub(paste("(", alta, ")*$", sep=""), "", as.character(refseq))
            if(nchar(subseq)>1){
                len <- i - nchar(alta) + 1
                break
            }
            i <- i + nchar(alta)
        }
        which <- GRanges(chr, IRanges(pos, pos+len))
    }else{
        if(is(check.range, "GRanges")){
            which <- check.range
        }else if(is(check.range,"IRanges")){
            which <- GRanges(chr, check.range)
        }else{
            stop("check.range should be a range data")
        }
        ##indelGR <- GRanges(chr, IRanges(pos+1, width=nchar(alta)))
        ##depth <- length(unique(findOverlaps(reads, indelGR)@queryHits))
    }
    return(which)
}


#' filter out reads based on mapq and LQ20
#'
#' @param reads The sequence of reads, a GAlignments object.
#' Assumed to have a column qual of PhredQuality objects or a column mapq of mapping quality.
#' @param baseQ13  Choose reads with the number of at most bad base quality (< 13). default NULL
#' @param mapQ  Remove all reads with quality less than mapq. Default 20. Set to NULL to keep
#' all reads.
#'
#' @return The altered reads.
#' @export
filterReads <- function(reads, baseQ13=NULL, mapq=20){
    if(!is.null(mapq)){
        reads <- reads[mcols(reads)$mapq >= mapq]
    }
    if(!is.null(baseQ13)){
        if(length(unique(qwidth(reads)))>1){
            qint <- as(mcols(reads)$qual, "IntegerList")
            idx.lq <- sum(qint < 13) <= baseQ13
            reads <- reads[idx.lq]
        }else{
            q <- as.matrix(mcols(reads)$qual)
            reads <- reads[rowSums(q < 20) <= baseQ13]
        }
    }
    return(reads)
}

trimSoftclip <- function(reads, trim=FALSE){
    ##fix positions for reads with softclip
    rstart <- start(reads)
    rend <- end(reads)

    if(trim){ #remove softclip parts
        seqs <- mcols(reads)$seq
        sidx <- grep("^\\d*S", cigar(reads))
        seqs[sidx] <- DNAStringSet(seqs[sidx],
                                   start=as.numeric(sub("(^[0-9]*)S.*", "\\1", cigar(reads)[sidx]))+1)
        sidx <- grep("\\d*S$", cigar(reads))
        width <- width(seqs[sidx]) - as.numeric(sub(".*([0-9]*[A-Z])([0-9]*)S$", "\\2", cigar(reads)[sidx]))
        seqs[sidx] <- DNAStringSet(seqs[sidx], start=1, width=width)
        mcols(reads)$seq <- seqs
    }else{ #keep softclip parts
        sidx <- grep("^\\d*S", cigar(reads))
        rstart[sidx] <- rstart[sidx] - as.numeric(sub("(^[0-9]*)S.*", "\\1", cigar(reads)[sidx]))
        sidx <- grep("\\d*S$", cigar(reads))
        rend[sidx] <- rend[sidx] + as.numeric(sub(".*([0-9]*[A-Z])([0-9]*)S$", "\\2", cigar(reads)[sidx]))
    }
    mcols(reads)$rstart <- rstart
    mcols(reads)$rend <- rend
    ## pos.start <- min(rstart)
    ## pos.end <- max(rend)
    return(reads)
}


## Alignment
#' Construct alternative reference sequence
#'
#' To construct an alternative referecen sequence.
#' This depend on:
#'     1. the standard and alternative sequence
#'     2. check region, in particular the start position
#'     3. the standard sequence
#'
#' ISSUE: in repeatRegion()/repReg() we had calculated alta. relationship?
#'
#' @export
#' 

RefAltAln <- function(vrange, reads, check.range, genome=BSgenome.Hsapiens.UCSC.hg38){
    if(!all(c("rstart", "rend") %in% colnames(mcols(reads)))){
        stop("run trimSoftclip first")
    }
    pos.start <- min(mcols(reads)$rstart)
    pos.end <- max(mcols(reads)$rend)
    refGR <- GRanges(seqnames(vrange), IRanges(pos.start, pos.end))
    ## Reference seqs
    Ref <- getReferenceSeq(refGR, genome)
    ## variant relative position
    pos <- start(vrange)
    pos.rstart <- pos - pos.start + 1
    pos.rend <- pos - pos.start + width(check.range)
    ## Alt seqs
    ## Alt <- unlist(c(subseq(Ref, 1, pos.rstart-1),
    ##                 DNAStringSet(sub(ref, alt(vrange), as.character(subseq(Ref, pos.rstart))))))
    alt.start <- DNAStringSet(subseq(Ref, 1, pos.rstart-1))
    alt.end <- lapply(unlist(alt(vrange)), function(x)sub(ref(vrange), x, subseq(Ref, pos.rstart)))
    Alt <- lapply(alt.end, function(x)unlist(c(alt.start, DNAStringSet(x))))
    
    ref.aln <- pairwiseAlignment(mcols(reads)$seq, Ref,
                                 type="overlap", gapOpening=5, gapExtension=2)
    alt.aln <- lapply(Alt, function(x){
        pairwiseAlignment(mcols(reads)$seq,
                          x, type="overlap", gapOpening=5, gapExtension=2)})
    ## alt.aln <- pairwiseAlignment(mcols(reads)$seq, Alt,
    ##                              type="overlap", gapOpening=5, gapExtension=2)
    pos.rend.alt <- pos.rend - (width(ref(vrange)) - width(unlist(alt(vrange))))
    
    return(list(Ref=ref.aln, Alt=alt.aln,
                rRange=IRanges(pos.rstart, pos.rend),
                aRange=IRanges(pos.rstart, unlist(pos.rend.alt))))
}

## mismatch check for Ref Or Alt
AlnCheck <- function(alignment, range, over.bp=1, window=5, m=1){
    pos.rstart <- start(range)
    pos.rend <- end(range)
    ## overlap
    f1.ol <- start(Views(alignment)) <= pos.rstart-(over.bp-1) &
        end(Views(alignment)) >= pos.rend+(over.bp-1)
    ## no indel
    f2.indel <- lengths(insertion(alignment))==0 & lengths(deletion(alignment))==0
    ## mismatch
    mt <- mismatchTable(alignment)
    if(!is.null(window)){
        mt <- mt[mt$SubjectStart >= pos.rstart-window & mt$SubjectStart <= pos.rend+window,,drop=FALSE]
    }
    if(m>0){
        m.count <- table(factor(mt$PatternId, levels=1:length(alignment)))
        f3.mt <- m.count <= m
        ## fix bug for SNP, start position can't be mismatch
        f3.mt[mt$PatternId[mt$SubjectStart == pos.rstart]] <- FALSE
    }

    f.aln <- f1.ol & f2.indel & f3.mt
    f.unsure <- f2.indel & f3.mt & !f1.ol #not overlap but well aligned
    f.aln[f.unsure] <- NA
    return(f.aln)
}


mCigar <- function(aln, F, reads){
    varLen <- nchar(colnames(F))
    varLen <- varLen[1] - varLen[-1]
    for(i in 1:(ncol(F)-1)){
        f1 <- which(F[,i+1])
        if(varLen[i] > 0){
            ## insertion
            match1 <- start(aln$aRange)[i] - start(Views(aln$Alt[[i]]))[f1] + 1
            match2 <- width(mcols(reads)$seq)[f1] - match1
            reads@cigar[f1] <- paste0(match1, "M", varLen[i], "D", match2, "M")
            ##start1 <- start(vrange) - match1 + 1
        }else if(varLen[i] < 0){
            ## deletion
            match1 <- start(aln$aRange)[i] - start(Views(aln$Alt[[i]]))[f1] + 1
            match2 <- width(mcols(reads)$seq)[f1] - match1 + varLen[i]
            reads@cigar[f1] <- paste0(match1, "M", -varLen[i], "I", match2, "M")
            ##start2 <- start(vrange) - match1 + 1
        }
    }
    return(reads)
}


#' Check ref/alt counts for given GAlignments
#'
#' 
GAlnCheck <- function(galn, pos, ref, alt){
    reads <- elementMetadata(galn)$seq
    ## mismatch
    if(nchar(ref)==1 & all(unlist(nchar(alt))==1)){
        p1 <- pos - start(galn) + 1
        f1 <- p1 > 0
        p1[!f1] <- 1
        F1 <- DNAStringSet(reads, start=p1, width=1) == ref
        F1[!f1] <- NA
    }

}
