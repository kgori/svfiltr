# Author: Max Stammnitz

#' @export
#' @importFrom "utils" read.table
annotate.genes <- function(x, genes_file){

    # 1. Load file with gene annotations
    genes <- read.table(genes_file, header = T, sep = "\t", comment.char = "")
    genes[, 3] <- gsub("chrM", "MT", genes[, 3])
    genes[, 3] <- gsub("chr", "", genes[, 3])

    # 2. Check First breakpoint
    BP1.Ranges <- GenomicRanges::GRanges(seqnames = S4Vectors::Rle(x[, 1]),
                                         ranges = IRanges::IRanges(start = as.numeric(x[, 3]),
                                                                   end = as.numeric(x[, 4])))
    Genes.Ranges <- GenomicRanges::GRanges(seqnames = S4Vectors::Rle(genes[, 3]),
                                           ranges = IRanges::IRanges(start = genes[, 7],
                                                                     end = genes[, 8]))
    OL <- IRanges::findOverlaps(BP1.Ranges, Genes.Ranges)
    OL <- as.matrix(OL)
    colnames(OL) <- c("Coordinates", "Genes")
    OL[, 2] <- paste(as.character(genes[OL[, 2], 2]))

    # 3. Add hits to end
    add <- rep("", nrow(x))
    add[as.numeric(OL[, 1])] <- as.character(OL[, 2])
    x <- cbind(x, add)
    colnames(x)[ncol(x)] <- "BP1"

    # 4. Check First breakpoint
    BP2.Ranges <- GenomicRanges::GRanges(seqnames = S4Vectors::Rle(x[, 5]),
                                         ranges = IRanges::IRanges(start = as.numeric(x[, 7]),
                                                                   end = as.numeric(x[, 8])))
    OL <- IRanges::findOverlaps(BP2.Ranges, Genes.Ranges)
    OL <- as.matrix(OL)
    colnames(OL) <- c("Coordinates", "Genes")
    OL[, 2] <- paste(as.character(genes[OL[, 2], 2]))

    # 5. Add hits to end
    add <- rep("", nrow(x))
    add[as.numeric(OL[, 1])] <- as.character(OL[, 2])
    x <- cbind(x, add)
    colnames(x)[ncol(x)] <- "BP2"

    # 6. Output
    return(x)
}
