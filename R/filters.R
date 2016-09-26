## Filter functions

#' Find rows of a dataframe where any of the values
#' is above a threshold
#' @param df A (numeric) dataframe.
#' @param threshold A number.
#' @return TRUE or FALSE.
#' @export
any_greater_than <- function(df, threshold) {
    apply(df, 1, function(x) any(x > threshold))
}

#' Find rows of a dataframe where any of the values
#' is below a threshold
#' @param df A (numeric) dataframe.
#' @param threshold A number.
#' @return TRUE or FALSE.
#' @export
any_less_than <- function(df, threshold) {
    apply(df, 1, function(x) any(x < threshold))
}

#' Find rows of a dataframe where all of the values
#' are above a threshold
#' @param df A (numeric) dataframe.
#' @param threshold A number.
#' @return TRUE or FALSE.
#' @export
all_greater_than <- function(df, threshold) {
    apply(df, 1, function(x) all(x > threshold))
}

#' Find rows of a dataframe where all of the values
#' are below a threshold
#' @param df A (numeric) dataframe.
#' @param threshold A number.
#' @return TRUE or FALSE.
#' @export
all_less_than <- function(df, threshold) {
    apply(df, 1, function(x) all(x < threshold))
}

#' Find rows of a dataframe where exactly one of the values
#' is above a threshold
#' @param df A (numeric) dataframe.
#' @param threshold A number.
#' @return TRUE or FALSE.
#' @export
one_greater_than <- function(df, threshold) {
    apply(df, 1, function(x) sum(x > threshold) == 1)
}

#' Find rows of a dataframe where exactly one of the values
#' is below a threshold
#' @param df A (numeric) dataframe.
#' @param threshold A number.
#' @return TRUE or FALSE.
#' @export
one_less_than <- function(df, threshold) {
    apply(df, 1, function(x) sum(x < threshold) == 1)
}

#' Attempt to find the boundary between noisy, low-coverage SVs
#' and reliable high-coverage SVs
#' @param df A dataframe.
#' @param index The columns of the dataframe to operate on.
#' @return TRUE or FALSE.
#' @export
estimate_boundaries <- function(df, index, exclude_above = 300) {
    names <- colnames(df)[index]
    boundaries <- list()

    for (i in 1:length(names)) {
        dat <- df[[names[i]]]

        # Ignore very high numbers of reads
        dat <- dat[dat < exclude_above]

        # Estimate cutoff with 2D kmeans
        cl <- kmeans(dat, centers = 2)
        bounds1 <- range(dat[which(cl$cluster == 1)])
        bounds2 <- range(dat[which(cl$cluster == 2)])

        if (all(bounds1 < bounds2)) {
            boundary <- min(bounds2)
        } else {
            boundary <- min(bounds1)
        }
        boundaries[names[i]] <- boundary
    }

    boundaries
}


#' Filter data frame so all samples selected by `index` have read support
#' greater than the boundary threshold, or 0
#' @param df A dataframe
#' @param boundaries Estimated boundaries from estimated_boundaries
#' @param index Columns of df to examine
#' @param lower_grace A 'grace' number of reads: below this many reads an SV is considered absent
#' @param upper_grace A 'grace' number of reads: if an SV is this many reads below the boundary, it still counts as present
#' @return Filtered dataframe
#' @importFrom "stats" kmeans
#' @export
apply_boundaries <- function(df, boundaries, index, lower_grace = 0, upper_grace = 0) {
    # Upper condition matches samples with read coverage above the boundary,
    # less the upper grace value
    upper_condition <- df[, index] > (unlist(boundaries) - upper_grace)

    # Lower condition matches samples with read coverage higher than the
    # lower grace value
    lower_condition <- df[, index] <= lower_grace

    # Return df matching either upper or lower condition
    df[rowSums(upper_condition | lower_condition) == length(index), ]
}


# Remove structural variants that self overlap
#' @importFrom "GenomicRanges" GRanges
#' @importFrom "IRanges" IRanges %outside%
#' @importFrom "S4Vectors" Rle
#' @export
self_overlap <- function(df) {
    # 1. Make Granges objects of left and right
    input_left_ranges <- GRanges(seqnames = S4Vectors::Rle(df[, 1]),
                                 ranges = IRanges(start = as.integer(df[, 3]),
                                                  end = as.integer(df[, 4])))

    input_right_ranges <- GRanges(seqnames = S4Vectors::Rle(df[, 5]),
                                  ranges = IRanges(start = as.integer(df[, 7]),
                                                   end = as.integer(df[, 8])))

    # 2. Pairwise overlap-testing
    df <- df[input_left_ranges %outside% input_right_ranges, ]

    # 3. Output
    return(df)
}

#' @export
#' @importFrom "utils" read.table
#' @importFrom "GenomicRanges" GRanges
#' @importFrom "IRanges" IRanges
#' @importFrom "S4Vectors" Rle
simple_repeats <- function(df, repeat_file){
    repeats <- read.table(repeat_file, header = F)

    # 1. Left Breakpoint
    input_left_ranges <- GRanges(seqnames = Rle(df[, 1]),
                                 ranges = IRanges(start = as.integer(df[, 3]),
                                                  end = as.integer(df[, 4])))

    repeat_ranges <- GRanges(seqnames = Rle(repeats[, 1]),
                             ranges = IRanges(start = as.integer(repeats[, 2]),
                                              end = as.integer(repeats[, 3])))

    # 2. Match with positions in Input
    overlaps_left <- IRanges::findOverlaps(input_left_ranges, repeat_ranges)
    overlaps_left <- as.matrix(overlaps_left)
    colnames(overlaps_left) <- c("Sets", "Repeats")

    # 3. Right Breakpoint
    input_right_ranges <- GRanges(seqnames = Rle(df[, 5]),
                                  ranges = IRanges(start = as.integer(df[, 7]),
                                                   end = as.integer(df[, 8])))

    # 4. Match with positions in Input
    overlaps_right <- IRanges::findOverlaps(input_right_ranges, repeat_ranges)
    overlaps_right <- as.matrix(overlaps_right)
    colnames(overlaps_right) <- c("Sets", "Repeats")

    # 5. Remove the ones which match
    # cat("\n Total removed: ",
    #     round(length(unique(c(overlaps_left[, 1], overlaps_right[, 1]))) / nrow(df), 4) * 100,
    #     "%")
    df <- df[-unique(c(overlaps_left[, 1], overlaps_right[, 1])), ]

    # 6. Output
    return(df)
}

#' Retain SVs that span a distance greater than threshold
#' @param df A dataframe
#' @param threshold Breakpoints must be at least this far apart, or on different chromosomes
#' @return A filtered dataframe
#' @importFrom "GenomicRanges" GRanges
#' @importFrom "IRanges" IRanges
#' @export
variant_distance <- function(df, threshold) {
    left <- GRanges(seqnames = df[, 1],
                    ranges = IRanges(start = as.integer(df[, 3]),
                                     end = as.integer(df[, 4])))
    right <- GRanges(seqnames = df[, 5],
                     ranges = IRanges(start = as.integer(df[, 7]),
                                      end = as.integer(df[, 8])))
    distances <- GenomicRanges::distance(left, right)
    df[which(distances > threshold | is.na(distances)), ]
}
