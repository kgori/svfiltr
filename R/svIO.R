extract_lumpy_read_count <- function(df, index) {
    pattern <- "^./.:([[:digit:]]+):([[:digit:]]+):([[:digit:]]+)$"
    col <- df[, index]
    prefix = colnames(df)[index]
    l <- regmatches(col, regexec(pattern, col))
    total <- as.integer(unlist(lapply(l, '[[', 2)))
    disc <- as.integer(unlist(lapply(l, '[[', 3)))
    split <- as.integer(unlist(lapply(l, '[[', 4)))
    total_col = paste(prefix, "-Total", sep="")
    disc_col = paste(prefix, "-Discordant", sep="")
    split_col = paste(prefix, "-Split", sep="")
    out <- data.frame(total_col = total, disc_col = disc, split_col = split)
    colnames(out) <- c(total_col, disc_col, split_col)
    out
}

#' Returns a dataframe from a hydra-multi run
#' @export
load_hydra <- function(filename, header = FALSE) {
    #!#!#!#!#!#!# TODO
    # Refactor this section - common to load_brass
    if (!file.exists(filename)) {
        stop("File not found")
    }

    # 0. Preliminary file exam - see how many samples we've got
    GZIPPED <- FALSE

    if (stringi::stri_endswith_fixed(filename, ".gz")) {
        conn <- gzfile(filename)
        GZIPPED <- TRUE
    } else {
        conn <- file(filename)
    }
    open(conn)
    line <- readLines(conn, 1)  # read first line
    headskipcount <- 0
    while (grepl("^#|^%", line)) {
        line <- readLines(conn, 1)
        headskipcount <- headskipcount + 1
    }
    split <- strsplit(line, "\t")[[1]]
    nFields <- length(split)
    #!#!#!#!#!#!# - Refactor

    nsamples <- nFields - 24
    coltypes <- readr::cols(readr::col_character(),
                            readr::col_integer(),
                            readr::col_integer(),
                            readr::col_character(),
                            readr::col_integer(),
                            readr::col_integer(),
                            readr::col_skip(),
                            readr::col_skip(),
                            readr::col_factor(c("+", "-")),
                            readr::col_factor(c("+", "-")),
                            readr::col_skip(),
                            readr::col_skip(),
                            readr::col_skip(),
                            readr::col_skip(),
                            readr::col_skip(),
                            readr::col_skip(),
                            readr::col_skip(),
                            readr::col_skip(),
                            readr::col_skip(),
                            readr::col_skip(),
                            readr::col_skip(),
                            readr::col_skip(),
                            readr::col_skip(),
                            readr::col_skip()
    )

    for (i in 1:nsamples) {
        coltypes[[1]][[1 + length(coltypes[[1]])]] <- readr::col_integer()
    }

    colnames <- c("Lower.Location",
                  "Lower.Start",
                  "Lower.End",
                  "Upper.Location",
                  "Upper.Start",
                  "Upper.End",
                  "Lower.Strand",
                  "Upper.Strand",
                  paste("Sample", 1:(nsamples), sep = ""))

    data <- readr::read_tsv(filename,
                            col_names = colnames,
                            col_types = coltypes,
                            skip = headskipcount,
                            progress = TRUE)

    as.data.frame(data)
}

# Returns a data frame from a LUMPY bedpe output file
#' @export
load_lumpy <- function(filename) {
    if (!file.exists(filename)) {
        stop("File not found")
    }

    # 0. Preliminary file exam - see how many samples we've got
    GZIPPED <- FALSE

    if (stringi::stri_endswith_fixed(filename, ".gz")) {
        conn <- gzfile(filename)
        GZIPPED <- TRUE
    } else {
        conn <- file(filename)
    }
    open(conn)
    line <- readLines(conn, 1)  # read first line
    headskipcount <- 0
    while (grepl("^##", line)) {
        line <- readLines(conn, 1)
        headskipcount <- headskipcount + 1
    }
    close(conn)

    split <- strsplit(line, "\t")[[1]]
    nFields <- length(split)
    nSamples <- nFields - 21

    ## Use what we've learnt to prepare the readr load_tsv function
    coltypes <- readr::cols(readr::col_character(),  # CHROM_A
                            readr::col_integer(),    # START_A
                            readr::col_integer(),    # END_A
                            readr::col_character(),  # CHROM_B
                            readr::col_integer(),    # START_B
                            readr::col_integer(),    # END_B
                            readr::col_skip(),       # ID
                            readr::col_skip(),       # QUAL
                            readr::col_factor(c("+", "-")), # STRAND_A
                            readr::col_factor(c("+", "-")), # STRAND_B
                            readr::col_factor(c("DUP",
                                                "DEL",
                                                "BND",
                                                "INV")),    # TYPE
                            readr::col_skip(),       # FILTER
                            readr::col_skip(),       # NAME_A
                            readr::col_skip(),       # REF_A
                            readr::col_skip(),       # ALT_A
                            readr::col_skip(),       # NAME_B
                            readr::col_skip(),       # REF_B
                            readr::col_skip(),       # ALT_B
                            readr::col_character(),  # INFO_A
                            readr::col_character(),  # INFO_B
                            readr::col_skip())       # FORMAT


    for (i in 1:nSamples) {
        coltypes[[1]][[1 + length(coltypes[[1]])]] <- readr::col_character()
    }

    remaining <- nFields - length(coltypes[[1]])
    if (remaining > 0) {
        for (i in 1:remaining) {
            coltypes[[1]][[1 + length(coltypes[[1]])]] <- readr::col_skip()
        }
    }

    data <- readr::read_tsv(filename,
                            # col_names = colnames,
                            col_types = coltypes,
                            skip = headskipcount,
                            progress = TRUE)
    data <- as.data.frame(data)
    colnames(data)[1:8] <- c("Lower.Location",
                             "Lower.Start",
                             "Lower.End",
                             "Upper.Location",
                             "Upper.Start",
                             "Upper.End",
                             "Lower.Strand",
                             "Upper.Strand")
    data[, 1:8] <- data[, c(1, 7, 2, 3, 4, 8, 5, 6)]
    colnames(data)[1:8] <- colnames(data)[c(1, 7, 2, 3, 4, 8, 5, 6)]

    # Extract read counts
    for (i in 12:(12+nSamples-1)) {
        data <- cbind(data, extract_lumpy_read_count(data, i))
    }
    data[, -(12:(12+nSamples-1))]
}

# Returns a data frame from a brass bedpe output file.  Filters out commented lines and read ID columns
#' @export
load_brass <- function(filename, header = FALSE) {
    if (!file.exists(filename)) {
        stop("File not found")
    }

    # 0. Preliminary file exam - see how many samples we've got
    GZIPPED <- FALSE

    if (stringi::stri_endswith_fixed(filename, ".gz")) {
        conn <- gzfile(filename)
        GZIPPED <- TRUE
    } else {
        conn <- file(filename)
    }
    open(conn)
    line <- readLines(conn, 1)  # read first line
    headskipcount <- 0
    while (grepl("^#|^%", line)) {
        line <- readLines(conn, 1)
        headskipcount <- headskipcount + 1
    }

    split <- strsplit(line, "\t")[[1]]
    nFields <- length(split)
    if (suppressWarnings(is.na(as.integer(split[nFields])))) {
        nSamples <- (nFields - 9) / 2
    } else {
        nSamples <- nFields - 8
    }

    tailskipcount <- 0

    # gzfile doesn't allow seeking (gzip is a stream compressor,
    # so need to decompress everything prior to seeked-for point)
    if (!GZIPPED) {
        seek(conn, where = -5000, origin = "end")
        line <- readLines(conn, 1)

        while (!identical(line, character(0))) {
            if (grepl("^#|^%", line)) {
                break
            }
            line <- readLines(conn, 1)
        }

        while (!identical(line, character(0))) {
            tailskipcount <- tailskipcount + 1
            line <- readLines(conn, 1)
        }
    }

    close(conn)

    ## Use what we've learnt to prepare the readr load_tsv function
    coltypes <- readr::cols(readr::col_character(),
                            readr::col_factor(c("+", "-")),
                            readr::col_integer(),
                            readr::col_integer(),
                            readr::col_character(),
                            readr::col_factor(c("+", "-")),
                            readr::col_integer(),
                            readr::col_integer())

    for (i in 9:(nSamples + 8)) {
        coltypes[[1]][[i]] <- readr::col_integer()
    }

    remaining <- nFields - length(coltypes[[1]])
    if (remaining > 0) {
        for (i in 1:remaining) {
            coltypes[[1]][[1 + length(coltypes[[1]])]] <- readr::col_skip()
        }
    }

    colnames <- c("Lower.Location",
                  "Lower.Strand",
                  "Lower.Start",
                  "Lower.End",
                  "Upper.Location",
                  "Upper.Strand",
                  "Upper.Start",
                  "Upper.End",
                  paste("Sample", 1:(nSamples), sep = ""))

    data <- readr::read_tsv(filename,
                            col_names = colnames,
                            col_types = coltypes,
                            skip = headskipcount,
                            progress = TRUE)
    res <- as.data.frame(data[1:(nrow(data) - tailskipcount), ])
    res[complete.cases(res),]
}

#' Load a filtered data frame
#' @export
load_filtered <- function(filename) {
    if (!file.exists(filename)) {
        stop("File not found")
    }
    conn <- file(filename)
    open(conn)
    line <- readLines(conn, 1)  # read first line
    close(conn)

    n_cols <- length(strsplit(line, "\t")[[1]])

    coltypes <- readr::cols(readr::col_character(),
                            readr::col_factor(c("+", "-")),
                            readr::col_integer(),
                            readr::col_integer(),
                            readr::col_character(),
                            readr::col_factor(c("+", "-")),
                            readr::col_integer(),
                            readr::col_integer())

    for (i in 9:(n_cols-4)) {
        coltypes[[1]][[i]] <- readr::col_integer()
    }

    coltypes[[1]][[n_cols-3]] <- readr::col_number()
    coltypes[[1]][[n_cols-2]] <- readr::col_number()
    coltypes[[1]][[n_cols-1]] <- readr::col_character()
    coltypes[[1]][[n_cols]] <- readr::col_character()

    res <- as.data.frame(readr::read_tsv(filename, col_types = coltypes))
    res[complete.cases(res),]
}

#' Load a filtered data frame
#' @export
load_filtered <- function(filename) {
    if (!file.exists(filename)) {
        stop("File not found")
    }
    conn <- file(filename)
    open(conn)
    line <- readLines(conn, 1)  # read first line
    close(conn)

    n_cols <- length(strsplit(line, "\t")[[1]])

    coltypes <- readr::cols(readr::col_character(),
                            readr::col_factor(c("+", "-")),
                            readr::col_integer(),
                            readr::col_integer(),
                            readr::col_character(),
                            readr::col_factor(c("+", "-")),
                            readr::col_integer(),
                            readr::col_integer())

    for (i in 9:(n_cols-4)) {
        coltypes[[1]][[i]] <- readr::col_integer()
    }

    coltypes[[1]][[n_cols-3]] <- readr::col_number()
    coltypes[[1]][[n_cols-2]] <- readr::col_number()
    coltypes[[1]][[n_cols-1]] <- readr::col_character()
    coltypes[[1]][[n_cols]] <- readr::col_character()

    as.data.frame(readr::read_tsv(filename, col_types = coltypes))
}
