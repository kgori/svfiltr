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
                            progress = T)
    as.data.frame(data)
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
