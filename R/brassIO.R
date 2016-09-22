# Returns a data frame from a brass bedpe output file.  Filters out commented lines and read ID columns
#' @export
load_brass <- function(filename, header = FALSE) {
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
