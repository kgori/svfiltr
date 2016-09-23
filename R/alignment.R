#' Convert data to binary
binary_phylip <- function(data, threshold = 5) {
    lines <- list()
    lines[1] <- paste(ncol(data), nrow(data))
    max_name_len <- max(nchar(colnames(data)))
    for (name in colnames(data)) {
        characters <- paste0(as.numeric(data[[name]] > threshold),
                             collapse = "")
        line <- stringr::str_c(name,
                               stringi::stri_dup(" ", max(10, max_name_len + 1) - nchar(name)),
                               characters)
        lines[length(lines) + 1] <- line
    }
    lines
}

#' Write data as binary sequence alignment in Phylip format
#' @export
write_phylip <- function(filename, data, threshold = 5) {
    phylip_lines <- binary_phylip(data, threshold)
    conn <- file(filename, "w")
    writeLines(as.character(phylip_lines), conn)
    close(conn)
}
