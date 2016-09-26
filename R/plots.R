#' @export
#' @importFrom "grDevices" rgb
#' @importFrom "circlize" circos.par circos.initialize circos.trackPlotRegion
#' @importFrom "circlize" circos.clear circos.genomicLink circos.axis get.cell.meta.data
circos_plot <- function(df){

    # 1. Define borders
    chromosome.ranges <- matrix(0, ncol=2, nrow=40)
    chromosome.ranges[,1] <- 1 # define all chromosomal start positions (= 1)
    chromosome.ranges[,2] <- c(122678785, 85426708, 91889043, 88276631, 88915250, 77573801, 80974532, 74330416,
                               61074082, 69331447, 74389097, 72498081, 63241923, 60966679, 64190966, 59632846,
                               64289059, 55844845, 53741614, 58134056, 50858623, 61439934, 52294480, 47698779,
                               51628933, 38964690, 45876710, 41182112, 41845238, 40214260, 39895921, 38810281,
                               31377067, 42124431, 26524999, 30810995, 30902991, 23914537, 123869142, 16727) # define all chromosomal lengths

    rownames(chromosome.ranges) <- c(1, 2, 3, 4, 5, 6, 7, 8,
                                     9, 10, 11, 12, 13, 14, 15, 16,
                                     17, 18, 19, 20, 21, 22, 23, 24,
                                     25, 26, 27, 28, 29, 30, 31, 32,
                                     33, 34, 35, 36, 37, 38, "X", "MT")

    # 2. Plot
    #color.ramps <- 1-c(1/df[,grep('READS', colnames(df))[1]]) # darker colour for SVs with more read supports (could still think of a nicer function...)
    circos.par("track.height" = 0.2, cell.padding = c(0, 0, 0, 0))

    # Initialise
    circos.initialize(factors = rownames(chromosome.ranges),
                      xlim = chromosome.ranges)

    # Build circos-'frame'
    circos.trackPlotRegion(ylim = c(0, 1),
                           panel.fun = function(x, y) {print(get.cell.meta.data("xlim"))},
                           track.height = 0.02,
                           bg.col = c(2:c(nrow(chromosome.ranges) + 1)),
                           bg.border = c(2:c(nrow(chromosome.ranges) + 1)),
                           track.index = 1)

    # Add chromosome names
    for (i in 1:nrow(chromosome.ranges)){
        circos.axis(h='top',sector.index = rownames(chromosome.ranges)[i],
                    major.at = chromosome.ranges[i,2]/2,
                    labels = rownames(chromosome.ranges)[i],
                    direction = "outside", major.tick.percentage = 1, labels.cex=3,
                    labels.away.percentage=1/1.2, minor.ticks = 4)
    }

    # Add links
    reg1 <- df[,c(1,3,4)]
    reg2 <- df[,c(5,7,8)]
    stopifnot(nrow(reg1) == nrow(reg2))
    rownames(reg1) <- seq(1, nrow(reg1))
    rownames(reg2) <- seq(1, nrow(reg2))

    circos.genomicLink(region1=reg1,region2=reg2,
                       col=rgb(1,0,0,1.0), #color.ramps*0.3),
                       lwd=5, rou=0.9)
    circos.clear()
}


#' @importFrom "graphics" abline par text
#' @export
plot_hists <- function(df, index, mfrow = c(3, 4), boundaries = NULL, ...) {
    names <- colnames(df)[index]
    par(mfrow = mfrow)
    for (i in 1:length(names)) {
        dat <- df[[names[i]]]

        # Ignore very high numbers of reads
        dat <- dat[dat < 1000]

        MASS::truehist(dat, main = names[i], xlab = "reads", ...)

        if (!is.null(boundaries)) {
            boundary <- boundaries[[names[i]]]

            # Get size of current graphics to work out text offsets
            curr_ymax <- par("usr")[4]
            curr_xrange <- par("usr")[2] - par("usr")[1]

            # Plot a vertical line at boundary, and add text label
            abline(v = boundary, col = "red")
            text(x = boundary + 0.1 * curr_xrange,
                 y = curr_ymax * 0.9,
                 col = "red",
                 labels = paste("=", boundary, sep = ""))
        }
    }
}
