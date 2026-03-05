# calculates the wave speed after simulations are done

wave_speed <- function(solocs.all){
    # columns are v, l, r, time, x, y, "unknown" (cell quality), S, E, I, R, C, Z

    # get introduction point
#     setnames(solocs.all, c('ctX', 'ctY'), c('x', 'y'))
    setorder(solocs.all, v, l, r, time, x, y)
    init.pts <- solocs.all[E > 0 | I > 0 | C > 0, .SD[1], by=.(v, l, r)]
    setnames(init.pts, c('x', 'y'), c('xi', 'yi'))

    # get incidence by cell by time
    # time during infection should be contiguous b/c simulations end if there are no E, I, or C in the population
    solocs.all <- solocs.all[E > 0 | I > 0 | C > 0, ]
    solocs.all[, tot.inf := E + I + C]

    # grab the distance of each cell from the origin (could be cached for speed)
    # connect initial points to infected points
    solocs.all <- unique(solocs.all[, .(v, l, r, time, tot.inf, x, y)][init.pts[, .(v, l, r, xi, yi)], on=.(v, l, r)])
    solocs.all[, dist := sqrt((x - xi)^2 + (y - yi)^2)]

    # get mean, median, sd of infected cell distances
    solocs.all[, avg := mean(dist), by=.(v, l, r, time)]
    solocs.all[, var := var(dist),  by=.(v, l, r, time)]
    solocs.all[, max := max(dist),  by=.(v, l, r, time)]
    solocs.all[, min := min(dist),  by=.(v, l, r, time)]

    # get speed of average distance
    avg.dist.diff <- unique(solocs.all[,.(v, l, r, time, avg)])
    avg.dist.diff[, avg.dist.diff := avg - shift(avg, 1), by=.(v, l, r)]
    # na values in shifted average distances at time 1 are from not having time 0 data to shift from, so 0 makes sense
    avg.dist.diff[is.na(avg.dist.diff) & time == 1, avg.dist.diff := 0]

    # connect average wave speed to output data table
    solocs.all <- solocs.all[avg.dist.diff, on=.NATURAL]

    return(solocs.all)
}
