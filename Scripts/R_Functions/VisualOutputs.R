
VisualOutputs <- function(out.list, variables, land_grid_list, parameters, lands_names){
    ## NEED TO:
    ## -- add land image to background of gifs?
    ## -- get ldsel table from NND calc for attributes of landscapes

    # Generates output plots for examples
    library(data.table)

    setDT(variables)

    # put outputs in variables
    lapply(names(out.list), function(x) eval(parse(text=paste0(x, '= as.data.table(out.list["',x,'"])'))))
    lapply(names(out.list), function(x) paste0(x, '= as.data.table(out.list["',x,'"])'))
    lapply(out.list, function(x) print(names(x)))

    # make things into data.tables and name columns for sanity
    tm.mat <- as.data.table(out.list["tm.mat"])
    setnames(tm.mat, unlist(lapply(strsplit(names(tm.mat), 'tm.mat.'), function(x) unlist(x)[2])))
    tm.mat <- tm.mat[!is.na(timestep),]
    summ.vals <- as.data.table(out.list["summ.vals"])
    setnames(summ.vals, unlist(lapply(strsplit(names(summ.vals), 'summ.vals.'), function(x) unlist(x)[2])))
    summ.vals <- summ.vals[!is.na(endtime)]
    incidence <- as.data.table(out.list["incidence"])
    setnames(incidence, unlist(lapply(strsplit(names(incidence), 'incidence.'), function(x) unlist(x)[2])))
    incidence[, max.time := max(timestep), by=.(var, rep, land)]
    detections <- as.data.table(out.list["detections"])
    setnames(detections, unlist(lapply(strsplit(names(detections), 'detections.'), function(x) unlist(x)[2])))
    detections[, max.time := 0]
    if(nrow(detections) > 0) {
        detections[,max.time := max(time), by=.(var, rep, land)]
        unq.det <- unique(detections[, .(var, land, rep, timestep, loc, max.time, code, detected)])
    }
    allzones <- as.data.table(out.list["allzone"])
    setnames(allzones, unlist(lapply(strsplit(names(allzones), 'allzone.'), function(x) unlist(x)[2])))
    solocs.all <- as.data.table(out.list['solocs.all'])
    setnames(solocs.all, unlist(lapply(strsplit(names(solocs.all), 'solocs.all.'), function(x) unlist(x)[2])))
#     setnames(solocs.all, c('x','y'), c('ctX','ctY'))
    solocs.all[,nlive := S+E+I+R]

    wv.speed <- wave_speed(solocs.all)

    # generate an image for environment quality grid
    ## only works for single grid in land_grid_list
    grid.centers.out <- rbindlist(lapply(seq(length(land_grid_list)), function(i){
        grid.key <- as.data.table(land_grid_list[[i]][[2]])
        setnames(grid.key, c('cell','tlX','tlY','trX','trY','ctX','ctY','rasval')[1:ncol(grid.key)]) ## see Make_Grid.R
        grid.centers <- cbind(grid.key[,.(cell, ctX, ctY, rasval)], land=i)
        return(grid.centers)
    }))

    ## get some plots of temporal data
    # time-based matrix outputs
    tm.pop.unq <- unique(tm.mat[,.(var, land)])
    rows.plt <- round(sqrt(length(unique(tm.pop.unq[,var]))))
    cols.plt <- ceiling(sqrt(length(unique(tm.pop.unq[,var]))))
    for (l in unique(tm.pop.unq[,land])){
        png(paste0('./test_outputs/tm.plots_land_', l, '.png'), width=450*cols.plt, height=500*rows.plt)
        ## make plot dimensions dynamics based on parameter inputs
        pdef <- par(mfrow=c(rows.plt, cols.plt),
                oma=c(8, 1.2, 0.2, 0.2),
                lwd=2.3, cex.lab=1.8, cex.axis=1.6, cex.main=2, cex.sub=1.4)
        # par(mfrow=c(1,1), oma=c(8,0.2,0.2,0.2),lwd=2.3, cex.lab=1.8, cex.axis=1.6, cex.main=2, cex.sub=1.4)
        # par(mfcol=c(length(unique(tm.mat[,var])), length(unique(tm.mat[,land]))))
        # set plot ranges based on maximum of everything that will be on the multiplot figure
        xrng = range(tm.mat[,timestep])
        yrng = log1p(range(tm.mat[,.(BB, S, E, I, R, C, Z)]))
        # loop over the parameter combinations for each plot panel
        seirczbb.temporal <- function(v, l, plt.i, rows.plt, cols.plt, dat=tm.mat, vlist = variables){
            vardat <- paste(vlist[v,], collapse=' ') # quick and dirty parameter inclusion
            # xlabel default
            xlabi <- ''
            # ylabel default
            ylabi <- ''
            # default plot panel margin sizes
            bmar <- lmar <- 2
            if (plt.i / rows.plt >= rows.plt){ # first row of panels
                xlabi <- 'time (weeks)'
                bmar <- 4
            }
            if (plt.i %% cols.plt == 1){ # first column of panels
                ylabi <- 'log(1+individuals)'
                lmar <- 4
            }
            pdef2 <- par(mar = c(bmar, lmar, 1, 1))
            # create a blank plot
            plot(0,0,xlim=xrng, ylim=yrng, col=NULL, ann=FALSE)
            mtext(xlabi, 1, line=2.5, cex=1.8)
            mtext(ylabi, 2, line=2.5, cex=1.8)
            mtext(paste('vars',v,'| land',l), 3, line=-1.5, cex=1.7, font=2)
            mtext(vardat, 3, line=-2.7)
            # subset with things that have only the land tile and variable combination
            subdat <- dat[var==v & land==l,]
            # count up reps
            reps <- unique(subdat[,rep])
            # draw a line of each type for each rep
            lapply(reps, function(r) lines(log1p(S) ~ timestep, data=subdat[rep==r | rep == 0,], col='olivedrab', lty=1))
            lapply(reps, function(r) lines(log1p(E) ~ timestep, data=subdat[rep==r | rep == 0,], col='orange', lty=1))
            lapply(reps, function(r) lines(log1p(I) ~ timestep, data=subdat[rep==r | rep == 0,], col='red', lty=1))
            lapply(reps, function(r) lines(log1p(R) ~ timestep, data=subdat[rep==r | rep == 0,], col='blue', lty=1))
            lapply(reps, function(r) lines(log1p(C) ~ timestep, data=subdat[rep==r | rep == 0,], col='purple', lty=1))
            lapply(reps, function(r) lines(log1p(Z) ~ timestep, data=subdat[rep==r | rep == 0,], col='black', lty=1))
            lapply(reps, function(r) lines(log1p(BB) ~ timestep, data=subdat[rep==r | rep == 0,], col='pink', lty=1))
            lapply(reps, function(r) abline(v=subdat[rep==r,max(timestep)], col='grey'))
            lapply(reps, function(r) abline(v=subdat[rep==0,max(timestep)], col='black', lty=3))
            par(pdef2)
        }
        mapply(seirczbb.temporal, v=tm.pop.unq[land==l, var], l=tm.pop.unq[land==l,land], plt.i=seq(nrow(tm.pop.unq[land==l,])), MoreArgs = list(rows.plt=rows.plt, cols.plt=cols.plt))
        # create a legend under the panels
        par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), mar=c(0, 0, 0, 0), new=TRUE)
        plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
        legend("bottom", legend=c('S','E','I','R','C','Z','Births'),#'rep1','rep2'),
            lty=c(rep(1,8), 2),
            col=c('olivedrab','orange','red','blue','purple','black','pink'),#'grey','grey'),
            horiz=TRUE, bty='n', cex=2.3, lwd=3)
        dev.off()
        par(pdef)
    }


    ## plots of the spatial data
    ### spatial habitat quality
#     lapply(seq(length(land_grid_list)), function(lnd){
#         png(paste0('./test_outputs/landscape_',lnd,'.png'))
#         rasvalmat <- as.matrix(dcast(grid.centers.out[land==lnd,], ctX ~ ctY, value.var = c('rasval')))
#         plot(range(grid.centers.out[land==lnd,ctX]), range(grid.centers.out[land==lnd,ctY]), type='n')
#         ## will need to set cex for points to have even squares based on output plot size
#         points(ctY ~ ctX, data=grid.centers.out[land==lnd,], col=terrain.colors(5, rev=TRUE)[grid.centers.out[land==lnd,rasval*10]], pch=15)
#         dev.off()
#     })
#     browser()

    ### timing and extent of incidence/detections
    # connect incidence with cell locations
    unq.eic <- unique(incidence[,.(var,land,rep,timestep,loc,max.time)])
#     mapply(function(v, l){
#         input.dat <- unq.eic[var==v & land==l & max.time >= 10,]
#         if(nrow(input.dat) != 0) {
# #         dev.new()
# #     browser()
#         usable.reps <- unique(input.dat[,rep])
#         png(paste0('./test_outputs/incidence_v', v, 'l', l, '.png'), width=1000, height=500*length(usable.reps))
#     #     if(nrow(input.dat) == 0) {browser(); input.dat <- unq.eic[var==v & land==l,]}
#             defpar <- par(mfrow=c(length(unique(input.dat[,rep])),2), oma=c(8,1.2,0.2,0.2))
#
#             lapply(usable.reps, function(x){
#                 inc.cell <- input.dat[rep==x,][grid.centers, on=.(loc = cell)][!is.na(timestep),]
#                 setorder(inc.cell, timestep)
#                 inc.cell.early <- inc.cell[, .SD[1], by=.(var, land, rep, loc)]
#                 plot(range(grid.centers[,ctX]), range(grid.centers[,ctY]), type='n', main=paste('incidence: vars',v,'| land',l), cex.lab = 1.4, xlab='x coordinate (Km)', ylab='y coordinate (Km)')
#                 points(ctY ~ ctX, data=inc.cell.early, col=rainbow(max(timestep))[timestep], pch=16,cex=0.4)
#                 points(ctY ~ ctX, data=inc.cell.early[timestep==1,], pch='X', cex=2)
# #                 browser()
#                 hist(inc.cell[,timestep], breaks=max(inc.cell[,timestep]), col=rainbow(max(inc.cell[,timestep])), main='Timing of Incidence', ylab='# New Incidence', xlab='Time (weeks)', cex.lab = 1.4)
#                 # Add line to indicate when first detection began
#                 abline(v=parameters['detectday'], lwd=3, lty=2, col='black')
#             })
#             # add legend
#             # create a legend under the panels
#             par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), mar=c(0, 0, 0, 0), new=TRUE)
#             plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
#             legend("bottom", legend=c('Intro point','First detect'),#'rep1','rep2'),
#                 lty=c(NA, 2),
#                 col=c('black','black'),#'grey','grey'),
#                 pch=c('X',NA),
#                 horiz=TRUE, bty='n', cex=2.3, lwd=3)
#             par(defpar)
#             dev.off()
#         }
#         },
#     v=tm.pop.unq[,var],
#     l=tm.pop.unq[,land]
#     )

    # connect detections with cell locations
#     unq.eic <- unique(detections[!is.na(loc) & detected != 0 & max.time >= 10,.(var,land,rep,timestep,loc,max.time,code,detected)])
#     mapply(function(v, l){
#         input.dat <- unq.eic[var==v & land==l & max.time >= 10,]
#         if(nrow(input.dat) != 0) {
#             input.dat <- aggregate(input.dat[,detected], list(rep=input.dat[,rep], timestep=input.dat[,timestep], loc=input.dat[,loc]), sum)
#             setDT(input.dat)
#             setnames(input.dat, 'x', 'detected')
#             setorderv(input.dat, cols = c('rep','timestep','loc'))
#         #         dev.new()
#             usable.reps <- unique(input.dat[,rep])
#             png(paste0('./test_outputs/detection_v', v, 'l', l, '.png'), width=1000, height=500*length(usable.reps))
#             defpar <- par(mfrow=c(length(unique(input.dat[,rep])),2))
#
#             lapply(unique(input.dat[,rep]), function(x){
#                 inc.cell <- input.dat[rep==x,][grid.centers, on=.(loc = cell)][!is.na(timestep),]
#                 setorder(inc.cell, timestep)
#     #             inc.cell.early <- inc.cell[, .SD[1], by=.(var, land, rep, loc)]
# #                 browser()
#                 plot(range(grid.centers[,ctX]), range(grid.centers[,ctY]), type='n', main=paste('detections: vars',v,'| land',l), cex.lab = 1.4, xlab='x coordinate (Km)', ylab='y coordinate (Km)')
#                 points(ctY ~ ctX, data=inc.cell, col=rainbow(max(timestep))[timestep], pch=16,cex=1.4)
#                 points(ctY ~ ctX, data=inc.cell[timestep==min(inc.cell[,timestep]),], pch='X', cex=2)
#                 # setup for a stacked barplot
#                 bp <- dcast(unq.eic[var==v & land==l & max.time >= 10 & rep==x,][grid.centers, on=.(loc = cell), nomatch=NULL], timestep ~ code, value.var='detected', fun.aggregate=sum)[!is.na(timestep),]
# #                 setnames(bp, c('0','1'), c('C','IE'))
# #                 bp[,tot := C + IE]
#                 col.raw = rainbow(max(bp[,timestep]))[bp[,timestep]]
#                 if (nrow(as.matrix(bp)) > 1){
#                     barplot(t(as.matrix(bp))[c(2,3),], space=0, xlab='timestep', col='white', ylab='detected (dark=dead)', main = x, names.arg=t(as.matrix(bp))[1,])
#                     # dark portion of bars indicates how many detected were found dead vs. culled
#                     # makes two-tone stacked barplots (https://stackoverflow.com/a/59411350)
#                     # uses darken() from 'colorspace' package; probably is a non-package way to shift color hex codes
#                     for (i in 1:ncol(t(as.matrix(bp)))){
#                         xx = t(as.matrix(bp))[2:3,]
#                         xx[,-i] = NA
#                         colnames(xx)[-i] = NA
#                         barplot(xx,col=c(darken(col.raw[i],0.4),col.raw[i]), add=T, axes=F, space=0)
#                     }
#                 }
#
# #                 plot(detected ~ timestep,data=inc.cell, main='Timing of Detections', ylab='# New Detections', xlab='Time (weeks)', cex.lab = 1.4, type='b', pch=16)
#             })
#             par(defpar)
#             dev.off()
#         }
#     },
#     v=tm.pop.unq[,var],
#     l=tm.pop.unq[,land]
#     )

    radius = parameters['Rad']
#     unq.combos <- unique(unq.eic[max.time > 10,.(var, land, rep)])
    if (radius != 0 & parameters[['detectday']] < parameters[['thyme']]){
        ## plotting incidence with zone of control over time for a bunch of plots
        # detection locations
        unq.detections <- unique(detections[,.(var,land,rep,time,loc,max.time,detected)])

        unq.incidence <- unique(incidence[,.(var,land,rep,timestep,loc,max.time)])
        unq.incidence[,is.inf := 1]
        unq.incidence[,loc.min := min(timestep), by=.(var, land, rep, loc)]
        unq.incidence[,loc.max := max(timestep), by=.(var, land, rep, loc)]
        unq.join <- CJ(timestep=seq(parameters[['thyme']]),
                    var=unique(unq.incidence[,var]),
                    land=unique(unq.incidence[,land]),
                    rep=unique(unq.incidence[,rep]),
                    loc=unique(unq.incidence[,loc]))
        unq.join <- unq.join[!unq.incidence, on=.(var, land, rep, timestep, loc)]
        unq.join[, is.inf := 0]
        unq.join <- unq.join[unique(unq.incidence[,.(var, land, rep, loc, loc.min, loc.max, max.time)]), on=.(var, land, rep, loc)]
        unq.join <- unq.join[timestep >= loc.min,]
    #         unq.join <- unq.join[timestep <= loc.max & timestep >= loc.min,]
        unq.incidence <- rbind(unq.incidence, unq.join)

        unq.incidence[,loc.min := NULL]
        unq.incidence[,loc.max := NULL]
        setorder(unq.incidence, var, land, rep, timestep, loc)
        ## gets all runs than last at least 50 weeks, infects at least 5 cells during the course of the run, and then takes the first and last rep of each var/land combo to use for gif example
        unq.combos <- unique(unq.incidence[max.time > 50,sum(is.inf), by=.(var, land, rep, timestep)][,max(V1), by=.(var, land, rep)][V1 > 5, .(var,land,rep)])[, .SD[.N], by=.(var, land)]

    #         unq.combos <- unq.combos[var == 7 & rep == 2,]
#         mapply(function(i, j, k){
#             sub.incidence <- unique(unq.incidence[var == i & land == j & rep == k & loc != 0,.(timestep, loc, is.inf)])[grid.centers.out[land==j,], on=.(loc = cell), nomatch=NULL] # gets rid of code column
#             if(nrow(unq.detections) > 0){ sub.detections <- unique(unq.detections[var == i & land == j & rep == k & detected != 0,.(time, loc)])[grid.centers.out[land==j,], on=.(loc = cell), nomatch=NULL] } # gets rid of code column
#             if(nrow(allzones) >0){sub.zones <- allzones[var == i & land == j & rep == k & loc != 0,][grid.centers.out[land==j,], on=.(loc=cell), nomatch=NULL]}
#             sub.solocs <- solocs.all[v == i & l == j & r == k,]
#             paneldim <- ceiling(sqrt(max(sub.incidence[,timestep])))
#             png(paste0('./test_outputs/sptmplot_vars_', i,'_land_', j, '_rep_', k,'.png'), width=2000, height=2000)
#             par(mfrow = c(paneldim, paneldim), oma=c(0,0,0,0), mar=c(0,0,0,0))
#             lapply(seq(max(sub.incidence[,timestep])), function(x){
#                 plot(y ~ x, data=sub.solocs[time == x,], pch='.', col='gray', xlim=c(0,100), ylim=c(0,100))
#                 if(nrow(allzones) > 0){points(ctY ~ ctX, data=sub.zones[timestep <= x,], pch='.', col='yellow')}#, cex=1.2, xlim=c(0,100), ylim=c(0,100))
#     #                 points(y ~ x, data=sub.solocs[timestep == x,], pch=1, col='gray')
#     #                 points(ctY ~ ctX, data=sub.incidence[timestep <= x & is.inf == 0,], pch='.', col='red')
#                 points(ctY ~ ctX, data=sub.incidence[timestep == x & is.inf == 1,], pch=3, col='red')
#                 if (x > parameters['detectday'] & exists('sub.detections')) points(ctY ~ ctX, data=sub.detections[time <= x,], pch=2, col='blue')
#             })
#             dev.off()

#             library(animation)
#             plot.step <- function(x){
#                 panels <- layout(matrix(c(1, 1, 1, 2, 3, 4), nrow=3, ncol=2), widths=c(3, 1), heights=c(1, 1, 1), respect=FALSE)
#                 plot(y ~ x, data=sub.solocs[time == x,], pch='.', col='gray', xlim = c(0, 100), ylim = c(0, 100), main=paste('week', x), cex=log1p(nlive))
#                 if(nrow(allzones) > 0) points(ctY ~ ctX, data=sub.zones[time <= x,], pch=3, col='yellow')
#                 points(ctY ~ ctX, data=sub.incidence[timestep == x & is.inf == 1,], pch=3, col='red')
#                 if (x > parameters['detectday'] & exists('sub.detections')) points(ctY ~ ctX, data=sub.detections[time <= x,], pch=2, col='blue')
#                 plot(allnlive ~ time, data=sub.solocs[, allnlive := sum(nlive), by=time][order(time)], main='Live pigs', type='l')
#                 points(allnlive ~ time, data=sub.solocs[time==x,], pch=1, cex=2)
#                 plot(allnlive ~ timestep, data=unique(sub.incidence[, allnlive := sum(is.inf), by=timestep][,.(timestep, allnlive)])[order(timestep)], main='Infected cells', type='l', col='red')
#                 points(allnlive ~ timestep, data=sub.incidence[timestep==x,], pch=3, col=2, cex=2)
#                 if (exists('sub.detections')){
#                     temp.detect <- unique(sub.detections[, allnlive := length(unique(loc)), by=time][,.(time, allnlive)])[order(time)][,cs:=cumsum(allnlive)]
#                     plot(cs ~ time, data=temp.detect, main='Detections', type='l', col='blue')
#                     points(cs ~ time, data=temp.detect[time==x, ], pch=2, col='blue', cex=2)
#                 }
#             }
#             out.gif <- paste0('testgif_vars_', i, '_land_', j, '_rep_', k, '.gif')
#             saveGIF({
#                 lapply(seq(min(sub.incidence[,timestep]), max(sub.incidence[,timestep])), plot.step)
#             }, movie.name = out.gif, ani.width=850, ani.height=600, interval=0.1, imgdir='./test_outputs')
#         }, i=unq.combos[,var], j=unq.combos[,land], k=unq.combos[,rep])

        # moves gifs to the proper folder
#         lapply(list.files(pattern='*.gif'), function(x) file.rename(x, paste0('./test_outputs/',x)))
    }



    ## incidence/time -- averaged over replicates (state in this is SEI, not FL/SC)
#     incidence <- edge.interactions[incidence, on=.(v=var, l=land, r=rep, time=timestep)]
#     setnames(incidence, c('v','l','r','time'), c('var','land','rep','timestep'))
# #     incidence[edge==0,edge.time:= min(timestep), by=.(var, land, rep)]
#     incidence <- incidence[is.na(edge),]#[,edge:=NULL]
#     inc.time <- dcast(incidence[,.N, by=.(var, land, rep, state, timestep)], timestep ~ var + land + state, value.var = 'N', fun.aggregate = mean, na.rm=TRUE, fill=NA)
#     # incidence by addition/subtraction by timestep (i.e. rate of change by week)
#     inc.time.rate <- inc.time[,lapply(.SD, function(x) x - data.table::shift(x, 1)), .SDcols=2:ncol(inc.time)]
#
#     combonames <- unlist(lapply(strsplit(names(inc.time.rate), '_[a-z]'), function(x) x[[1]]))
#     name.index <- seq(3)
#     pltdim <- unlist(unique(incidence[,.(var, land)])[,lapply(.SD, function(x) length(unique(x)))])
#     pltdim[1] <- ceiling(pltdim[1]/pltdim[2])
#
#     mapply(function(v, l){
#         if(v==1 & l != 1) {
#             dev.off()
#             png(paste0('./test_outputs/inc.time_land_',l,'.png'), width=pltdim[1]*450, height = pltdim[2]*500)
#             par(mfrow=pltdim, mar=c(2,2,0,0), oma=c(2,2,0.2,0.2))
#         } else if (v==1 & l == 1){
#             png(paste0('./test_outputs/inc.time_land_',l,'.png'), width=pltdim[1]*450, height = pltdim[2]*500)
#             par(mfrow=pltdim, mar=c(2,2,0,0), oma=c(2,2,0.2,0.2))
#         }
#         in.rate <- inc.time.rate[,.SD, .SDcols=patterns(paste0('^', v, '_', l))]
#         if(length(names(in.rate)) < 3){
#             present.names <- unlist(lapply(strsplit(names(in.rate), paste0(v, '_', l, '_')), function(x) x[[2]]))
#             name.index <- which(c('carcass','exposed','infected') %in% present.names)
#         }
#         vnames <- unlist(variables[v, .(state, variant, density)])
#         vnames <- paste('st:',vnames[1], 'strn:',vnames[2], 'dens:', vnames[3])
#         matplot(in.rate, type='l', lty=1, col=name.index, main=vnames)
#         abline(h=0, col='grey')
#     }, v=unique(incidence[,.(var, land)])[,var]
#     , l=unique(incidence[,.(var, land)])[,land])
#     dev.off()
#
#     ## incidence/time -- averaged over replicates but by cumulative incidence
#     inc.time.cum <- inc.time[,lapply(.SD, cumsum), .SDcols=2:ncol(inc.time)]
#     mapply(function(v, l){
#         if(v==1 & l != 1) {
#             dev.off()
#             png(paste0('./test_outputs/inc.cum_land_',l,'.png'), width=pltdim[1]*450, height = pltdim[2]*500)
#             par(mfrow=pltdim, mar=c(2,2,0,0), oma=c(2,2,0.2,0.2))
#         } else if (v==1 & l == 1){
#             png(paste0('./test_outputs/inc.cum_land_',l,'.png'), width=pltdim[1]*450, height = pltdim[2]*500)
#             par(mfrow=pltdim, mar=c(2,2,0,0), oma=c(2,2,0.2,0.2))
#         }
#         in.rate <- inc.time.cum[,.SD, .SDcols=patterns(paste0('^', v, '_', l))]
#         if(length(names(in.rate)) < 3){
#             present.names <- unlist(lapply(strsplit(names(in.rate), paste0(v, '_', l, '_')), function(x) x[[2]]))
#             name.index <- which(c('carcass','exposed','infected') %in% present.names)
#         }
#         vnames <- unlist(variables[v, .(state, variant, density)])
#         vnames <- paste('st:',vnames[1], 'strn:',vnames[2], 'dens:', vnames[3])
#         matplot(in.rate, type='l', lty=1, col=name.index, main=vnames)
#         abline(h=0, col='grey')
#     }, v=unique(incidence[,.(var, land)])[,var]
#     , l=unique(incidence[,.(var, land)])[,land])
#     dev.off()

    ## maximum incidence rate (boxplots)


#     tm.mat[,EI := E+I]
#     tm.mat[,dEI := EI-data.table::shift(EI,1), by=.(var,land,rep)]
#     tm.mat[is.na(dEI) & timestep==1,dEI := 0]
#     tm.mat[,incidenceEI := dEI/S]
#     max.inc <- tm.mat[,max(incidenceEI), by=.(var,land,rep)]
#     matplot(exp(dcast(tm.mat, timestep ~ var + land + rep, value.var = 'incidenceEI', fun.aggregate=mean, na.rm=TRUE, fill=NA)[,-c(1)]), type='l', lty=1)

    edge.interactions <- unique(wv.speed[x==0.25 | y==0.25 | x==99.75 | y==99.75, .(v,l,r,time)])[,edge := 0]
    tm.mat.edge <- edge.interactions[tm.mat, on=.(v=var, l=land, r=rep, time=timestep)]
    setnames(tm.mat.edge, c('v','l','r','time'), c('var','land','rep','timestep'))
    tm.mat.edge[is.na(edge), seq.test := 1:.N, by=.(var,land,rep)]
    tm.mat.edge <- tm.mat.edge[timestep == seq.test,][,seq.test := NULL]

#     tm.mat.edge[,EIC := E+I+C]
#     tm.mat.edge[,dEIC := EIC-data.table::shift(EIC,1), by=.(var,land,rep)]
#     tm.mat.edge[is.na(dEIC) & timestep==1,dEIC := 0]
#     tm.mat.edge[,incidence := dEIC/S]
#     max.inc <- tm.mat.edge[,max(incidence), by=.(var,land,rep)]
#     matplot((dcast(tm.mat.edge, timestep ~ var + land + rep, value.var = 'incidence', fun.aggregate=mean, na.rm=TRUE, fill=NA)[,-c(1)]), type='l', lty=1)

    tm.mat.edge[,EI := E+I]
    tm.mat.edge[,dEI := EI-data.table::shift(EI,1), by=.(var,land,rep)]
    tm.mat.edge[is.na(dEI) & timestep==1,dEI := 0]
    tm.mat.edge[,incidenceEI := dEI/S]
#     matplot((dcast(tm.mat.edge, timestep ~ var + land + rep, value.var = 'incidenceEI', fun.aggregate=mean, na.rm=TRUE, fill=NA)[,-c(1)]), type='l', lty=1, main="Incidence", ylab="(E+I_t - E+I_t-1)/S_t", xlab='Week')

    max.inc <- tm.mat.edge[,max(incidenceEI), by=.(var,land,rep)]
    max.inc <- max.inc[variables[,id := 1:.N], on=.(var=id)]
#     boxplot(V1 ~ contact + variant + as.factor(land) + as.factor(density), data=max.inc, ylab='maximum incidence')

    ## wavespeed stuff
    ## 4 plots for temporal dynamics -- incidence/time, max distance/time, peak speed, area/time
    png('./test_outputs/wavespeed.png', width=1200, height=1400)
    par.og <- par(mfrow=c(2,2), oma=c(8,1.2,0.2,0.2), mar=c(4,4,1,1))

    tm.mat.inc <- dcast(tm.mat.edge, timestep ~ var + land, value.var = 'incidenceEI', fun.aggregate=mean, na.rm=TRUE, fill=NA)[,-c(1)]
    matplot(tm.mat.inc, type='b', lty=1, main="Incidence", ylab="(E+I_t - E+I_t-1)/S_t", xlab='Week', col=unlist(tstrsplit(names(tm.mat.inc), '_', keep=1)), pch=as.numeric(unlist(tstrsplit(names(tm.mat.inc), '_', keep=2))))

    wv.speed <- edge.interactions[wv.speed, on=.NATURAL]
    wv.speed <- wv.speed[is.na(edge),]
#     avg.wave <- dcast(wv.speed, time~ v+l, value.var='avg', fun.aggregate=mean, fill=NA)
#     matplot(avg.wave[,2:ncol(avg.wave)], type='b', lty=1, col=rep(1:nrow(variables), each=length(land_grid_list)), pch=seq(length(land_grid_list)),
#             ylab='avg dist from source', xlab='week', main='Avg dist from introduction')

    max.dist.wave <- dcast(wv.speed, time~ v+l, value.var='max', fun.aggregate=max, fill=NA)
    yrng <- range(max.dist.wave[,2:ncol(max.dist.wave)], na.rm=TRUE)
    matplot(max.dist.wave[,2:ncol(max.dist.wave)], type='b', lty=1, col=rep(1:nrow(variables), each=length(land_grid_list)), pch=seq(length(land_grid_list)),
            ylab='max dist from source', xlab='week', main='Max dist from introduction', ylim=yrng)
#     matplot(max.dist.wave.edge[,2:ncol(max.dist.wave.edge)], type='l', lty=2, col=rep(1:nrow(variables), each=length(land_grid_list)), pch=seq(length(land_grid_list)), add=TRUE)

    spd.wave <- dcast(unique(wv.speed[is.na(edge),.(v,l,r,time,avg.dist.diff)]), time~ v+l, value.var='avg.dist.diff', fun.aggregate=mean, fill=NA)
    matplot(spd.wave[,2:ncol(spd.wave)], type='b', lty=1, col=rep(1:nrow(variables), each=length(land_grid_list)), pch=seq(length(land_grid_list)),
            ylab='km/week', xlab='week', main='Avg. wave speed from introduction point')
    abline(h=0, col='grey')

#     pk.dist <- dcast(wv.speed[order(v,l,r,time,-tot.inf)][,.SD[1], by=.(v,l,r,time)], time~ v+l, value.var=c('tot.inf','dist'), fun.aggregate=mean, fill=NA)
#     matplot(log1p(pk.dist[,.SD,.SDcols=patterns('^dist')]), log1p(pk.dist[,.SD,.SDcols=patterns('^tot')]), type='b', col=rep(1:nrow(variables), each=length(land_grid_list)), lty=1, pch=seq(length(land_grid_list)),
#             ylab='peak intensity (E+I+C)', xlab='peak distance', main='Peak intensity vs. distance by timestep')

    inf.cells1 <- unique(wv.speed[is.na(edge),tot.inf.cells := length(tot.inf), by=.(v,l,r,time)][,.(v,l,r,time,tot.inf.cells)])
    inf.cells <- dcast(inf.cells1, time ~ v+l, value.var='tot.inf.cells', fun.aggregate=mean, fill=NA)
    matplot(log1p(inf.cells[,2:ncol(inf.cells)]), type='b', lty=1, col=rep(1:nrow(variables), each=length(land_grid_list)), pch=seq(length(land_grid_list)),
            ylab='log(1+tot. infected cells)', xlab='week', main='Infected cells over time')

#     inf.indiv1 <- unique(wv.speed[,tot.inf.indiv := sum(tot.inf), by=.(v,l,r,time)][,.(v,l,r,time,tot.inf.indiv)])
#     inf.indiv <- dcast(inf.indiv1, time ~ v+l, value.var='tot.inf.indiv', fun.aggregate=mean, fill=NA)
#     matplot(log1p(inf.indiv[,2:ncol(inf.indiv)]), type='b', lty=1, col=rep(1:nrow(variables), each=length(land_grid_list)), pch=seq(length(land_grid_list)),
#             ylab='log(1+tot. infected individuals E+I+C)', xlab='week', main='Infected individuals over time')

    par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), mar=c(0, 0, 0, 0), new=TRUE)
    plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
    legend("bottom", legend=c('vars:', 1:nrow(variables), 'land:', seq(length(land_grid_list))),
           lty=c(NA, rep(1,nrow(variables)), NA, rep(NA,length(land_grid_list))),
           col=c(NA, seq(nrow(variables)), NA, rep(1, length(land_grid_list))),
           pch=c(rep(NA, 2+nrow(variables)), seq(length(land_grid_list))),
           horiz=TRUE, bty='n', cex=1.3, lwd=2)

    dev.off()

    par(par.og)

#     ## Plot movement and density input parameters, with z response of spatial spread rate and maximum incidence
#     # get names of variables that are not connected to state (leaves state, density, B1, and any other user-defined variables)
#     wild.vars <- names(variables)[names(variables) %in% c('ss','B2','shape','rate','F1','F2_int','F2_B','F2i_int','F2i_B') == FALSE]
#     # grab unique combinations of the remaining variables
#     mvden.var <- unique(variables[,..wild.vars])
#
#     # connect to spatial spread rate
#     vl.avg.spd <- apply(spd.wave[,2:ncol(spd.wave)], 2, mean, na.rm=TRUE)
#     vl.avg.spd2 <- data.table(names(vl.avg.spd), vl.avg.spd)
#     vl.avg.spd2[,v := as.numeric(unlist(lapply(strsplit(V1, '_'), function(x) x[[1]])))]
#     vl.avg.spd2[,l := as.numeric(unlist(lapply(strsplit(V1, '_'), function(x) x[[2]])))]
#
#     # connect to maximum incidence
#     inf.indiv2 <- dcast(inf.indiv1, time ~ v+l, value.var='tot.inf.indiv', fun.aggregate=max, fill=NA)
#     vl.max.inc <- apply(inf.indiv2[,2:ncol(inf.indiv2)], 2, max, na.rm=TRUE)
#     vl.max.inc2 <- data.table(names(vl.max.inc), vl.max.inc)
#     vl.max.inc2[,v := as.numeric(unlist(lapply(strsplit(V1, '_'), function(x) x[[1]])))]
#     vl.max.inc2[,l := as.numeric(unlist(lapply(strsplit(V1, '_'), function(x) x[[2]])))]
#
#     vl.spd.inc <- vl.max.inc2[vl.avg.spd2, on=.NATURAL][,.(v, l, vl.max.inc, vl.avg.spd)]
#
#     var.spd.inc <- vl.spd.inc[variables[,v:=seq(nrow(variables))], on=.NATURAL]
#
#     png('./test_outputs/inc_spread.png', width=1200, height=1200)
#     par(mfrow=c(2,2))
#     plot(B1 ~ density, data=var.spd.inc[state=='FL'], col=rgb(vl.max.inc/max(vl.max.inc), 0, 0), pch=15, cex=10, main='FL max incidence')
#     plot(B1 ~ density, data=var.spd.inc[state=='SC'], col=rgb(vl.max.inc/max(vl.max.inc), 0, 0), pch=15, cex=10, main='SC max incidence')
#     plot(B1 ~ density, data=var.spd.inc[state=='FL'], col=rgb(vl.avg.spd/max(vl.avg.spd), 0, 0), pch=15, cex=10, main='FL avg spread rate')
#     plot(B1 ~ density, data=var.spd.inc[state=='SC'], col=rgb(vl.avg.spd/max(vl.avg.spd), 0, 0), pch=15, cex=10, main='SC avg spread rate')
#     dev.off()

    ## Proportion of simulations that ASF establishes

    ldsel.nnd <- readRDS('./Landscape_Setup/NND_Lands/4_Output/ldsel.rds')
    setDT(ldsel.nnd)
    setnames(ldsel.nnd, 'index', 'tileno')
    lands_names[,plandno := as.numeric(unlist(tstrsplit(unlist(tstrsplit(file, '_', keep=2)), fixed=TRUE, '.', keep=1)))]
    ldsel.nnd <- ldsel.nnd[lands_names, on=.(tileno=plandno)]

    # takes anything with (1) 2 or more cases at a given time, OR (2) at least 10 cases total
    ## should be establishment criteria from madison's paper
    tm.mat[,max.inf := max(I), by=.(var, land, rep)]
    tm.mat[,tot.exp := sum(E), by=.(var, land, rep)]
    tm.mat[,est := 0]
    tm.mat[max.inf > 1 & tot.exp >= 10, est := 1]
#     established <-  unique(tm.mat[,est.no := sum(est, na.rm=TRUE), by=.(var, land)][,est.pct := est.no/.N, by=.(var, land)][,.(var, land, est.no, est.pct)])
    established <- unique(unique(tm.mat[,.(var, land, rep, est)])[,est := sum(est)/.N, by=.(var, land)][,-'rep'])
    established <- established[variables[,.(contact, variant, density, id)], on=.(var=id)]
    est.lands <- established[ldsel.nnd, on=.(land=land), nomatch=NULL]


    browser()
#     established[, cat.name := paste(land, state, variant, density, sep='_')]
#     est.wide.pct <- dcast(established, var ~ land, value.var='est.pct')
#     est.wide.no <- dcast(established, var ~ land, value.var='est.no')
    png('./test_outputs/PCt_estab.png', width=1000, height=800)
    barplot(est ~ as.factor(land) + as.factor(var), data=established, beside=TRUE, legend=TRUE, names=established[,paste(contact, variant, density, sep='_')], ylim=c(0,1), main='% Establishment by land', xlab="Movement_Strain_Density", ylab='% established')
    dev.off()
    ## will want to have this as heatmap for landscape attributes for each of the 8 variable combinations

    unq.parms <- unique(est.lands[,.(var, contact, variant, density)])

    png('./test_outputs/land_parcombos_estab.png', width=1000, height=800)
    par(mfrow=c(length(unique(unq.parms[,density])), nrow(unq.parms)/length(unique(unq.parms[,density]))))
    lapply(seq(nrow(unq.parms)), function(y){
        sub.parms <- unq.parms[y,]
        plot(nnd_med ~ disp, data=est.lands[var == sub.parms[,var] & contact == sub.parms[,contact] & variant == sub.parms[,variant] & density == sub.parms[,density]],
             pch=15, col=rgb(est,0,0), cex=4, main=paste(unlist(sub.parms), collapse='_'))
    })
    dev.off()

    ## max observed incidence for var/land combos
    max.incidence <- tm.mat.edge[, max(incidenceEI, na.rm=TRUE), by=.(var, land)][unq.parms, on='var'][ldsel.nnd, on='land']
    setnames(max.incidence, 'V1','max.inc')
    max.incidence[, max.inc.norm := max.inc/max(max.inc)]
    png('./test_outputs/land_parcombos_maxincd.png', width=1000, height=800)
    par(mfrow=c(length(unique(unq.parms[,density])), nrow(unq.parms)/length(unique(unq.parms[,density]))))
    lapply(seq(nrow(unq.parms)), function(y){
        sub.parms <- unq.parms[y,]
        plot(nnd_med ~ disp, data=max.incidence[var == sub.parms[,var] & contact == sub.parms[,contact] & variant == sub.parms[,variant] & density == sub.parms[,density]],
             pch=15, col=rgb(max.inc.norm,0,0), cex=4, main=paste(unlist(sub.parms), collapse='_'))
    })
    dev.off()

    ## max.dist.wave
    wvdist.t60 <- unique(wv.speed[is.na(edge), maxtime := max(time), by=.(v,l)][time <=60,][time == 60 | time == maxtime, .(v,l,r,time, max)])[unq.parms, on=.(v=var)][ldsel.nnd, on=.(l=land)]
    wvdist.t60[,max.norm := max/max(max)]
    png('./test_outputs/land_parcombos_maxdist.png', width=1000, height=800)
    par(mfrow=c(length(unique(unq.parms[,density])), nrow(unq.parms)/length(unique(unq.parms[,density]))))
    lapply(seq(nrow(unq.parms)), function(y){
        sub.parms <- unq.parms[y,]
        plot(nnd_med ~ disp, data=wvdist.t60[v == sub.parms[,var] & contact == sub.parms[,contact] & variant == sub.parms[,variant] & density == sub.parms[,density]],
             pch=15, col=rgb(max.norm,0,0), cex=4, main=paste(unlist(sub.parms), collapse='_'))
    })
    dev.off()

    ## spd.wave
    wvspd.avgdistdiff <- unique(wv.speed[is.na(edge), mean(avg.dist.diff), by=.(v,l)])[unq.parms, on=.(v=var)][ldsel.nnd, on=.(l=land)]
    setnames(wvspd.avgdistdiff, 'V1','avg.spd')
    wvspd.avgdistdiff[,avg.spd.norm := avg.spd/max(avg.spd)]
    png('./test_outputs/land_parcombos_wvspd.png', width=1000, height=800)
    par(mfrow=c(length(unique(unq.parms[,density])), nrow(unq.parms)/length(unique(unq.parms[,density]))))
    lapply(seq(nrow(unq.parms)), function(y){
        sub.parms <- unq.parms[y,]
        plot(nnd_med ~ disp, data=wvspd.avgdistdiff[v == sub.parms[,var] & contact == sub.parms[,contact] & variant == sub.parms[,variant] & density == sub.parms[,density]],
             pch=15, col=rgb(avg.spd.norm,0,0), cex=4, main=paste(unlist(sub.parms), collapse='_'))
    })
    dev.off()

    ## inf.cells
    max.inf.cells <- wv.speed[is.na(edge), max(tot.inf.cells), by=.(v,l)][unq.parms, on=.(v=var)][ldsel.nnd, on=.(l=land)]
    setnames(max.inf.cells, 'V1', 'max.inf.cells')
    max.inf.cells[,max.inf.cells.norm := max.inf.cells/max(max.inf.cells)]
    png('./test_outputs/land_parcombos_infcells.png', width=1000, height=800)
    par(mfrow=c(length(unique(unq.parms[,density])), nrow(unq.parms)/length(unique(unq.parms[,density]))))
    lapply(seq(nrow(unq.parms)), function(y){
        sub.parms <- unq.parms[y,]
        plot(nnd_med ~ disp, data=max.inf.cells[v == sub.parms[,var] & contact == sub.parms[,contact] & variant == sub.parms[,variant] & density == sub.parms[,density]],
             pch=15, col=rgb(max.inf.cells.norm,0,0), cex=4, main=paste(unlist(sub.parms), collapse='_'))
    })
    dev.off()

    ## landscape attributes
    # land_grid_list
    # 1: number of cells
    # 2: table:
    #     1: cell id
    #     2: dim edge low
    #     3: dim edge high
    #     4: dim edge low
    #     5: dim edge high
    #     6: center point x?
    #     7: center point y?
    #     8: cell preference value
#     lapply(seq(length(land_grid_list)), function(x){
#         in.rast <- rast(land_grid_list[[x]][[3]], type='xyz')
#         corr <- autocor(in.rast, global=TRUE)
#         return(corr)
#     })
}
