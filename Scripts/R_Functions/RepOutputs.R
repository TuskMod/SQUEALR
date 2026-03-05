## repackage outputs after simulation/burn-in run

rep_outputs <- function(out.list, v, l, r, parameters, out.opts, prevrep.in = as.list(rep(NA, 5))){

    list2env(parameters, .GlobalEnv)

    # if it is the burnin, the r=0 below will overwrite NA defaults in this irrelevant
    # if it is after the burn-in, these will have existing values that new values will be added to
    tm.mat <- prevrep.in[[1]]
    summ.vals <- prevrep.in[[2]]
    incidence <- prevrep.in[[3]]
    detections <- prevrep.in[[4]]
    allzone <- prevrep.in[[5]]

    ## test these outputs if out.opts doesn't include them
    ## id by variable combination, landscape, rep, and timestep for one row per timestep data
    end.tm <- out.list$endtime
    id.r <- c(v,l,r)
    tm.mat.r <- cbind(matrix(id.r, nrow=end.tm, ncol=3, byrow=TRUE), seq(end.tm))
    colnames(tm.mat.r) <- c("var","land","rep","timestep")

    #Handle effective removal rate (timestep output)
    tm.mat.r <- cbind(tm.mat.r, out.list$Ct[seq(end.tm)])
    colnames(tm.mat.r)[ncol(tm.mat.r)] <- "Ct"

    # Births
    tm.mat.r <- cbind(tm.mat.r, out.list$BB[seq(end.tm)])
    colnames(tm.mat.r)[ncol(tm.mat.r)] <- 'BB'

    #Handle sounderlocs (optional...)
    ## seems like other things were supposed to happen in sounderlocsSummarize, if we want to use those this will have to change
    if ("sounderlocs" %in% out.opts){
        solocs.r <- sounderlocsSummarize(out.list$sounderlocs, r)[[1]]
        solocs.all <- cbind(v, l, r, out.list$sounderlocs)
        tm.mat.r <- cbind(tm.mat.r, solocs.r[3:8])
    }

    # detections (optional...)
    ## has a row for each timestep AND detection type, with timestep, code (1=live,0=dead), number of individuals detected, and position
    ## still need to test sample = 1
    if (end.tm > detectday & 'alldetections' %in% out.opts){
        detections.r <- out.list$alldetections
        detections.r <- detections.r[detections.r[,1] <= end.tm & detections.r[,1] >= detectday,]
        n.det <- nrow(detections.r)
        detections.r <- suppressWarnings(cbind(matrix(id.r, ncol=3, nrow=n.det, byrow=TRUE), detections.r)) # gave a warning if there were no detections; very annoying
        colnames(detections.r) <- c('var', 'land', 'rep', 'timestep', 'code', 'detected', 'loc')
        allzone.r <- out.list$allzonecells
        colnames(allzone.r) <- c('var', 'land', 'rep', 'timestep', 'loc')
    }

    # incidence -- more rows than timesteps, separate output (optional...)
    if ('incidence' %in% out.opts){ # don't want this if it's a burn-in output
        incidence.r <- out.list$incidence
        n.inc = nrow(incidence.r)
        incidence.r <- suppressWarnings(cbind(matrix(id.r, ncol=3, nrow=n.inc, byrow=TRUE), incidence.r))
        colnames(incidence.r) <- c('var', 'land', 'rep', 'timestep', 'state', 'loc')
    }

    # single value per vlr combination outputs
    # Tinc # sumTculled # Mspread # IConDD # ICatDD # TincToDD # TincFromDD # DET (total detections) #
    summ.vals.r <- matrix(c(id.r, unlist(out.list[c(1, 2, 4:9, length(out.list))])), nrow=1)

    colnames(summ.vals.r) <- c('var','land','rep',names(out.list[c(1, 2, 4:9, length(out.list))]))

    # Put new results in output objects OR add new results to existing output objects
    if(r==1 & l==1 & v==1){
        tm.mat <- tm.mat.r
        summ.vals <- summ.vals.r
        if ("incidence" %in% out.opts) incidence <- incidence.r
        if ("alldetections" %in% out.opts) {
            detections <- matrix(NA, ncol=7, nrow=0)
            colnames(detections) <- c('var', 'land', 'rep', 'timestep', 'code', 'detected', 'loc')
            allzone <- matrix(NA, ncol=5, nrow=0)
            colnames(allzone) <- c('var', 'land', 'rep', 'timestep', 'loc')
        }

    } else {
        tm.mat <- rbind(tm.mat, tm.mat.r)
        tm.mat <- tm.mat[!is.na(tm.mat[,1]),, drop=FALSE]
        summ.vals <- rbind(summ.vals, summ.vals.r)
        summ.vals <- summ.vals[!is.na(summ.vals[,1]),, drop=FALSE]
        if ("incidence" %in% out.opts) {
            incidence <- rbind(incidence, incidence.r)
            incidence <- incidence[!is.na(incidence[,1]),, drop=FALSE]
        }
        if (end.tm > detectday & "alldetections" %in% out.opts) {
            detections <- rbind(detections, detections.r)
            detections <- detections[!is.na(detections[,1]),, drop=FALSE]
            if (is.null(nrow(allzone))) {
                # if it hasn't been established yet
                allzone <- allzone.r
            } else {
                allzone <- rbind(allzone, allzone.r)
                allzone <- allzone[!is.na(allzone[,1]),, drop=FALSE]
            }
            colnames(allzone) <- c('var', 'land', 'rep', 'timestep', 'loc')
        }
    }

    return(list(tm.mat, summ.vals, incidence, detections, allzone, solocs.all))
}
