
CullingOneRun <- function(pop, idNEW, idZONE, Intensity, alphaC, centroids, Rad, inc, i, POSlive, POSdead, POSlive_locs, POSdead_locs, NEGlive, NEGdead, DetP, cullstyle){
    cells <- nrow(centroids)

    #if there were infected pigs detected in the last time step
    if(length(idNEW) > 0){
        #initiate vector to store pairs of each infected grid cell ID with every other cell ID
        #nrow for each combo, 3 cols
        #col 1-each grid cell ID with detected infection
        #col 2-each paired grid cell
        #col 3-distance between infected grid cell ID with detection and each grid cell
        ## i *think* this is causing a problem with integer overflow when things get too large (nrow is doing that somewhere in this file)
        pairedIDs <- matrix(nrow = cells * length(idNEW), ncol=3)
        #get matrix of each infected grid cell ID paired with every other cell ID
        if(length(idNEW) == 1){
            pairedIDs[seq(cells), 1] <- idNEW
            pairedIDs[seq(cells), 2] <- seq(cells)
        } else if(length(idNEW) > 1){
            for(j in 1:length(idNEW)){
                if(j == 1){
                    pairedIDs[seq(cells), 1] <- idNEW[j]
                    pairedIDs[seq(cells), 2] <- seq(cells)
                }
                if(j > 1){
                    cells=nrow(centroids)
                    pairedIDs[(1 + (j - 1) * cells):(j * cells), 1] <- idNEW[j]
                    pairedIDs[(1 + (j - 1) * cells):(j * cells), 2] <- seq(cells)
                }
            }

        }
        #get distance between each infected grid cell ID paired with every other cell ID
        #below results in matrix with one column of each grid cell id and the second column the distance to each infected
        ## I think this could be handled with a cached matrix of distances between cells that is calculated once
        pairedIDs[,3] <- sqrt((centroids[pairedIDs[, 1], 1] - centroids[pairedIDs[, 2], 1])^2 + (centroids[pairedIDs[, 1], 2] - centroids[pairedIDs[, 2], 2])^2)
        #get all grid cells where the distance between that grid cell and an infected detected pig is less than the pre-determined radius
        idout <- pairedIDs[pairedIDs[,3] <= Rad,]

    } else { #else, if no infections detected in last time step...no new grid cells added
        idout <- matrix(nrow=0, ncol=3)
    }

    #fullZONE contains all grid cells with detected infections, from last time step and all before
    # first column is where it is detected, 2nd column is distance from detected point by each other cell
    fullZONE <- rbind(idZONE, idout)

    #get all unique grid cells in the zone
    allINzone <- unique(fullZONE[, 2])
    allINzone <- allINzone[!is.na(allINzone)] ## clean out the rows of NA values, ain't nobody got time for that

    #get total number of grid cells in the zone
    Uall <- length(allINzone)

    #get total area of the zone
    ZONEkm2 <- Uall * inc^2

    #get which rows in the zone
    #total number of sounders=soundINzone
    soundINzone <- which(pop[,3] %in% allINzone)

    #get total number of pigs (live and dead) in zone
    pigsinzone <- sum(pop[soundINzone, 1],pop[soundINzone, 12],pop[soundINzone ,13])

    #total number of exposed/infected/infectious carcass pigs in zone
    EICinzone <- sum(pop[soundINzone, 9], pop[soundINzone, 10], pop[soundINzone, 12])

    #get total number of pigs outside the zone
    pigsoutzone <- sum(pop[-soundINzone, 1], pop[-soundINzone, 12], pop[-soundINzone, 13])

    #get total number of infected pigs outside the zone
    EICoutzone <- sum(pop[-soundINzone, 9], pop[-soundINzone, 10], pop[-soundINzone, 12])

    #get total number of  individuals (inside and outside zone)
    totalpigs <- sum(pop[, 1], pop[, 12], pop[, 13])

    #get total number of infected individuals (inside and outside zone)
    totalEIC <- EICinzone + EICoutzone

    #For output to get landscape-level effective removal rate
    Ct <- pigsinzone / totalpigs

    #####################################
    ###### Begin Culling Algorithm ######
    #####################################
    #if there are pigs to cull...
    # if(is.na(pigsinzone)){pigsinzone=0} ## if there are na's here, it means there are issues further up that probably shouldn't be glossed over...
    if(pigsinzone > 0){
        #get number of pigs for each grid cell in zone
        #and get their status, SEIRCZ
        #initiate empty matrix, nrow for each grid cell, 7 for each of SEIRCZ
        popINzone <- pop[soundINzone,, drop=FALSE]
        # make usable
        fullZONE.dt <- as.data.table(fullZONE)
        # clean out in-zone points that have multiple detection cells; keep the id and distance of the closest detection
        fullZONE.dt <- unique(fullZONE.dt[fullZONE.dt[, min(V3), by=V2], on=.(V2, V3=V1)])
        setnames(fullZONE.dt, c('V1','V2','V3'), c('det.cell','cell','dist'))
        # convert popINzone to a data table and aggregate state variables by summation
        popINzone.agg <- as.data.table(popINzone[,c(3, 8:13) ,drop=FALSE])[, lapply(.SD, sum), by=cell, .SDcols=2:7]
        # connect popINzone table with cell distances; clean out rows by cell id which have multiple detection points at the same distance (just keep one, doesn't matter which for this)
        fullZONEpigs <- fullZONE.dt[popINzone.agg, on=.(cell = cell)][, .SD[1], by=cell]

        fullZONEpigs[, cell.total := S+E+I+R+C+Z]
        fullZONEpigs[, cprob := 1 - 1/((1 + ..alphaC)^(cell.total/(..inc^2)))]
        norm.val <- TwoDt(0, 3, Rad/2)
        fullZONEpigs[,dist.weight := TwoDt(dist, 3, ..Rad/2)/..norm.val]

        fullZONEpigs[,cull.cell := rbinom(nrow(fullZONEpigs), cell.total, cprob * dist.weight)]
        setorder(fullZONEpigs, dist, -cell.total)
        fullZONEpigs[cull.cell > 0, cum.cull := round(cumsum(cell.total)/fullZONEpigs[, sum(cell.total)], 2)]
        # if the cumulative proportion of culling starting from the center hits the target 0.05 exactly, stick with that
        if (any(fullZONEpigs[, cum.cull] == 0.05, na.rm=TRUE)){
            target.value <- 0.05
        } else if (max(fullZONEpigs[, cum.cull], na.rm=TRUE) < 0.05){
            # if the cumulative proportion of culling doesn't reach 0.05, just take everything marked for culling by cull.cell
            target.value <- max(fullZONEpigs[,cum.cull], na.rm=TRUE)
        } else if (min(fullZONEpigs[,cum.cull], na.rm=TRUE) > 0.05 & sum(fullZONEpigs[,cell.total]) < 20){
            # if the lowest cumulative proportion of cullable cells is too big AND it is because the total pigs in the zone is less than what would be necessary to make ANY cell contain  less than 5% of the total pigs, just take the cell closest to the source as the culled locations
            # keeps the culling going when it starts to become the victim of its own success, since with an increasing area and decreasing density you will eventually have a single pig making up more than 5% of the total in the cell
            target.value <- min(fullZONEpigs[,cum.cull], na.rm=TRUE)
        } else {
            # if the cumulative culling proprotion goes above 0.05 and no cell lands it exactly at 0.05, take up to the first cell that pushes
            # the cumulative culling proportion above 0.05
            target.value <- fullZONEpigs[cum.cull > 0.05,][1,cum.cull]
        }
        # mark the selected cells for culling
        fullZONEpigs[cull.cell > 0 & cum.cull <= target.value, actual.cull := 1]

        if (sum(fullZONEpigs[,actual.cull], na.rm=TRUE) > 0 & nrow(fullZONEpigs) > 0){
            culled <- sum(fullZONEpigs[actual.cull == 1, cell.total])
        #get which cells culled
            removalcells <- as.matrix(fullZONEpigs[actual.cull == 1, cell])
        #get which pigs culled
            removalpigs <- as.matrix(fullZONEpigs[actual.cull == 1, ])

        ###### Update surveillance data ######

        #POSlive_i is a matrix with a row for each timestep
        #column one of poslive is the number of exposed/infected pigs detected at that timestep
        #sum removalpigs column 6,7
            POSlive_i <- sum(removalpigs[,5], removalpigs[,6])
            if(!is.null(DetP) & DetP != 1){
                POSlive_i_sel <- rbinom(nrow(removalpigs), rowSums(removalpigs[, 5:6]), DetP)
                POSlive_i <- sum(removalpigs[POSlive_i_sel, 5:6])
            } else if(DetP == 1){
                POSlive_i_sel <- as.numeric(which(rowSums(removalpigs[, 5:6, drop=FALSE]) > 0))
            }
        #POSdead is a matrix with a row for each timestep
        #column one of poslive is the number of infected carcasses detected at that timestep
        #sum removalpigs column 9
            POSdead_i <- sum(removalpigs[, 8])
            if(!is.null(DetP) & DetP != 1){
                POSdead_i_sel <- rbinom(nrow(removalpigs), removalpigs[, 8], DetP)
                POSdead_i <- sum(removalpigs[POSdead_i_sel, 8])
            } else if (DetP==1){
                POSdead_i_sel  <-  as.numeric(which(removalpigs[, 8] > 0))
            }
        #list of length thyme, each timestep is vector of grid cell locations where live infected pigs detected at that ts
            if(POSlive_i > 0){
                POSlive_locs_i <- removalpigs[POSlive_i_sel, 1]

            } else { POSlive_locs_i <- 0 }

            if(POSdead_i > 0){
                POSdead_locs_i <- removalpigs[POSdead_i_sel, 1]
            } else { POSdead_locs_i <- 0 }

            #vector of nrow timestop, count of total SR removed
            NEGlive_i <- sum(removalpigs[,4], removalpigs[,7])
            if(!is.null(DetP) & DetP != 1){
                NEGlive_i_missed <- length(POSlive_i_sel[POSlive_i_sel==0])
                NEGlive_i <- NEGlive_i + NEGlive_i_missed
            }
            #vector of nrow timestep, count of total Z removed
            NEGdead_i <- sum(removalpigs[,9])
            if(!is.null(DetP) & DetP != 1){
                NEGdead_i_missed <- length(POSdead_i_sel[POSdead_i_sel == 0])
                NEGdead_i <- NEGdead_i + NEGdead_i_missed
            }

            #remove removed sounders from pop
            removalrows <- which(pop[,3] %in% removalcells)
            removedpop <- pop[-removalrows,, drop=FALSE]
        } else{
            POSlive_i <- 0
            POSdead_i <- 0
            POSlive_locs_i <- NA
            POSdead_locs_i <- NA
            NEGlive_i <- 0
            NEGdead_i <- 0
            culled <- 0
            Ct <- 0
            removedpop <- pop
        }
        #idZONE:
        #grid cell ids that had a positive detection, grid cell ids that are within the zone, distance
        idZONE <- fullZONE

        #send updated objects to output list
        output.list <- vector(mode="list", length=13)
        output.list[[1]] <- POSlive_i
        output.list[[2]] <- POSdead_i
        output.list[[3]] <- POSlive_locs_i
        output.list[[4]] <- POSdead_locs_i
        output.list[[5]] <- NEGlive_i
        output.list[[6]] <- NEGdead_i
        output.list[[7]] <- idZONE
        output.list[[8]] <- removalcells
        output.list[[9]] <- culled
        output.list[[10]] <- ZONEkm2
        output.list[[11]] <- removedpop
        output.list[[12]] <- Ct
        output.list[[13]] <- allINzone
    } else {
        output.list <- vector(mode="list", length=13)
        output.list[[1]] <- 0
        output.list[[2]] <- 0
        output.list[[3]] <- 0
        output.list[[4]] <- 0
        output.list[[5]] <- 0
        output.list[[6]] <- 0
        output.list[[7]] <- idZONE
        output.list[[8]] <- 0
        output.list[[9]] <- 0
        output.list[[10]] <- 0
        output.list[[11]] <- pop
        output.list[[12]] <- 0
        output.list[[13]] <- allINzone
    }

    return(output.list)

} #function closing bracket

