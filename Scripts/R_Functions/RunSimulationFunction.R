
### Initialize grid(s): ---------------
    #Method for class 'list'
        #Initialize_Grids(object, parameters=parameters, movement=c(parameters$shape,parameters$rate))
    #Arguments
                #movement
                    #default is vector with gamma distribution shape and rate fed from parameters file

      #grid.opts- 
        #"homogenous" or "heterogeneous", default for class numeric is "homogenous" and default for type SpatRaster and SpatRasterCollection is "ras"
          #ras- 
            #use input raster to create grid
          #homogeneous-
            #creates grid with 7 columns, no land class variables
          #heterogeneous-
            #creates a neutral random landscape model with X lc variables
    #Value
      #a nested list of grid parameters
RunSimulationReplicates <- function(land_grid_list, parameters, variables, cpp_functions, reps, mv.parms){

    names(variables)[names(variables) == "density"] <- "dens"

    # Filter variables (parameters with >1 value) out of parameters
    ## selecting variables is done in SetVarParms.R
    parameters <- parameters[names(parameters) %in% names(variables) == FALSE]
    parameters <- parameters[-grep('.\\_\\_.', names(parameters))]

    #Pull needed parms from parameters for all reps
    list2env(parameters, .GlobalEnv)

    # looping table for mapply
    lvtable <- expand.grid(vars = seq(nrow(variables)), land = 1:length(land_grid_list), rep = seq(reps))
#     lvtable <- expand.grid(vars = seq(nrow(variables)), land = seq(length(land_grid_list)), rep = seq(reps))

    # movement parameters from NND landscape selection
    setDT(mv.parms)


    # loops over combinations of variables, lands, and reps
    rep.list <- mapply(function(v.val, l.val, r.val){

        # add in vars
        vars <- variables[v.val,]
        vars <- as.list(vars)
        list2env(vars, .GlobalEnv)

        #calc vals based on variables
        N0 <- dens*area
        K <- N0*1.5
        parameters <- c(parameters, vars)
        parameters$K <- K

        # movement parameters
        parameters$shape <- as.numeric(mv.parms[l.val, gamma.shape])
        parameters$rate <- as.numeric(mv.parms[l.val, gamma.rate])

        #loop through landscapes
        centroids <- land_grid_list[[l.val]]$centroids
        grid <- land_grid_list[[l.val]]$grid

        # create sounders in starting locations according to N0 and ss parameters
        pop <- InitializeSounders(centroids, grid, c(N0, ss), pop_init_grid_opts)
        # add an infected individual near the center of the simulation space
        pop <- InitializeInfection(pop, centroids, grid, parameters)
        # pre-create outputs to catch output data
        outputs <- Initialize_Outputs(parameters)
        # Run simulation
        out.list <- SimulateOneRun(outputs, pop, centroids, grid, parameters, cpp_functions, K, v.val, l.val, r.val)
        # Handle outputs
        rep.out <- rep_outputs(out.list, v.val, l.val, r.val, parameters, out.opts)
        return(rep.out)
    },
    v.val=lvtable[,1], l.val=lvtable[,2], r.val = lvtable[,3])

    tm.mat <- rbindlist(rep.list[1,])
    summ.vals <- rbindlist(lapply(rep.list[2,], as.data.table))
    incidence <- rbindlist(rep.list[3,which(lapply(rep.list[3,], ncol) > 1)])
    detections <- rbindlist(lapply(rep.list[4,][!is.na(rep.list[4,])], as.data.table))
    allzone <- rbindlist(lapply(rep.list[5,][!is.na(rep.list[5,])], as.data.table))
    solocs.all <- rbindlist(rep.list[6,])

    return(list('tm.mat' = tm.mat, 'summ.vals' = summ.vals, 'incidence' = incidence, 'detections' = detections, 'allzone' = allzone, 'solocs.all' = solocs.all))#, 'wv.speed' = wv.speed))
}










