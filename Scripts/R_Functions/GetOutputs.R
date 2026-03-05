
GetOutputs <- function(pop, centroids, BB, Incidence, Tculled, ICtrue, out, detectday, Ct, out.opts, input.opts){
  
  #List of outputs created here:	
	#Tinc #sum of all exposures over simulation 
	#sum(Tculled)
	#idT #last day there is an infectious individual
	#Mspread #max spread of infection
	#IConDD #number of I, C, and E on detection day
	#ICatEnd #number of I,C,E on last day 
	#TincFromDD #sum of all exposures starting day after detection day
	#TincToDD #sum of all exposures up until detection day
	#DET #total number of detections
	#allzone # all cells in monitoring zone

#to add:
  #toggle option: locs, dataframe of x/y locs for each sounder at each timestep, with SEIRCZ status
  #toggle option: cell numbers of all in zone at each timestep 
  #toggle options: POSlive, POSdead, POSlive_locs, POSdead_locs, but as dataframes
  #toggle options: num new infections at each step
  #toggle options: R0 at each step
  
    ### Summarize main outputs ###

    #Tinc, this just sums all of the exposures
    Tinc <- sum(Incidence)

    #sum all culled
    sumTculled <- sum(Tculled)

    if(any(ICtrue!=0)){
        #Find last day there was an infectious individual
        idT <- which(ICtrue!=0)[length(which(ICtrue!=0))]
        #IConDD #number of I, C, and E on detection day
        IConDD <- ICtrue[detectday]
        #ICatEnd #number of I,C,E on last day
        ICatEnd <- ICtrue[idT]
    } else {
        idT <- 1
        IConDD <- 0
        ICatEnd <- 0
    }
    #Find max spread of infection
    Mspread <- max(out[,2])

    #TincToDD #sum of all exposures up until detection day
    TincToDD <- sum(Incidence[idT:detectday])

    #TincFromDD #sum of all exposures starting day after detection day
    TincFromDD <- sum(Incidence[detectday:idT])

    #DET #total number of detections
    DET <- NA

    ### Compile main outputs ###

    #send all main outputs to list
    list.all <- list("Tinc" = Tinc,
                "sumTculled" = sumTculled,
                "BB" = BB,
                "Mspread" = Mspread,
                "IConDD" = IConDD,
                "ICatEnd" = ICatEnd,
                "TincToDD" = TincToDD,
                "TincFromDD" = TincFromDD,
                "DET" = DET,
                "Ct" = Ct,
                "pop" = pop)

    ### Summarize optional outputs ###

    #loc.list: is pop[,c(3,8:13)] for each timestep in list of length thyme
    #cell num, S, E, I, R, C, Z
    #convert loc.list into dataframe

    if("sounderlocs" %in% out.opts){
        loc.list <- input.opts$loc.list
        input.rng <- 1:length(loc.list)
        input.rng <- input.rng[lapply(loc.list, length) != 0]

        locs.df <- rbindlist(lapply(input.rng, function(i){
            # get locations
            locs.i <- as.data.table(centroids[matrix(loc.list[[i]], ncol=7)[,1], , drop=FALSE])
            # add timestep and state variable values to locs.i
            locs.i <- cbind(i, locs.i, loc.list[[i]][, 2:7])
            # set names for output
            colnames(locs.i) <- c("time", "x", "y", "unknown", 'S', 'E', 'I', 'R', 'C', 'Z')
            return(locs.i)
        }))
    }

    #Get POSlive, dead, and all locs in datafame format from list
    if("alldetections" %in% out.opts){
        #all should be same length, thyme
        POSlive <- input.opts$POSlive
        POSlive_locs <- input.opts$POSlive_locs
        POSdead <- input.opts$POSdead
        POSdead_locs <- input.opts$POSdead_locs
        allzone <- input.opts$allzonecells

        detections <- rbindlist(lapply(1:length(POSlive), function(i){
            live.detections.i <- data.table(i, 1, POSlive[[i]], POSlive_locs[[i]])
            # Create dead detections after live detections are populated
            dead.detections.i <- data.table(i, 0, POSdead[[i]], POSdead_locs[[i]])
            # Combine live and dead detections into detections.i
            detections.i <- rbind(live.detections.i, dead.detections.i)
            return(detections.i)
        }))

        # Now, outside of the loop, check if sample == 1 to combine live detections with sampled pigs
        if (sample == 1) {
            # Get the sampled pigs for the timestep
            pigs_sampled_timestep <- input.opts$pigs_sampled_timestep
            ## check these hard-coded values
            sampled_pigs_column <- matrix(unlist(pigs_sampled_timestep), nrow = 52, ncol = 1)  # Convert to a matrix (1 row, 52 columns)
            sampled_pigs_column_dup <- matrix(rep(sampled_pigs_column, each = 2), ncol = 1)
            # Combine the live detections with the sampled pigs as a new column
            detections <- cbind(detections, sampled_pigs_column_dup)
        }
    }
    #Get Incidence and R0 vals summarized in data frame
    if("incidence" %in% out.opts){
        #Timestep
        #Num new exposures from infected (incidence)
        #Num infected (prev timestep)
        #Rt (total, homog) at each time step
        inc.mat <- input.opts$incidence
        inc.df <- as.data.frame(inc.mat)
        colnames(inc.df) <- c("timestep", "state", "loc")
        inc.df[inc.df$state == 10, 2] <- "infected" #any infected that could have exposed susceptibles
        inc.df[inc.df$state == 9, 2] <- "exposed" #exposed this timestep
        inc.df[inc.df$state == 12, 2] <- "carcass" #exposed this timestep
    }

    ### Append optional outputs ###

    if("sounderlocs" %in% out.opts){
        templist <- vector(mode="list", length=1)
        templist[[1]] <- locs.df
        list.all <- append(list.all, templist)
        names(list.all)[length(list.all)] <- "sounderlocs"
    }

    if("idzone" %in% out.opts){
        templist <- vector(mode="list", length=1)
        templist[[1]] <- input.opts$idzone.mat
        list.all <- append(list.all, templist)
        names(list.all)[length(list.all)] <- "idzone"
    }

    if("alldetections" %in% out.opts){
        templist <- vector(mode="list", length=1)
        templist[[1]] <- detections
        list.all <- append(list.all, templist)
        names(list.all)[length(list.all)] <- "alldetections"

        templist <- vector(mode="list", length=1)
        templist[[1]] <- input.opts$allzonecells
        list.all <- append(list.all, templist)
        names(list.all)[length(list.all)] <- "allzonecells"
    }

    if("incidence" %in% out.opts){
        templist <- vector(mode="list", length=1)
        templist[[1]] <- inc.df
        list.all <- append(list.all, templist)
        names(list.all)[length(list.all)] <- "incidence"
    }

    return(list.all)
}
