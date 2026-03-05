InitializeInfection <- function(pop, centroids, grid, parameters){
	
######################################
######## Initialize Infection ######## 
######################################
#num_inf_0=1 #how many pigs to infect starting off
	num_inf_0 <- parameters$num_inf_0

	#find the midpoint of the grid where the infected sounder will end up
	midpoint <- c(median(centroids[, 1]), median(centroids[, 2]))
	id <- which(centroids[, 1] >= midpoint[1] & centroids[, 2] >= midpoint[2])[1] #location on grid closest to midpoint

    # generate a sounder of one infected individual at num_inf_0 points (usually 1)
    infected <- InitializeSounders(centroids, grid, c(id, num_inf_0), pop_init_type="init_single", pop_init_grid_opts="homogeneous")
    # (manually change state values so S=0 and I = 1)
    infected[, 8] <- 0
    infected[, 10] <- 1
    ## could have the option to change initial infected sounder size
    #combine infected pig with pop matrix
    pop <- rbind(pop, infected)

	return(pop)
}
