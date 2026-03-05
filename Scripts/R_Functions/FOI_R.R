FOI_R <- function(pop, centroids, cells, B1, B2, F1, F2_int, F2_B, F2i_int, F2i_B){
    #Create I/C matrices. Just want vector with locations of I/C, with nums of I/C in right place
    I <- matrix(0, nrow=cells, ncol=1)
    C <- matrix(0, nrow=cells, ncol=1)
#     I[pop[which(pop[, 10] != 0), 3], 1] <- pop[which(pop[, 10] != 0), 10]
#     C[pop[which(pop[, 12] != 0), 3], 1] <- pop[which(pop[, 12] != 0), 12]
    ## The way this is structured, if there are more than one sounder in a cell, it will take the last nonzero value of the state variable from pop and use that as the only value for that point in I or C vectors.
    ## E.g.
#     > pop[pop[,3] == 17701,]
#          Nlive pref  cell     dist  ctrx  ctry pcell S E I R C Z
#     [1,]     4 0.51 17701 1.807427 49.25 43.25 18905 0 0 4 0 0 0
#     [2,]     2 0.51 17701 1.309366 51.25 43.25 17104 1 0 1 0 0 0
#     [3,]     2 0.51 17701 0.000000 47.25 45.75 17701 2 0 0 0 0 0
#
#     > I[pop[pop[,3] == 17701,3]]
#     [1] 1 1 1

#     > I[17700:17703]
#     [1] 0 1 0 0
    ## The I = 4 row is always erased because it is above the row with I=2
    ## Instead, we can add up the state variables for each cell to get FOI
    pop.temp <- as.data.table(pop)
    # sum up I and c by cell
    pop.temp <- pop.temp[,lapply(.SD, sum), .SDcols=c(10,12), by=.(cell)]
    # add to I and C matrices to plug into existing functions
    I[pop.temp[,cell], 1] <- pop.temp[,I]
    C[pop.temp[,cell], 1] <- pop.temp[,C]

    Pse <- matrix(0, nrow=cells, ncol=1)
    temp <- I + C
    id <- which(temp > 0)
    W <- matrix(0, nrow=cells, ncol=2)

    W[, 1] <- F1*I
    W[, 2] <- F1*C

    if(length(id) == 1){
        X1 <- matrix(nrow=length(id), ncol=1)
        Y1 <- matrix(nrow=length(id), ncol=1)

        X1 <- centroids[id, 1]
        Y1 <- centroids[id, 2]

        dist <- sqrt((centroids[, 1] - X1)^2 + (centroids[, 2] - Y1)^2)
        dist <- as.data.frame(dist)
        colnames(dist) <- "X"
        # direct (live-live) contact with infectious pigs in one cell and susceptible pigs in another
        prob <- F2_int + F2_B * dist
        # indirect (dead-live) contact with infectious dead pigs in one cell and susceptible pigs in another
        probi <- F2i_int + F2i_B * dist
        ## not sure about this below... doesn't this mean sounders in the same cell cannot pass to each other??
        prob[dist == 0] <- 0
        probi[dist == 0] <- 0

        B <- B1 * I[id[1]] * prob + B2 * C[id[1]] * probi
        Pse <-  1 - exp( - W[, 1] - W[, 2] - t(B))
    }

    if(length(id)>1){
        B <- Fast_FOI_function((id - 1), centroids[, 1, drop=FALSE], centroids[, 2, drop=FALSE],
                            cells, F2_int, F2_B, F2i_int, F2i_B, I[id, , drop=FALSE], C[id, , drop=FALSE], B1, B2)
        Pse <-  1 - exp( - W[, 1] - W[, 2] - t(colSums(B)))
    }

    return(Pse[pop[,3]])
  
}
