#outputs: pop,Incidence,BB
StateChanges <- function(pop, centroids, cells, parameters, Incidence, BB, i){

    list2env(parameters, .GlobalEnv)

    # should have list2env for parameters/outputs here for consistency (instead of 1K+ inputs)
    # births
    Sdpb <- matrix(nrow=nrow(pop),ncol=1)
    Sdpb[,1] <- 0

    # natural deaths
    Sdpd <- matrix(nrow=nrow(pop), ncol=1)
    Edpd <- matrix(nrow=nrow(pop), ncol=1)
    Rdpd <- matrix(nrow=nrow(pop), ncol=1)
    Sdpd[,1] <- 0
    Edpd[,1] <- 0
    Rdpd[,1] <- 0

    # disease state change recording
    Eep <- matrix(nrow=nrow(pop), ncol=1)
    Eep[, 1] <- 0
    Iep <- matrix(nrow=nrow(pop), ncol=1)
    Iep[, 1] <- 0
    Rep <- matrix(nrow=nrow(pop), ncol=1)
    Rep[, 1] <- 0
    Cep <- matrix(nrow=nrow(pop), ncol=1)
    Cep[, 1] <- 0

    # Carcass decay recording
    Ccd <- matrix(nrow=nrow(pop), ncol=1)
    Ccd[, 1] <- 0
    Zcd <- matrix(nrow=nrow(pop), ncol=1)
    Zcd[, 1] <- 0


    ########### Determine Births ###########
    # subset sounder sets with live, uninfected individuals
    idN <- pop[pop[, 8, drop=FALSE] > 0 | pop[, 9, drop=FALSE] > 0 | pop[, 11, drop=FALSE] > 0,, drop=FALSE]

    # Number of live, uninfected individuals
    liveind <- sum(colSums(pop)[c(8, 9, 11)])

    # get row indices of live individuals
    liverows <- which(pop[, 8, drop=FALSE] > 0 | pop[, 9, drop=FALSE] > 0 | pop[, 11, drop=FALSE] > 0) #rownums with live indiv

    # density-dependent birth rate
    ## this will need adjustment for culling so unculled areas don't get a growth boost artificially
    Brate <- Pbd * rowSums(pop[, c(8, 9, 11)]) * (1 - sum(pop[, c(8, 9, 11)]) / K)
    Sdpb <- pmin(10 * rowSums(pop[, c(8, 9, 11)]), rnbinom(nrow(pop), mu = Brate, size = 0.01))

    # get total births, using Brate as mean in a poisson
    BB[i] <- sum(Sdpb)

    ######## Determine disease state change probabilities ########
    # Susceptible to Exposed
    Pse <- FOI_R(pop, centroids, cells, B1, B2, F1, F2_int, F2_B, F2i_int, F2i_B) #cpp parallel version, 22x faster than R version
    Pse[Pse < 0] <- 0

    # Exposed to infected -- based on the incubation period
    Pei <- 1 - exp(-1 / (rpois(nrow(pop), incub) / 7))

    # Infected to contageous corpse or recovered, depending on Pir (recovery rate)
    Pic <- 1 - exp(-1 / (rpois(nrow(pop), infpd) / 7))

    ######## Conduct the State Transitions ########
    # susceptible state changes
    Sdpd <- rbinom(nrow(pop), pop[, 8], death)
    Eep <- rbinom(nrow(pop), pop[, 8] - Sdpd, Pse)

    # exposed state changes
    Edpd <- rbinom(nrow(pop), pop[, 9], death)
    # exposed -> infected
    Iep <- rbinom(nrow(pop), pop[, 9] - Edpd, Pei)

    # infected state changes
    Rep <- rbinom(nrow(pop), pop[, 10], Pir * Pic)
    # infected -> contageous carcass
    Cep <- rbinom(nrow(pop), pop[, 10] - Rep, (1-Pir)*Pic)

    # recovered state changes
    Rdpd <- rbinom(nrow(pop), pop[, 11], death)

    # infected carcass state changes
    Ccd <- rbinom(nrow(pop), pop[, 12], Pcr)

    # uninfected carcass state changes
    Zcd <- rbinom(nrow(pop), pop[, 13], Pcr)

    # collect incidence info
    Incidence[i] <- Incidence[i] + sum(Eep)

    ######## Update pop matrix ########
    pop[,8] <- pop[,8] - Eep + Sdpb - Sdpd #S
    pop[,9] <- pop[,9] - Iep + Eep - Edpd #E
    pop[,10] <- pop[,10] - Rep - Cep + Iep #I
    pop[,11] <- pop[,11] + Rep - Rdpd #R
    pop[,12] <- pop[,12] + Cep - Ccd #C
    pop[,13] <- pop[,13] + Sdpd + Rdpd + Edpd - Zcd #Z

    #move dead individuals (C or Z) into their own rows
    deadguys <- pop[pop[, 12] > 0 | pop[, 13] > 0,, drop=FALSE]

    # if there are deadguys....
    if(nrow(deadguys) != 0){
        # remove abundance and all live guy counts from deadguy set
        deadguys[, 1] <- 0
        deadguys[, 8] <- 0
        deadguys[, 9] <- 0
        deadguys[, 10] <- 0
        deadguys[, 11] <- 0

        # set all deadguys in pop rows to zero
        pop[, 12][pop[, 12] > 0] <- 0
        pop[, 13][pop[, 13] > 0] <- 0

        # add deadguys to pop matrix
        pop <- rbind(pop, deadguys)
    }

    # Update abundance numbers (live individuals only count in abundance)
    pop[,1] <- rowSums(pop[, 8:11, drop=FALSE])

    # eliminate anything where there are no pigs or carcasses at all
    pop <- pop[rowSums(pop[, 8:13, drop=FALSE]) > 0,]
    pop <- pop[, 1:13]

    return(list(pop, Incidence, BB, "Eep"=sum(Eep), "Sdpb"=sum(Sdpb), "Sdpd"=sum(Sdpd), "Iep"=sum(Iep), "Edp"=sum(Edpd), "Rep"=sum(Rep), "Cep"=sum(Cep), "Rdpd"=sum(Rdpd), "Ccd"=sum(Ccd), "Zcd"=sum(Zcd)))

}
