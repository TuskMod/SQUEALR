## split sounders if they are too large

sounderSplit <- function(pop, ss){
    # calculate probability of a sounder splitting, based on ss (sounder size) parameter
    # uses 1+ss so that 50% chance of split is when sounder is 1 higher than ss
    # may need to adjust steepness of function; look at StateChanges.R death probability function to see where parameter might fit
    split.prob <- exp(-(1+2*ss) + pop[,1])/(1+exp(-(1+2*ss) + pop[,1]))
    # if population gets too high too quickly, split.prob becomes NaN causing an error. In this case, probability of split becomes 1
    split.prob[is.nan(split.prob)] <- 1
    # sounders of size 2*ss or less will not split
    split.prob[pop[,1] <= 2*ss] <- 0
    # choose which sounders will split
    split.bin <- rbinom(split.prob, 1, split.prob)
    # split them roughly in half
    if(sum(split.bin) != 0) {
        split.size <- floor(pop[split.bin==1, 8:11, drop=FALSE]/2)
        # get a matrix of the new sounders
        new.sounders <- cbind(pop[split.bin==1, 1:7, drop=FALSE], split.size, C=0, Z=0)
        # subtract the departed sounders from the original sounders
        pop[split.bin==1,8:11] <- pop[split.bin==1,8:11, drop=FALSE] - split.size
        # add the new sounders to the modified original list
        pop <- rbind(pop, new.sounders)
        pop[,1] <- rowSums(pop[,8:11, drop=FALSE])
    }
    return(pop)
}
