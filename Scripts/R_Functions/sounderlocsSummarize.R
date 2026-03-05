
sounderlocsSummarize <- function(sounderlocs, r){
    #default will be just basic SEIRCZ summary per timestep
    SEIRCZ.only <- sounderlocs[, c('time', 'S', 'E', 'I', 'R', 'C', 'Z')] # to be compatible with homogeneous landscapes, have to choose by name b/c they skip a column in centroids
    SEIRCZ.rep <- SEIRCZ.only %>%
    dplyr::group_by(time) %>%
    dplyr::summarise(S=sum(S), E=sum(E), I=sum(I), R=sum(R), C=sum(C), Z=sum(Z)) %>% as.data.frame()
    SEIRCZ.rep$rep <- r
    SEIRCZ_total <- SEIRCZ.rep[, c(8, 1:7)] # want rep in front
    list.out <- list("SEIRCZ_total" = SEIRCZ_total)
    return(list.out)
  
}


