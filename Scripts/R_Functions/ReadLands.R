ReadLands <- function(path){
    fs <- list.files(path,full.names=TRUE)
    plands_list <- vector(mode="list",length=length(fs))
    for(f in 1:length(fs)){
        plands_list[[f]] <- terra::rast(fs[f])
    }
    plands_sprc <- terra::sprc(plands_list)
    return(plands_sprc)
}
