#GetSurfaceParms(parameters00,plands_sprc[1])
GetSurfaceParms<-function(parameters, land){
    if (parameters$grid.opts == 'ras'){
        inc=res(land)[1]/1000 #get resolution in km
        km_len=dim(land)[1]*inc
        area=km_len^2

        #add to parameters list
#         parameters=append(parameters,inc)
#         names(parameters)[length(parameters)]<-"inc"
        parameters$inc <- inc

        parameters <- append(parameters,km_len)
        names(parameters)[length(parameters)] <- "km_len"
#         parameters$km_len <- km_len

        parameters <- append(parameters,area)
        names(parameters)[length(parameters)] <- "area"
#         parameters$area <- area
    } else {
        km_len <- parameters$inc * parameters$len
        area <- km_len^2
        parameters <- append(parameters, km_len)
        names(parameters)[length(parameters)] <- "km_len"
        parameters <- append(parameters, area)
        names(parameters)[length(parameters)] <- "area"
    }

	return(parameters)

}
