# Calculates probability of occurring at a given distance (based on seed dispersal kernel; spatially balanced with a fat tail)

TwoDt <- function(dist, shape, scale){
    ((shape-1)/(pi * (scale^2)))*((1+((dist^2)/(scale^2)))^-shape)
}
