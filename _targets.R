# _targets.R

#Friday May 16, get parameters file set up to be able to loop
    #create function to make matrices with a row for each parameter setting
        #inputs: parameters, list of inputs, structured so that each unit is its own setting
            #check that list of inputs, names match parms with 'input' label.. otherwise do warning/stop
            #create matrix with all combinations of parameters
            #output variable parameter matrix

# Targets setup --------------------
setwd(this.path::this.dir())

#load libraries for targets script
#install.packages("geotargets", repos = c("https://njtierney.r-universe.dev", "https://cran.r-project.org"))
library(targets)
library(tarchetypes)
library(geotargets)
# library(crew)

# This hardcodes the absolute path in _targets.yaml, so to make this more
# portable, we rewrite it every time this pipeline is run (and we don't track
# _targets.yaml with git)
tar_config_set(
    store = file.path(this.path::this.dir(),("_targets")),
    script = file.path(this.path::this.dir(),("_targets.R"))
)

#Source functions in pipeline
lapply(list.files(file.path("Scripts","R_Functions"), full.names = TRUE, recursive = TRUE), source)

#set options
options(clustermq.scheduler="multicore")

#Load packages
tar_option_set(packages = c("Rcpp",
                            "pracma",
                            "rdist",
                            "tidyverse",
                            "RcppArmadillo",
                            "RcppParallel",
                            "stringr",
                            "dplyr",
                            "sf",
                            "raster",
                            "terra",
                            "NLMR", ## NLMR barely still works, recommend using something else (CRAN friendly) for package dev.
                            "EnvStats",
                            "clustermq",
                            "deSolve",
                            "colorspace",
                            "data.table")#,
#                seed = 12345, ## can set the seed for reproducibility, or NA for non-reproducible totally stochastic -- see targets manual section 9.2
#                 error = 'stop') # for troubleshooting
)

# Pipeline ---------------------------------------------------------

list(
  
    ## Input raw data files -----

    ### Input parameters file: -----------
    tar_target(parameters_txt, file.path("Parameters.txt"), format="file"),

    ### Input landscapes directory: -----------
    tar_target(lands_path, file.path("Landscape_Setup","NND_Lands","4_Output", "sel_plands"), format="file"),
#     tar_target(lands_path, file.path("Input","lands"), format="file"),

    ## Get names of land tiles in order so we can associate lands with their data previously calculated
    tar_target(lands_names, data.table(land = seq(1, length(list.files(lands_path))), file = list.files(lands_path))),

    ## Read and format input data -----
    tar_terra_sprc(plands_sprc, ReadLands(lands_path)),

    ### Read and format parameters file: -----------
    tar_target(parameters0, FormatSetParameters(parameters_txt)),

    tar_target(variables, SetVarParms(parameters0)),

    ## Input cpp scripts as files to enable tracking -----
    tar_target(Fast_FOI_Matrix_script, file.path("Scripts","cpp_Functions","Fast_FOI_Matrix.cpp"), format="file"),
    tar_target(FindCellfromCentroid_script, file.path("Scripts","cpp_Functions","FindCellfromCentroid.cpp"), format="file"),
    tar_target(Movement_Fast_Generalized_script, file.path("Scripts","cpp_Functions","Movement_Fast_Generalized.cpp"), format="file"),
    tar_target(Movement_Fast_RSFavail_script, file.path("Scripts","cpp_Functions","Movement_Fast_RSFavail.cpp"), format="file"),
    tar_target(SpatialZones_fast_script, file.path("Scripts","cpp_Functions","SpatialZones_fast.cpp"), format="file"),
  
  ## Initialize surface -----
  ### Initialize grid(s): ---------------
    #Method for class 'SpatRaster'
        #Initialize_Grids(object)
    #Method for class 'SpatRasterCollection': 
        #Initialize_Grids(object)
    #Method for class 'numeric'
        #Initialize_Grids(object,grid.opts="homogenous")
    #Arguments
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

    ## makes the grid.opt parameters functional to choose lands variation... shifted grid.opts = "heterogeneous" to be randomized landscape (was "random" before, but unlisted)
    ## ugly, but functional:
    ## should probably move this to the InitializeGrids file
    tar_target(land_grid_list, {if (parameters0$pop_init_grid_opts == 'homogeneous'){
                                  if(parameters0$grid.opts != 'ras'){ # if grid.opts is homogeneous or heterogeneous
                                    # make a grid either uniform or random with even initial pig locations
                                    InitializeGrids(c(parameters0$len, parameters0$inc), parameters0$grid.opt)
                                  } else if (parameters0$grid.opts == 'ras'){ # if there is an input raster
                                    InitializeGrids(plands_sprc, parameters0$grid.opts)
                                  }
                                } else if (parameters0$pop_init_grid_opts == 'heterogeneous'){
                                    # make a grid with uneven pig initial locations...
                                    if (parameters0$grid.opts == 'homogeneous') {
                                      # can't do neutral plane with random pig distribution
                                      stop('Cannot run homogeneous grid.opts with heterogeneous pop_init_grid_opts')
                                    } else if (parameters0$grid.opts == 'heterogeneous'){
                                      # random pig distribution with random landscape
                                      InitializeGrids(c(parameters0$len, parameters0$inc), parameters0$grid.opt)
                                    } else if (parameters0$grid.opts == 'ras'){
                                      # random pig distribution with raster landscape
                                      InitializeGrids(plands_sprc, parameters0$grid.opts)
                                    }
                                }
                              }),

    ### Get surface parameters: ---------------
    tar_target(parameters, GetSurfaceParms(parameters0, plands_sprc[1])),

    ## Get landscape-specific movement parameters
    tar_target(mv.params, if(parameters$grid.opts=='ras'){ readRDS('./Landscape_Setup/NND_Lands/4_Output/ldsel.rds') } else { return(NA)}),

    ## Run Model ---------------
    tar_target(out.list,
        RunSimulationReplicates(land_grid_list = land_grid_list,
                                parameters = parameters,
                                variables = variables,
                                cpp_functions = list(Fast_FOI_Matrix_script, Movement_Fast_Generalized_script),
                                reps = parameters$nrep,
                                mv.parms = mv.params
        )
#         , cue = tar_cue(seed = FALSE) # allows existing simulation outputs to stand despite having stochastic elements, so long as inputs are the same
    ),
      ## Copy paste everything in the {} including the {} to run simulations using targets outputs without running targets so you can read the error messages and outputs! :)
      ## {lapply(list.files('./Scripts/R_Functions/', full.names=TRUE), source) ;RunSimulationReplicates(tar_read(land_grid_list), tar_read(parameters0), tar_read(variables), list(tar_read(Fast_FOI_Matrix_script), tar_read(Movement_Fast_Generalized_script)), tar_read(parameters)$nrep, tar_read(mv.params))}#, tar_read(burn.list)) }


    tar_target(plot_outputs, VisualOutputs(out.list, variables, land_grid_list, parameters, lands_names))
      ## Copy paste everything in the {} including the {} to run simulations using targets outputs without running targets so you can read the error messages and outputs! :)
      ## {lapply(list.files('./Scripts/R_Functions/', full.names=TRUE), source) ; VisualOutputs(tar_read(out.list), tar_read(variables), tar_read(land_grid_list), tar_read(parameters), tar_read(lands_names)) }


) # end targets list

## run tar_visnetwork() to visualize targets dependencies

