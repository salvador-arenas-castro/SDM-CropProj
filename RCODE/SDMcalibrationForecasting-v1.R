
# Load R libraries
library(raster)
library(biomod2)
library(magrittr)

source("./RCODE/AncillaryFunctions-v1.R")

## -------------------------------------------------------------------------------------- ##
## Parameters and initialization ----
## -------------------------------------------------------------------------------------- ##

# NOTE: THESE PARAMETERS MUST BE SET PRIOR TO RUNNING MODELS

# Input table with x,y species records located in DATA/ResponseVariable
SOURCE_INPUT_TABLE <- "OLIVE_VARIETIES_DATA_all.csv"

# Raster file format for predicotr data. This will be used to scan the 
# DATA/PredictorVariables folder and build one raster stack object for each 
# projection scenario/date
RASTER_FILE_FORMAT <- "asc"

# Output folder path (default points to OUTPUTS sub-folder)
# Do not use relative paths otherwise it will cause biomod2 to fail
# when looping through each species
OUTPUT_FOLDER <- "C:/MyFiles/R-dev/SDM-CropProj/OUTPUTS/MODS"

# Species names vector. If NULL then names will be extracted from 
#the NAME column in the table
SPECIES_NAMES <- NULL

# Coordinate system of the data (must be the same for species records and 
# predictor variables in raster format). The projection name is a PROJ4-like string.
# The default is set to WGS 1984 geographic coordinate system
# Set as NULL if not needed or if data already has a defined coordinate system.
COORD_SYSTEM <- "+init=epsg:4326"

# Number of pseudo-absence sets. The default is 1 for this example but it 
# must be set higher (e.g., 3,5,10...).
NR_PSEUDOABS_SETS <- 1

# Number of random pseudo-absences by set that will be generated
NR_PA_PER_SET <- 100

# biomod2's model hyperparameters
# By default GAM's smooth degree is changed to k=4 to avoid overly complex models
# Check ?BIOMOD_ModelingOptions for more details
BIOMOD_MODEL_OPTS <- BIOMOD_ModelingOptions(GAM = list(k = 4))

# A character defining models to  run. Check ?BIOMOD_Modeling for details
MODELS_TO_RUN <- c('GLM','GBM','GAM','CTA','ANN',
                   'FDA','MARS','RF','MAXENT.Phillips.2')

# Number of evaluation rounds for holdout cross-validation. Default is 10 
# for testing purposes only but should be slightly higher (e.g., 20, 30, ...). 
# Be aware that very large values will make model calibration and projection 
# very time-consuming
NR_EVAL_ROUNDS <- 10

# Percentage of data used for training the algorithms
PERC_TRAIN <- 80

# Species prevalence. Either NULL (default) or a 0-1 numeric used to build 
# 'weighted response weights'
PREVALENCE <- 0.5

# Number of permutation to estimate variable importance. A value of zero means 
# that it is not estimated. Set higher (10, 20) to calculate variable 
# importance values.
NR_VARIMP_ROUNDS <- 5

# Evaluation metric name used to subset partial models for building the 
# ensemble. Default is "TSS" the true-skill statistic. Other options are 
# "ROC" or "KAPPA"
EVAL_METRIC_NAME <- "TSS"

# A value between 0-1 typically closer to 1 (e.g. 0.9) as to defined the 
# top-best quantile of best models. Instead of defining a value for selecting 
# the models to be combined. 
# This option assesses the distribution of evaluation metrics in 
# EVAL_METRIC_NAME. If NULL this will not be used. Default is 0.8 meaning 
# that the top 20% models will be used in the ensemble.
EVAL_QUANTILE_THRESH <- 0.8

# Threshold value for the selected evaluation metric used to select 
# partial models for building the ensemble. Default: 0.75.
# This will be override by EVAL_QUANTILE_THRESH if that parameter is 
# not null.
EVAL_METRIC_THRESH <- 0.75


## -------------------------------------------------------------------------------------- ##
## Read the presence data ----
## -------------------------------------------------------------------------------------- ##

DataSpecies <- read.csv(paste("./DATA/ResponseVariable/",SOURCE_INPUT_TABLE,sep=""),
                        sep=",", stringsAsFactors = FALSE)

# Names for the species
if(is.null(SPECIES_NAMES)){
  spNames <- unique(DataSpecies$NAME)
}else{
  spNames <- SPECIES_NAMES
}


## -------------------------------------------------------------------------------------- ##
## Load environmental data by date/scenario ----
## -------------------------------------------------------------------------------------- ##

# Load environmental predictors from DATA/PredictorVariables/Current folder
current    <- stack(list.files(path="./DATA/PredictorVariables/Current", 
                               pattern='.asc$',full.names=TRUE))               
crs(current) <- CRS(COORD_SYSTEM)

# List the directories holding environmental data used for future projections
projDirs <- list.dirs("./DATA/PredictorVariables/Future", recursive = FALSE)

# Loop through each directory in projDirs and create a RasterStack object
for(projDir in projDirs){
  
  projName <- paste("future", basename(projDir), sep = "")
  
  tmpFileList <- list.files(path    = projDir, 
                            pattern = paste(".",RASTER_FILE_FORMAT, "$", sep = ""), 
                            full.names = TRUE)
  
  tmpRstStack <- stack(tmpFileList)
  
  if(!is.null(COORD_SYSTEM)){
    crs(tmpRstStack) <- CRS(COORD_SYSTEM)
  }
  
  assign(projName, tmpRstStack)
  
}

# Names for the 'projective' raster stacks
projNames <- c("current", paste("future",basename(projDirs),sep=""))


## -------------------------------------------------------------------------------------- ##
## Calibrate models by species ----
## -------------------------------------------------------------------------------------- ##

for(spName in spNames){
  
  # Change the working directory to the output path of the repo
  setwd(OUTPUT_FOLDER)
  
  # Subset the data and create a SpatialPoints object with input data
  myRespCoord <- DataSpecies[DataSpecies$NAME == spName, c("X","Y")]
  
  # Presence records as a SpatialPoints object for each species
  spPointData <- SpatialPoints(myRespCoord, proj4string=CRS(COORD_SYSTEM))
  
  # Prepare data for running biomod2
  myBiomodData <- BIOMOD_FormatingData(resp.var       = spPointData,
                                       expl.var       = current,
                                       resp.name      = spName,
                                       PA.nb.rep      = NR_PSEUDOABS_SETS,
                                       PA.nb.absences = NR_PA_PER_SET,
                                       PA.strategy    = 'random')
  
  # Variable importance # of permutations = 0 (no evaluation)
  # Evaluation metrics = ROC, TSS, KAPPA
  # No rescaling of model predictions !!
  # Full models = TRUE

  myBiomodModelOut <- BIOMOD_Modeling(
    myBiomodData,
    models         = MODELS_TO_RUN,
    models.options = BIOMOD_MODEL_OPTS,
    NbRunEval      = NR_EVAL_ROUNDS,
    DataSplit      = PERC_TRAIN,
    Prevalence     = PREVALENCE,
    VarImport      = NR_VARIMP_ROUNDS,
    models.eval.meth  = c('TSS','ROC','KAPPA'),
    SaveObj           = TRUE,
    rescal.all.models = FALSE,
    do.full.models    = TRUE)
  
  
  # Get model evaluation values
  myBiomodModelEval <- get_evaluations(myBiomodModelOut)
  
  # Print ROC scores
  print(myBiomodModelEval["ROC","Testing.data",,,])
  
  # Get boxplot stats
  print(fivenum(as.numeric(myBiomodModelEval["ROC","Testing.data",,,])))
  
  # Save evaluation metrics from the arrays
  evalDF.ROC <- as.data.frame(myBiomodModelEval["ROC","Testing.data",,,])
  evalDF.TSS <- as.data.frame(myBiomodModelEval["TSS","Testing.data",,,])
  evalDF.KAPPA <- as.data.frame(myBiomodModelEval["KAPPA","Testing.data",,,])
  
  write.csv(evalDF.ROC, file = paste(getwd(),"/",spName,"/",
                                     spName,"_evalDF_ROC.csv",sep=""))
  write.csv(evalDF.TSS, file = paste(getwd(),"/",spName,"/",
                                     spName,"_evalDF_TSS.csv",sep=""))
  write.csv(evalDF.KAPPA, file = paste(getwd(),"/",spName,"/",
                                       spName,"_evalDF_KAPPA.csv",sep=""))
  
  # Calculate variable importance across all PA sets, eval rouns and algorithms 
  varImportance <- get_variables_importance(myBiomodModelOut)
  varImportanceByVariableAVG <- apply(varImportance,1,mean, na.rm=TRUE)
  varImportanceByVariableSTD <- apply(varImportance,1,sd, na.rm=TRUE)
  vimpDF <- data.frame(cnames    = names(varImportanceByVariableAVG),
                       vimpAVG   = varImportanceByVariableAVG, 
                       varImpSTD = varImportanceByVariableSTD) %>% 
    arrange(desc(vimpAVG))
  
  print(vimpDF)
  write.csv(vimpDF, file = paste(getwd(),"/",spName,"/",spName,"_varImportance.csv",sep=""))
  
  ## -------------------------------------------------------------------------------------- ##
  ## Perform ensemble modelling ----
  ## -------------------------------------------------------------------------------------- ##
  
  if(!is.null(EVAL_QUANTILE_THRESH)){
    if(EVAL_METRIC_NAME == "ROC"){
      modSelCutoff <- quantile(evalDF.ROC, probs=EVAL_QUANTILE_THRESH, na.rm=TRUE)
    }else if(EVAL_METRIC_NAME == "KAPPA"){
      modSelCutoff <- quantile(evalDF.KAPPA, probs=EVAL_QUANTILE_THRESH, na.rm=TRUE)
    }else if(EVAL_METRIC_NAME == "TSS"){
      modSelCutoff <- quantile(evalDF.TSS, probs=EVAL_QUANTILE_THRESH, na.rm=TRUE)
    }else{
      stop("Unknown option in EVAL_METRIC_THRESH! Please check this.")
    }
  }else{
    modSelCutoff <- EVAL_METRIC_THRESH
  }
  
  # Adjust threshold value if the quantile is equal to 1
  if(modSelCutoff == 1){
    modSelCutoff <- 0.99
  }
  
  # Perform the ensemble modeling step by selecting target models and define how to combine them
  myBiomodEM <- BIOMOD_EnsembleModeling(modeling.output               = myBiomodModelOut,
                                        chosen.models                 = 'all',
                                        em.by                         = 'all',
                                        eval.metric                   = EVAL_METRIC_NAME,
                                        eval.metric.quality.threshold = modSelCutoff,
                                        prob.mean                     = TRUE,
                                        prob.cv                       = FALSE,
                                        prob.ci                       = FALSE,
                                        prob.ci.alpha                 = 0.05,
                                        prob.median                   = TRUE,
                                        committee.averaging           = FALSE,
                                        prob.mean.weight              = TRUE,
                                        prob.mean.weight.decay        = 'proportional')
  
  
  # Get evaluation scores for the Ensemble Modelling stage
  emEvalDF <- as.data.frame(get_evaluations(myBiomodEM))
  write.csv(emEvalDF, file = paste(getwd(),"/",spName,"/",
                                   spName,"_EnsMod_evalDF_AllMetrics.csv",sep=""))
  
  
  ## -------------------------------------------------------------------------------------- ##
  ## Obtain spatiotemporal projections ----
  ## -------------------------------------------------------------------------------------- ##
  
  # Models to consider in the ensemble and projection
  modelsToUse <- get_kept_models(myBiomodEM, 1)
  
  
  for(projName in projNames){
    
    # Obtain spatiotemporal projections
    myBiomodProj <- BIOMOD_Projection(modeling.output = myBiomodModelOut,
                                      new.env         = get(projName),
                                      proj.name       = projName, ## Name of the projection from above variable proj.name
                                      selected.models = modelsToUse,
                                      filtered.meth   = NULL,
                                      binary.meth     = NULL,
                                      compress        = 'gzip',
                                      clamping.mask   = TRUE,
                                      output.format   = '.grd',
                                      do.stack        = TRUE)
    
    
    # Perform the ensembling of projections
    myBiomodEF <- BIOMOD_EnsembleForecasting(projection.output = myBiomodProj,
                                             binary.meth       = c('ROC','TSS','KAPPA'),
                                             EM.output         = myBiomodEM,
                                             output.format     = '.grd')
    
    # Convert all output raster files to GeoTIFF
    inFolder <- paste(getwd(),"/",spName,"/proj_",projName,sep="")
    outFolder <- paste(inFolder,"/","GeoTIFF", sep="")
    dir.create(outFolder)
    
    # Convert raster data in native raster package to GeoTIFF format
    convertToGeoTIFF(inFolder, outFolder)
    
  }
  
  save.image(file=paste(spName,"ModObjects.RData",sep="_"))
  
}


