
library(raster)
library(biomod2)
library(ggplot2)
library(dplyr)
library(tidyr)
library(rlang)
library(crayon)
library(readr)
library(lindia)
library(jtools)

## -------------------------------------------------------------------------------------- ##
## EXPECTED INPUTS ----
## -------------------------------------------------------------------------------------- ##

# Three main inputs are required to run this script:
#
# (i) A set of modeling outputs placed in the OUTPUTS/MODS folder and  generated
# By the SDMcalibrationForecasting-v1.R script.
#
# (ii) Total Annual productivity data (either single or multi-year) by region of interest 
# formatted as a csv table placed in the DATA/productivityData folder. This table should 
# contain the following fields/columns: 
#     (a) region names
#     (b) region unique ID codes (in integer values)
#     (c) year
#     (d) annual production 
#
# (iii) A GeoTIFF Raster file with region data. Each region is identified by an 
# integer value, the same used as the unique ID in the Total Annual productivity 
# table. This is necessary to enable joining the data across these sources. To 
# enable area calculations this should be in a projected coordinate system or  (check 
# below PROJ_RASTER_DATA and PROJ_COORD_SYSTEM parameter).
#


## -------------------------------------------------------------------------------------- ##
## INPUT PARAMETERS ----
## -------------------------------------------------------------------------------------- ##

# Annual productivity table as a CSV file
PRODUCTIVITY_TABLE_NAME <- "Olive_ProductionSurface_Andalusia_2010_2014.csv"

# Name of the field/column with annual/multi-annual productivity data. 
# This column should be inside the table located in the DATA/ProductivityData folder
RESP_VAR_FIELD_NAME <- "PRODUCTION_tons"

# Name of the field that contains the region unique identifier code 
# (an ineteger values). This code must be the same used in the region raster 
# file inside the DATA/Regions folder.
REGION_ID_FIELD_NAME <- "PID"

# Field/column name containing the region names
REGION_NAME_FIELD_NAME <- "PROVINCE"

# Is productivity data multi-annual? If TRUE the annual data will be firstly 
# aggregated using the average and then used as input for regression with SDM 
# suitable area as predictor
IS_MULTI_YEAR <- TRUE

# Name of the field/column with the reference year. 
YEAR_FIELD_NAME <- "YEAR"

# Raster band index to use from each species binary SDM. This represents 
# which combination of ensemble measure / threshold metric is used to construct 
# the species/variey raster stack. The file has a total of 9 bands which represent the 
# two-factor combination of: (1) average ensemble, (2) median ensemble, and (3) 
# the weighted average ensemble vs. (A) TSS threshold, (B) ROC threshold and, (C) 
# Kappa threshold. The 9 bands are in the following order: 1A, 1B, 1C, 2A, 2B, 2C, 
# 3A, 3B, 3C.
SDM_RST_BAND_INDEX <- 1

# Coordinate system of the data (must be the same for species records and 
# predictor variables in raster format). The projection name is a PROJ4-like string.
# The default is set to WGS 1984 geographic coordinate system
# Set as NULL if not needed or if data already has a defined coordinate system.
COORD_SYSTEM <- "+init=epsg:4326"

# Set this to  TRUE if the raster data needs to be re-projected to a different 
# coordinate reference system as to allow area calculations. This may be required if 
# the original data is in a geographic coordinate system. 
# Set to FALSE if data does not require a re-projection to a projected system.
PROJ_RASTER_DATA <- TRUE

# A CRS PROJ4 string that will be used to re-project the raster input data.
PROJ_COORD_SYSTEM <- "+init=epsg:32630"

# Spatial resolution (pixel size) to attribute to the re-projected raster 
# (usually in meters)
RST_SPATIAL_RES <- 1000

# Print progress messages? If TRUE messages and progress bars will be displayed.
VERBOSE <- TRUE

## -------------------------------------------------------------------------------------- ##
## Load input data and verify it ----
## -------------------------------------------------------------------------------------- ##

if(VERBOSE){
  message(green("Loading annual productivity data..."))
}

# Load annual productivity data
fpath <- paste("./DATA/ProductivityData/",PRODUCTIVITY_TABLE_NAME,sep="")

if(!file.exists(fpath)){
  stop("Unable to find the annual productivity data file! Please check input data")
}else{
  
  prodDataByRegion <- read_csv(fpath)
  
  # Perform multi-annual aggregation of productivity data?
  if(IS_MULTI_YEAR){
    prodDataByRegion <- prodDataByRegion %>% 
      group_by(across(.cols = c(all_of(REGION_NAME_FIELD_NAME),
                                all_of(REGION_ID_FIELD_NAME)))) %>% 
      summarise(across(.cols=all_of(RESP_VAR_FIELD_NAME), .fns=mean), 
                .groups="drop_last")
  }
}

if(VERBOSE){
 message(green("Done!\n")) 
}

## ------------------------------------------------------------------------ ##

if(VERBOSE){
  message(green("Loading region raster data..."))
}

# Load raster data for the regions
regionRstPath <- list.files("./DATA/Regions", pattern=".tif$|.TIF$", 
                            full.names = TRUE)

if(is.null(regionRstPath) || length(regionRstPath)==0){
  stop("Region raster does not exists! Please check input files")
}else{
  regionRst <- raster(regionRstPath)
  
  if(!is.null(COORD_SYSTEM)){
    crs(regionRst) <- crs(regionRst)
  }
}

if(PROJ_RASTER_DATA){
  
  regionRst <- projectRaster(regionRst, 
                             res = RST_SPATIAL_RES, 
                             crs = crs(PROJ_COORD_SYSTEM),
                             method = "ngb")
  
  writeRaster(regionRst, "./OUTPUTS/PROD_MODS/regions.tif", overwrite = TRUE)
  
}

if(VERBOSE){
  message(green("Done!\n")) 
}

## ------------------------------------------------------------------------ ##

if(VERBOSE){
  message(green("Gathering data files to build raster stacks..."))
}

# List sub-folders in MODS containing outputs generated by biomod2 package
spOutputPaths <- list.dirs("./OUTPUTS/MODS", recursive = FALSE, 
                           full.names = TRUE)

# Species/varieties names taken from MODS sub-folders
spNames <- basename(spOutputPaths)

i <- 0
for(spName in spNames){
  
  i <- i + 1 
  
  # Generate raster stack of binary predictions, one for each projection date
  targetFolder   <- spOutputPaths[i]
  targetProjDirs <- list.dirs(targetFolder, recursive = FALSE)
  targetProjDirs <- targetProjDirs[grepl("proj_",targetProjDirs)]
  #print(targetProjDirs)
  
  if(i==1){
    
    # Projection names in each species/variety folder
    projNames <- basename(targetProjDirs)
    
    # Initialize a list containing the SDM raster files by projection date 
    projFilesList <- list()
    for(k in 1:length(projNames)){
      projFilesList[[projNames[k]]] <- vector(mode="character", 
                                              length=length(spNames))
    }
  }
  
  j <- 0
  for(targetProjDir in targetProjDirs){
    j <- j + 1
    # Generate raster file names to build a raster stack per projection date
    projFilesList[[j]][i] <- paste(targetProjDir,"/GeoTIFF/",
                                   projNames[j],"_",spName,"_ensemble_TSSbin.tif",sep="")
  }
}


if(VERBOSE){
  message(green("Done!\n")) 
}


## ------------------------------------------------------------------------ ##

if(VERBOSE){
  message(green("Loading data and making calculations of suitable area 
                by projection date/scenario..."))
}

# Generate the raster stacks for each projection date (e.g., current, future2040, ...) 
# and holding all species/varieties. Then process it ti extract suitable area for all 
# projections

# Remove proj_ prefix from projection names
projNamesShort <- gsub("proj_","",projNames)

# Calculate pixel area for calculations
pixArea <- res(regionRst)[1]*res(regionRst)[2]


j <- 0
for(projName in projNames){
  
  j <- j + 1
  filesToRead <- projFilesList[[projName]]

  for(i in 1:length(filesToRead)){
    
    # Read files per species/variety
    rst <- stack(filesToRead[i])
    
    if(i==1){
      # Extract only the required band index for ensemble 
      # type/threshold metric combination
      tmpRstStack <- rst[[SDM_RST_BAND_INDEX]] 
    }else{
      # Stack all layers across each species/variety
      tmpRstStack <- stack(tmpRstStack, 
                           rst[[SDM_RST_BAND_INDEX]])
    }
  }
  
  # Set names to each raster layer in the stack corresponding 
  # to species/variety names
  names(tmpRstStack) <- spNames
  
  # Re-project raster data?
  if(PROJ_RASTER_DATA){
    tmpRstStack <- projectRaster(from   = tmpRstStack, 
                                 to     = regionRst, # Use regions as template
                                 method = "ngb")
  }
  
  # Create the raster stack by assignment
  assign(projNamesShort[j], tmpRstStack)
  
  # Write data to file
  writeRaster(tmpRstStack, paste("./OUTPUTS/PROD_MODS/",projNamesShort[j],"_SpBinEnvSuit.tif",sep=""), 
                                 overwrite = TRUE)
  
  
  # Calculate the total overlapped suitable area for all species/varieties
  tmpSum <- calc(tmpRstStack, fun = sum)
  tmpSum[tmpSum > 0] <- 1
  objName <- paste(projNamesShort[j],"all",sep="_")
  assign(objName, tmpSum)
  
  # Write data to file
  writeRaster(tmpSum, paste("./OUTPUTS/PROD_MODS/",projNamesShort[j],
                            "_AllSpBinEnvSuit.tif",sep=""), 
              overwrite = TRUE)
  
  # Load raster data and summarize the total area for each projection date/scenario
  tmpRstStackDF <- stack(regionRst, tmpSum) %>% 
    values() %>% 
    as.data.frame() %>% 
    na.omit() %>% 
    `colnames<-`(c("RID", "SDM")) %>% # RID - "internal" region ID / SDM - overlapped suitability (0,1)
    group_by(RID) %>% 
    summarise(totalArea = sum(SDM) * pixArea, .groups="drop") %>% 
    mutate(projDateScnName = projNamesShort[j])
  
  if(j==1){
    areaSuitableByRegionProj <- tmpRstStackDF
  }else{
    areaSuitableByRegionProj <- bind_rows(areaSuitableByRegionProj,
                                      tmpRstStackDF)
  }
  
  if(VERBOSE){
    message(yellow("-> Finished projection date/scenario:",projNamesShort[j]))
  }
  
}

## ------------------------------------------------------------------------ ##

if(VERBOSE){
  message(green("Building dataset and performing log-log regression...\n")) 
}

# Make data for regression including all projection dates/scenarios (current and future)
#
prodDF_cal <- areaSuitableByRegionProj %>% 
  filter(projDateScnName == "current") %>% 
  # Join annual productivity data for current conditions
  left_join(prodDataByRegion, 
            by = c("RID" = REGION_ID_FIELD_NAME)) %>% 
  # Join future data with region names/codes
  bind_rows(areaSuitableByRegionProj %>% 
              left_join(prodDataByRegion %>% select(!!REGION_ID_FIELD_NAME,
                                                    !!REGION_NAME_FIELD_NAME), 
                        by = c("RID" = REGION_ID_FIELD_NAME)) %>% 
              filter(projDateScnName != "current")) %>%
  # Add log10 transformed columns for log-log regression
  mutate(totalArea_log10 = log10(totalArea)) %>%
  mutate(prod_log10 = log10(.data[[RESP_VAR_FIELD_NAME]])) %>% 
  as.data.frame


# Do log-log regression for current conditions
lmod <- lm(prod_log10 ~ totalArea_log10, 
           data = prodDF_cal %>% filter(projDateScnName == "current"))


# Show and write model summary statistics

summ(lmod)

sink(file = "./OUTPUTS/PROD_MODS/modelSummary.txt")
summ(lmod)
sink()

summ(lmod, confint = TRUE)

sink(file = "./OUTPUTS/PROD_MODS/modelSummary_ConfIntervals.txt")
summ(lmod, confint = TRUE)
sink()

# Diagnostic plots for the log-log regression model
gd <- gg_diagnose(lmod, plot.all = TRUE)
ggsave("./OUTPUTS/PROD_MODS/RegressionDiagnosticsPlots.png", plot = gd, height = 10, width = 8)


## ------------------------------------------------------------------------ ##

if(VERBOSE){
  message(green("Doing model projections across all projection dates/scenarios...\n")) 
}

# Do projections based on the log-log model across all projection dates/scenarios
# and get confidence intervals
pred_DF <- predict(lmod, newdata = prodDF_cal, interval="confidence") %>% 
  as.data.frame %>% 
  bind_cols(prodDF_cal)

write_csv(pred_DF, path = "./OUTPUTS/PROD_MODS/model_predictions.csv", na = "")

g1 <- ggplot(pred_DF) + 
  geom_bar(aes(y=fit, x=projDateScnName), stat="identity") + 
  geom_errorbar(aes(y=fit, x=projDateScnName,ymin = lwr, ymax = upr),
                width = 0.2) + 
  facet_wrap(.data[[REGION_NAME_FIELD_NAME]]~.) + 
  labs(title="Trends in production across scenarios and regions") + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) + 
  xlab("Projection scenario/date") + 
  ylab("Projected productivity values (log scale)")
  
plot(g1)

ggsave("./OUTPUTS/PROD_MODS/ProductivityTrendByRegionScenario.png", plot = g1, height = 8, width = 8)

