#!/usr/bin/env Rscript 

# RUN AT COMMAND LINE
#   pyenv activate gee

#TODO: make this run using control files (study, env) instead of running once per study & env
#TODO: default to using the band name as the output column name of the env variable
#TODO: add optional parameter to not start task (for debugging)
#TODO: handle exception where point/env data can't be found
#TODO: right now, all datasets need to be in the same folder (datP)
# ==== Breezy setup ====

'
Annotates gee point datasets over multiple entities and variables.

Usage:
anno_gee.r <dat> <out> <study> <env_id> <colname> <band> 
anno_gee.r (-h | --help)

Control files:

Parameters:
  dat: folder containing feature collections to annotate
  out: path to output directory and file on gcs. do not include file extension, url info or bucket
  study: the focal study id
  env_id: gge id
  colname: what to call the var
  band: band form the gee image
  
Options:
-h --help     Show this screen.
-v --version     Show version.

' -> doc

#---- Input Parameters ----#
if(interactive()) {

  .pd <- '~/projects/dynamic_niches'
  .wd <- file.path(.pd)
  .seed <- NULL
  rd <- here::here
  
  #Required parameters
  #.envPF <- 'projects/map-of-life/diegoes/dist2road_USA_full' # 'NASA/ORNL/DAYMET_V4'
  .env_id <- 'NASA/MEASURES/GFCC/TC/v3'
  .datP <- 'users/syanco/cranes/tracks/4124813'
  #.outPF <- 'benc/projects/ms2/poc/mosey_env/ mosey_env1/anno/east_portugal_dist2urban'
  .outP <- 'annotated'
  #.outPF <- 'benc/projects/mosey_env/annotated/76367850_dist2water_month_test'
  .study <- 4124813
  #Optional parameters
  .band <- 0 #4
  
  #.colEnv <- 'dist2road' #  'tmax'
  .col_name <- 'pfor'
  .groups <- NULL #only run for these groups. handy for testing
  .npts <- NULL
} else {
  suppressWarnings(
    suppressPackageStartupMessages({
        library(docopt)
        library(rprojroot)
  }))

  ag <- docopt(doc, version = '0.1\n')

  .wd <- getwd()
  .script <-  suppressWarnings(thisfile())
  .seed <- ag$seed
  rd <- is_rstudio_project$make_fix_file(.script)
  
  source('src/funs/input_parse.r')
  
  #Required parameters
  #.envPF <- ag$env
  .datP <- ag$dat
  .outP <- ag$out
  .study <- ag$study
  .env_id <- ag$env_id
  .col_name <- ag$colname
  .band <- as.numeric(ag$band)
  # .band <- ag$band
  
  
    #Optional parameters
  #.band <- ifelse(is.null(ag$band),0,as.integer(ag$band))
  #.colEnv <- ifelse(is.null(ag$envcol),'env',ag$envcol)
  
  #ifelse can't return NULL, which is annoying
  if(is.null(ag$groups)) {.groups <- NULL } else {.groups <- as.integer(parseCSL(ag$groups))}
  if(is.null(ag$npts)) {.npts <- NULL} else {.npts <- as.integer(ag$npts)}
    
}

#---- Initialize Environment ----#

message("Initializing environment...")

if(!is.null(.seed)) {message(paste('Random seed set to',.seed)); set.seed(as.numeric(.seed))}

t0 <- Sys.time()

source('src/init/startup.r')

#Source all files in the auto load funs directory
list.files('src/funs/auto',full.names=TRUE) %>% walk(source)
#source(rd('src/main/dist2water_month.r'))

#For some reason I need to set these before I load rgee
# I didn't have to do this before
# UPDATE: after installing gcloud via brew, I don't have to do this anymore
#reticulate::use_python("/Users/benc/.pyenv/versions/gee/bin/python", required = TRUE)
#reticulate::use_virtualenv('/Users/benc/.pyenv/versions/gee',required=TRUE)

suppressWarnings(
  suppressPackageStartupMessages({
    library(rgee)
  }))

#Initialize gee
# Sys.setenv("RETICULATE_PYTHON" = "/home/syanco/mambaforge/envs/annotate/bin/python")
ee_check()
ee_Initialize()
  
#TODO: do a check to make sure rgee initialized correctly

#---- Local parameters ----#
.entity <- 'study'
#TODO: I might not need some of these variables anymore
#TODO: have a list called 'cols' instead of vars for each column name
.colImageId <- 'image_id';
.colMillis <- 'millis'
.colTimestamp <- 'timestamp'
.colGrp <- 'grp'
.bucket <- '1000cranes-bucket'
#.groupSize <- 500e3
#datN <- basename(.datPF)
#assetType <- ee$data$getAsset(.envPF)$type
#assetType <- 'IMAGE_COLLECTION'



#---- Control files ----#


entities <- read_csv(file.path(.wd,'ctfs',glue('{.entity}.csv'))) %>% filter(run==1)
# envs <- read_csv(file.path(.wd,'ctfs/env.csv')) %>% filter(run==1)

#---- Perform analysis ----#

# for(i in 1:nrow(entities)) {

  # i <- 1
  # entity <- entities[i,]
  
  message(glue('Processing {.datP}'))
  #.datPF <- 'users/benscarlson/projects/ms2/poc/mosey_env/mosey_env1/east_portugal'
  
  # pts <- ee$FeatureCollection(file.path(.datP,entity$study_id))
  pts <- ee$FeatureCollection(file.path(.datP))
  
  # #LOOP over envs
  # for(j in 1:nrow(envs)) {
  # 
  #   # j <- 1
  #   env <- envs[j,]
    
    message(glue('Setting up annotation tasks for {.col_name}'))
    
    #Check if the layer is a computed layer. If so load it.
    #Otherwise load gee asset
    if(file.exists(file.path("~", .env_id))) {
      source(file.path("~",.env_id))
      layer <- getLayer()
      assetType <- getAssetType()
    } else {
      assetType <- ee$data$getAsset(.env_id)$type
      
      if(assetType=='IMAGE') {
        layer <- ee$Image(.env_id)$select(list(.band))
      } else if (assetType=='IMAGE_COLLECTION') {
        layer <- ee$ImageCollection(.env_id)
      }  else {
        stop(glue('Invalid asset type: {assetType}'))
      }
    }
    
    #For testing
    # pts <- ee$FeatureCollection(pts$toList(10))
    # pts$aggregate_array(.colGrp)$getInfo()
    
    
    if(is.null(.groups)) {
      #Groups run from 0...n, so to get number of groups need to add 1
      maxgrp <- pts$aggregate_max(.colGrp)$getInfo()
      groups <- 0:maxgrp
    } else {
      groups <- .groups
      message(glue('Running only for the following group numbers: {groups}'))
    }
    
    if(length(groups)>1) message(glue('Splitting annotation into {length(groups)} tasks'))
    
    #====
    
    #Note groups start at 0
    for(group in groups) {
      
      #group <- 0
      ptsGrp <- pts$filter(ee$Filter$eq(.colGrp,group))
    
      #Make sure there are points in the group. Can result in 0 records if .group
      # does not exist in the dataset.
      invisible(assert_that(ptsGrp$size()$getInfo() > 0))
      
      if(!is.null(.npts)) {
        message(glue('Running for a subset of {.npts} points'))
        ptsGrp <- ee$FeatureCollection(ptsGrp$toList(.npts))
      }
    
      #Code for testing to reduce size in groups
      #ptsGrp <- ee$FeatureCollection(ee$Feature(ptsGrp$first()))
      #00030000000000056e8f does not return a value for dist2road. Note the property is not there after annotation.
      #ptsGrp <- ptsGrp$filter(ee$Filter$eq('system:index','00030000000000056e8f'))
      
    
      if(assetType=='IMAGE') {
        #Note that the band is selected above, when loading the layer
        anno <- layer$reduceRegions(
          reducer = ee$Reducer$median()$setOutputs(list(.col_name)),
          collection = ptsGrp,
          scale = layer$projection()$nominalScale()
        )
        
      } else if(assetType=='IMAGE_COLLECTION') {
        
        #env <- dist2water_month()
        #env <- ee$ImageCollection(.envPF)
        
        ptsGrp <- ptsGrp$map(function(f) {
          mil = ee$Date(f$get(.colTimestamp))$millis()
          f <- f$set(.colMillis,mil)
          return(f)
        })
        
        filter <- ee$Filter$And(
          ee$Filter$lessThanOrEquals('system:time_start', NULL, .colMillis),
          ee$Filter$greaterThan('system:time_end', NULL, .colMillis)
        )
        
        joined <- ee$Join$saveAll('features')$apply(layer, ptsGrp, filter)
        
        # imgx <- ee$Image(joined$first())
        # imgx$projection()$nominalScale()$getInfo()
        
        anno <- joined$map(function(img) {
          #img <- ee$Image(joined$first())
          img <- ee$Image(img)$select(list(.band))
          
          #View(img$projection()$getInfo())
          #img$projection()$nominalScale()$getInfo()
          
          fc <- ee$FeatureCollection(ee$List(img$get('features')))
          
          vals <- img$reduceRegions(
            #TODO: if I'm just extracting the pixel value, should I use
            # ee$Reducer$first() instead ?
            reducer=ee$Reducer$median()$setOutputs(list(.col_name)),
            #scale=30,
            scale=img$projection()$nominalScale(),
            collection=fc)
            #tileScale=2)
          
          vals <- vals$map(function(f) {
            f$set(.colImageId,img$get('system:index'))
          })
          
          return(vals)
        })$flatten()
        
      } else {
        stop(glue('Invalid asset type: {assetType}'))
      }
      
      anno <- anno$sort('anno_id')
      #View(anno$getInfo()); quit()
      #glue('{.outPF}_{group}') #OLD, remove
      
      fileN <- glue('{.study}_{.col_name}')
      
      task <- ee$batch$Export$table$toCloudStorage(
        collection=anno,
        description=fileN,
        bucket=.bucket,
        fileNamePrefix=file.path(.outP,fileN),
        fileFormat='csv',
        selectors=list('anno_id',.col_name)
      )
      
      task$start()
      message(glue('Task for group {group} started'))
      message(glue('Results will be saved to gs://{.bucket}/{.outP}/{fileN}.csv'))
    
    } #loop over groups

  # } #loop over variables
  
# } #loop over entities
#---- Finalize script ----#
message(glue('Script complete in {diffmin(t0)} minutes'))

#OLD CODE

#This method of assigning dynamic groups is way too slow.
#Assign groups if dataset size is greater than max group size
# if(pts$size()$getInfo() > .groupSize) {
#   #Get the sorted anno ids
#   #Convert to string so they can be dict keys
#   annoIdList <- ee$List(pts$aggregate_array('anno_id'))$
#     sort()$
#     #slice(0,10)$
#     map(ee_utils_pyfunc(function(val) {
#       return(ee$Number(val)$format())
#     }))
#   
#   #Get group number by dividing by group size and casting to into to truncate
#   groupList <- ee$List$sequence(0,annoIdList$size()$subtract(1))$
#     map(ee_utils_pyfunc(function(val) {
#       return(ee$Number(val)$divide(.groupSize)$int())
#     }))
#   
#   #Create dictionary mapping from anno_id to group
#   groupDict <- ee$Dictionary$fromLists(annoIdList,groupList)
#   
#   #Now assign the group to the feature collection
#   pts = pts$map(ee_utils_pyfunc(function(pt) {
#     annoid <- ee$Number(pt$get('anno_id'))$format()
#     return(pt$set('anno_grp',groupDict$get(annoid)))
#   }))
#   
# } else {
#   #If dataset size is less than groupSize, then assign everything to group 0
#   pts <- pts$map(ee_utils_pyfunc(function(pt) {return(pt.set('anno_grp',0))}))
# }