#' Based on Dynamic World v1
#' Data available 2015-06-27 - present
#' Implementation of a computed layer
#' 
#' Should have the following two functions:
#'  getAssetType() Returns IMAGE or IMAGE_COLLECTION
#'  getLayer() Returns an image or image collection representing the layer
#'    Note that images from the layer need to return the correct espg and scale


# # for testing:
# counties <- ee$FeatureCollection("TIGER/2016/Counties")
# filtered <- counties$filter(ee$Filter$eq("NAMELSAD", "Dane COunty"))
# geometry <- filtered$geometry()

getAssetType <- function() {return('IMAGE_COLLECTION')}

getLayer <- function() {
  
  dw <- ee$ImageCollection("GOOGLE/DYNAMICWORLD/V1")
  
  propIC <- dw$map(function(img){
    
    # img <- dw$filterDate("2020-01-01", "2020-01-01")$
    #   filterBounds(geometry)
    
    # mode composit of habitat label
    classification <- img$select('label')
    dwComposite <- classification$reduce(ee$Reducer$mode())
    
    # Binary classifier (trees or not trees)
    trees <- dwComposite$eq(1)
    
    # perform the proportion calculation
    
    # get total pixels in the radius
    tot <- trees$reduceNeighborhood(
      reducer = ee$Reducer$count(),
      kernel = ee$Kernel$circle(
        radius = 1000, 
        units = 'meters', 
        normalize = FALSE
      )
    )
    
    # get only tree pixels
    treesMasked <- trees$selfMask()
    treesCount <- treesMasked$reduceNeighborhood(
      reducer = ee$Reducer$count(),
      kernel = ee$Kernel$circle(
        radius = 1000, 
        units = 'meters', 
        normalize = FALSE
      )
    )
 
    # calculate the proportion
    prop <- treesCount$divide(tot);
    
    # add the image properties
    prop <- prop$
      copyProperties(img)$
      set('system:time_end',
          ee$Date(img$get('system:time_start'))$
            advance(1,'day')$
            millis())$
      set('system:time_start',
          ee$Date(img$get('system:time_start'))$
            advance(1, 'day')$
            millis())
    
    return(prop)
  })  
  
  return(propIC)
}

#View(distIC$first()$getInfo())