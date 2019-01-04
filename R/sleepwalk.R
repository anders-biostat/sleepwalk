`%||%` <- function(x, y) {
  if(is.null(x)) {
    y
  } else {
    x
  }
}

#' @export
sleepwalk <- function( embeddings, featureMatrices = NULL, maxdists = NULL, pointSize = 1.5, 
                       distances = NULL, same = c( "objects", "features" ) ) {
  same = match.arg( same )
  
  if(is.null(featureMatrices)) {
    if(same == "features")
      stop("In the `same features` mode feature matrices must be defined")
    if(is.null(distances))
      stop("One of the two arguments must be defined: 'featureMatrices', 'distances'")
    stopifnot(nrow(distances) == ncol(distances))
  }
  
  stopifnot( is.numeric(pointSize) && length(pointSize) == 1 )
  
  #if there is only one embedding
  if(!is.null(dim(embeddings))) 
    embeddings <- list(embeddings)
  stopifnot( is.list(embeddings) )
  
  if(!is.null(dim(featureMatrices)))
    featureMatrices <- list(featureMatrices)
  if(!is.null(dim(distances)))
    distances <- list(distances)
  
  stopifnot( is.list(featureMatrices %||% distances) )

  stopifnot( length(embeddings) <= 9 )

  if( same == "objects" ) {
    stopifnot( length(embeddings) == length(featureMatrices %||% distances) | 
                 length(featureMatrices %||% distances) == 1 )
  }
  else {
    stopifnot( length(embeddings) == length(featureMatrices) )
  }  
  
  oneFM <- NULL
  if(length(embeddings) != length(featureMatrices %||% distances))
    oneFM <- 1
  
  for( i in 1:length(embeddings) ) {
    stopifnot( length( dim( embeddings[[i]] ) ) == 2 )
    stopifnot( length( dim( (featureMatrices %||% distances)[[oneFM %||% i]] ) ) == 2 )
    stopifnot( ncol( embeddings[[i]] ) == 2 )
    stopifnot( nrow( embeddings[[i]] ) == nrow( (featureMatrices %||% distances)[[oneFM %||% i]] ) )
    if( same == "objects" ) 
      stopifnot( nrow( embeddings[[i]] ) == nrow( embeddings[[1]] ) )
    else
      stopifnot( ncol( featureMatrices[[i]] ) == ncol( featureMatrices[[1]] ) )
  }
      
  #estimate maxdists from the data
  if(is.null(maxdists)) {
    if(!is.null(featureMatrices)) {
      maxdists <- sapply(1:length(featureMatrices), function(i) {
        message(paste0("Estimating 'maxdist' for feature matrix "), i)
        pairs <- cbind(sample(nrow(featureMatrices[[i]]), 1500, TRUE), 
                       sample(nrow(featureMatrices[[i]]), 1500, TRUE))
        median(sqrt(rowSums((featureMatrices[[i]][pairs[, 1], ] - featureMatrices[[i]][pairs[, 2], ])^2))) 
      })
    } else {
      maxdists <- sapply(distances, median)
    }
  }
  stopifnot( is.numeric(maxdists) )
    
  if( same == "objects" ) {
    stopifnot( length(maxdists) == length(featureMatrices %||% distances) )
  }
  else {
    stopifnot( length(maxdists) == length(featureMatrices) | length(maxdists) == 1 )
  }
  
  JsRCom::openPage( FALSE, system.file( package="sleepwalk" ), "sleepwalk.html" )
  if( same == "objects" ) 
     JsRCom::sendData( "mode", "A" )
  else
     JsRCom::sendData( "mode", "B" )
  
  JsRCom::sendData( "n_charts", length(embeddings) )
  JsRCom::sendData( "maxdist", maxdists, TRUE )
  JsRCom::sendData( "embedding", embeddings, TRUE )
  if(!is.null(featureMatrices)) {
    JsRCom::sendData( "featureMatrix", featureMatrices, TRUE )
  } else {
    JsRCom::sendData( "distance", distances, TRUE )
  }
  JsRCom::sendData( "pointSize", pointSize )
  JsRCom::sendCommand( "set_up_chart()" )
}
