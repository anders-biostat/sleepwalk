`%||%` <- function(x, y) {
  if(is.null(x)) {
    y
  } else {
    x
  }
}

#' @export
sleepwalk <- function( embeddings, featureMatrices, maxdists = NULL, pointSize = 1.5, same = c( "objects", "features" ) ) {
  same = match.arg( same )
  
  #if there is only one embedding
  if(!is.null(dim(embeddings))) 
    embeddings <- list(embeddings)
  if(!is.null(dim(featureMatrices)))
    featureMatrices <- list(featureMatrices)    

  stopifnot( is.list(embeddings) )
  stopifnot( is.list(featureMatrices) )
  
  stopifnot( is.numeric(pointSize) && length(pointSize) == 1 )
    
  stopifnot( length(embeddings) <= 9 )
  
  #estimate maxdists from the data
  if(is.null(maxdists)) {
    maxdists <- sapply(1:length(featureMatrices), function(i) {
      message(paste0("Estimating 'maxdist' for feature matrix "), i)
      pairs <- cbind(sample(nrow(featureMatrices[[i]]), 1500, TRUE), 
                     sample(nrow(featureMatrices[[i]]), 1500, TRUE))
      quantile(sqrt(rowSums((featureMatrices[[i]][pairs[, 1], ] - featureMatrices[[i]][pairs[, 2], ])^2)), 0.5) 
    })
  }
  stopifnot( is.numeric(maxdists) )
    
  if( same == "objects" ) {
    stopifnot( length(embeddings) == length(featureMatrices) | length(featureMatrices) == 1 )
    stopifnot( length(maxdists) == length(featureMatrices) )
  }
  else {
    stopifnot( length(embeddings) == length(featureMatrices) )
    stopifnot( length(maxdists) == length(featureMatrices) | length(maxdists) == 1 )
  }
  
  oneFM <- NULL
  if(length(embeddings) != length(featureMatrices))
    oneFM <- 1
  
  for( i in 1:length(embeddings) ) {
     stopifnot( length( dim( embeddings[[i]] ) ) == 2 )
     stopifnot( length( dim( featureMatrices[[oneFM %||% i]] ) ) == 2 )
     stopifnot( ncol( embeddings[[i]] ) == 2 )
     stopifnot( nrow( embeddings[[i]] ) == nrow( featureMatrices[[oneFM %||% i]] ) )
     if( same == "objects" ) 
        stopifnot( nrow( embeddings[[i]] ) == nrow( embeddings[[1]] ) )
     else
       stopifnot( ncol( featureMatrices[[i]] ) == ncol( featureMatrices[[1]] ) )
  }

  JsRCom::openPage( FALSE, system.file( package="sleepwalk" ), "sleepwalk.html" )
  if( same == "objects" ) 
     JsRCom::sendData( "mode", "A" )
  else
     JsRCom::sendData( "mode", "B" )
  JsRCom::sendData( "n_charts", length(embeddings) )
  JsRCom::sendData( "maxdist", maxdists, TRUE )
  JsRCom::sendData( "embedding", embeddings, TRUE )
  JsRCom::sendData( "featureMatrix", featureMatrices, TRUE )
  JsRCom::sendData( "pointSize", pointSize )
  JsRCom::sendCommand( "set_up_chart()" )
}
