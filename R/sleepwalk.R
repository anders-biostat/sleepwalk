sleepwalk <- function( embedding, featureMatrix, maxdist, pointSize = 1.5 ) {
  stopifnot( length( dim( embedding ) ) == 2 )
  stopifnot( length( dim( featureMatrix ) ) == 2 )
  stopifnot( ncol( embedding ) == 2 )
  stopifnot( nrow( embedding ) == nrow( featureMatrix ) )
  stopifnot( is.numeric( maxdist ) && length( maxdist ) == 1 )
  stopifnot( is.numeric(pointSize) && length(pointSize) == 1 )

  sleepwalkMulti( list( embedding ), list( featureMatrix ), maxdist, pointSize = pointSize )
}

sleepwalkMulti <- function( embeddings, featureMatrices, maxdists, pointSize = 1.5, same = c( "objects", "features" ) ) {
  same = match.arg( same )
  
  stopifnot( length(embeddings) <= 9 )
  stopifnot( length(embeddings) == length(featureMatrices) )
  if( same == "objects" ) 
     stopifnot( length(maxdists) == length(featureMatrices) )
  else
     stopifnot( length(maxdists) == 1 )
  stopifnot( is.list(embeddings) )
  stopifnot( is.list(featureMatrices) )
  stopifnot( is.numeric(pointSize) && length(pointSize) == 1 )
  
  for( i in 1:length(embeddings) ) {
     stopifnot( length( dim( embeddings[[i]] ) ) == 2 )
     stopifnot( length( dim( featureMatrices[[i]] ) ) == 2 )
     stopifnot( ncol( embeddings[[i]] ) == 2 )
     stopifnot( nrow( embeddings[[i]] ) == nrow( featureMatrices[[i]] ) )
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
