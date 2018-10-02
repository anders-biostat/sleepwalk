sleepwalk <- function( embedding, featureMatrix, maxdist ) {
  stopifnot( length( dim( embedding ) ) == 2 )
  stopifnot( length( dim( featureMatrix ) ) == 2 )
  stopifnot( ncol( embedding ) == 2 )
  stopifnot( nrow( embedding ) == nrow( featureMatrix ) )
  stopifnot( is.numeric( maxdist ) && length( maxdist ) == 1 )

  sleepwalkMulti( list( embedding ), list( featureMatrix ), maxdist )
}

sleepwalkMulti <- function( embeddings, featureMatrices, maxdists ) {
  stopifnot( length(embeddings) <= 9 )
  stopifnot( length(embeddings) == length(featureMatrices) )
  stopifnot( is.list(embeddings) )
  stopifnot( is.list(featureMatrices) )

  for( i in 1:length(embeddings) ) {
     stopifnot( length( dim( embeddings[[i]] ) ) == 2 )
     stopifnot( length( dim( featureMatrices[[i]] ) ) == 2 )
     stopifnot( ncol( embeddings[[i]] ) == 2 )
     stopifnot( nrow( embeddings[[i]] ) == nrow( featureMatrices[[i]] ) )
     stopifnot( nrow( embeddings[[i]] ) == nrow( embeddings[[1]] ) )
  }
      
  openPage( FALSE, system.file( package="sleepwalk" ), "sleepwalk.html" )
  sendData( "n_charts", length(embeddings) )
  sendData( "maxdist", maxdists )
  sendData( "embedding", embeddings )
  sendData( "featureMatrix", featureMatrices )
  sendCommand( "set_up_chart()" )
}
