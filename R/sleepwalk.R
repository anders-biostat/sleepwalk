`%||%` <- function(x, y) {
  if(is.null(x)) {
    y
  } else {
    x
  }
}

#' Interactively explore one or several 2D embeddings
#' 
#' A function to interactively explore a 2D embedding of some higher-dimensional
#' point cloud, as produced by a dimension reduction method such as MDS, t-SNE, or the like.
#' 
#' @param embeddings either an \eqn{n x 2} embedding matrix (where \eqn{n} is a number of points) or
#' a list of \eqn{n_i x 2} matices - one for each embedding. If \code{same = "objects"} all embedding
#' matrices must have the same number of rows.
#' @param featureMatrices either an \eqn{n x m} matrix of points coordinates in the feature-dimension
#' space or a list of such matrices - one for each embedding. The displayed distances will be calculated 
#' as Euclidean distances of the rows of these matrices. Alternatively, if \code{same = "objects"}
#' it is possible to provide the distances directly via the \code{distances} argument. 
#' If \code{same = "features"} then all the points must be from the same feature space and therefore
#' have the same number of columns. It is possible to use one feature matrix for all the embeddings.
#' @param maxdists a vector of the maximum distances (in feature space) for each provided feature or
#' distance matrix that should still be covered by the colour 
#' scale; higher distances are shown in light gray. This values can be changed later interactively.
#' If not provided, maximum distances will be estimated automatically as median value of the 
#' distances.
#' @param pointSize size of the points on the plots.
#' @param distances distances (in feature spase) between points that should be displayed as colours.
#' This is an alternative to \code{featureMatrices} if \code{same = "objects"}.
#' @param same defines what kind of distances to show; must be one of \code{"objects", "features"}.
#' \code{same = "objects"} is used when all the embeddings show the same set of points. In this case,
#' the distance from the selected point on each of embedding to all other points of the same embedding
#' is shown. The same or different feature of distance matrix can be used for that. \code{same = "features"}
#' is used to compare different sets of points (e.g. samples from different patients, or different batches) 
#' in the same feature space. In this case the distance is calculated from the selected point to all other 
#' points (including those in other embeddings).
#' @param saveToFile path to the .html file where to save the plots. The resulting page will be fully interactive
#' and contain all the data. If this is \code{NULL}, than the plots will be shown as the web page in your 
#' default browser. Note, that if you try to save that page, using your browser's functionality,
#' it'll become static.
#' 
#' The function opens a browser window and displays the embeddings as point clouds. When the user
#' moves the mouse over a point, all data points change colour such that their colour indicates
#' the feature-space distance to the point under the mouse cursor. This allows to quickly and
#' intuitively check how tight clusters are, how faithful the embedding is, and how similar
#' the clusters are.
#' 
#' @return None.
#' 
#' @author Simon Anders, Svetlana Ovchinnikova
#' 
#' @importFrom jsonlite toJSON
#' @export
sleepwalk <- function( embeddings, featureMatrices = NULL, maxdists = NULL, pointSize = 1.5, 
                       distances = NULL, same = c( "objects", "features" ), saveToFile = NULL ) {
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
    embeddings[[i]] <- as.matrix(embeddings[[i]])
    if(!is.null(featureMatrices)) {
      featureMatrices[[oneFM %||% i]] <- as.matrix(featureMatrices[[oneFM %||% i]])
    } else {
      distances[[oneFM %||% i]] <- as.matrix(distances[[oneFM %||% i]])
    }
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
  
  if(is.null(saveToFile)) {
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
  } else {
    content <- readLines(paste0(system.file( package="sleepwalk" ), "/", "sleepwalk.html"), warn = F)
    
    while(sum(grepl("script src", content)) != 0) {
      i <- which(grepl("script src", content))[1]
      fName <- gsub("<script src=\"(.*?)\"></script>", "\\1", content[i])
      script <- readLines(paste0(system.file( package="sleepwalk" ), "/", fName), warn = F)
      content <- c(content[1:(i - 1)], "<script>", script, "</script>", content[(i + 1):length(content)])
    }
    
    newLines <- c(
      paste0("mode = ", ifelse(same == "objects", "'A'", "'B'"), ";"),
      paste0("n_charts = ", length(embeddings), ";"),
      paste0("maxdist = ", toJSON(maxdists), ";"),
      paste0("embedding = ", toJSON(embeddings), ";"),
      ifelse(!is.null(featureMatrices), 
        paste0("featureMatrix = ", toJSON(featureMatrices), ";"),
        paste0("distance = ", toJSON(distances), ";")),
      paste0("pointSize = ", pointSize, ";"),
      "set_up_chart();"
    )
    
    content <- c(content[1:(length(content) - 3)], newLines, "</script>", "</body>", "</html>")
    
    writeLines(content, saveToFile)
  }
}

#' On selection
#' 
#' This function is called each time any points are selected or deselected.
#' You can customise it by redefining.
#' 
#' @param points a vector of indices of the selected points.
#' @param emb an index of the embedding, where the points have been selected.
#' 
#' @export
slw_on_selection <- function(points, emb) {
  message(paste0("You've selected ", length(points), " points from the embedding ", emb, "."))
  message(paste0("The indices of the selected points are now stored in the variable 'selPoints'."))
  message(paste0("You can also redefine this function 'slw_on_selection' that is called each time any points are selected."))
  message(paste0("It's first argument is a vector of indices of all the selected points, and the second one is the index of ",
                 "the embedding, where they were selected."))
}
