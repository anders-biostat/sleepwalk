`%||%` <- function(x, y) {
  if(is.null(x)) {
    y
  } else {
    x
  }
}

#' Variable to store the information about the current \code{sleepwalk} session.
#' @export
.slw <- new.env()

#' Interactively explore one or several 2D embeddings
#' 
#' A function to interactively explore a 2D embedding of some higher-dimensional
#' point cloud, as produced by a dimension reduction method such as MDS, t-SNE, or the like.
#'
#' The function opens a browser window and displays the embeddings as point clouds. When the user
#' moves the mouse over a point, the point gets selected and all data points change colour such
#' that their colour indicates
#' the feature-space distance to the point under the mouse cursor. This allows to quickly and
#' intuitively check how tight clusters are, how faithful the embedding is, and how similar
#' the clusters are to each other.
#' 
#' 
#' @param embeddings either an \eqn{n x 2} embedding matrix (where \eqn{n} is a number of points) or
#' a list of \eqn{n_i x 2} matrices - one for each embedding. If \code{same = "objects"} all embedding
#' matrices must have the same number of rows.
#' @param featureMatrices either an \eqn{n x m} matrix of point coordinates in the feature-dimension
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
#' @param titles a vector of titles for each embedding. Must be the same length as the list of 
#' \code{embeddings}.
#' @param distances distances (in feature space) between points that should be displayed as colours.
#' This is an alternative to \code{featureMatrices} if \code{same = "objects"}.
#' @param same defines what kind of distances to show; must be either \code{"objects"} or \code{"features"}.
#' Use \code{same = "objects"} when all the embeddings show the same set of points. In this case,
#' each embedding is colored to show the distance of the selected point to all other points.
#' The same or different features can be supplied as \code{featureMatrices}, to use the same or different distances
#' in the different embeddings.
#' \code{same = "features"}
#' is used to compare different sets of points (e.g. samples from different patients, or different batches) 
#' in the same feature space. In this case the distance is calculated from the selected point to all other 
#' points (including those in other embeddings).
#' @param compare defines what kind of comparison to perform; must be either \code{"embeddings"} or
#' \code{"distances"}. If \code{compare == "embeddings"}, then in each of the displayed embeddings
#' all the points will be coloured the same way, even if different feature or distance matrices are provided. 
#' This allows one to immediately identify the corresponding points and neighbourhoods in each of the embeddings.
#' If \code{commpare == "distances"}, point colours for each embedding are calculated independently. This allows, 
#' for instance, to compare different metrics or show an additional layer of information, when exploring an
#' embedding. This parameter has no effect if \code{same == "features"}.
#' @param ncol number of columns in the table, where all the embeddings are placed.
#' @param nrow number of rows in the table, where all the embeddings are placed.
#' @param saveToFile path to the .html file where to save the plots. The resulting page will be fully interactive
#' and contain all the data. If this is \code{NULL}, than the plots will be shown as the web page in your 
#' default browser. Note, that if you try to save that page, using your browser's functionality,
#' it'll become static.
#' @param on_selection a callback function that is called every time the user selects a group of points in
#' the web browser. From the \code{sleepwalk} app it gets two arguments: The first one is a vector of indices of
#' all the selected points and the second one is an index of an embedding from where the points were selected.
#' @param mode defines whether to use Canvas or SVG to display points. Using Canvas is faster and allows to plot 
#' more points simultaneously, yet we currently consider SVG mode to be more stable and vigorously tested. In future
#' versions SVG mode will be deprecated. Must be one of \code{canvas} or \code{svg}.
#' 
#' @return None.
#' 
#' @author Simon Anders, Svetlana Ovchinnikova
#' 
#' @references \url{https://doi.org/10.1101/603589}
#'
#' @examples 
#' #generate cockscrew-shaped 3D data with 3 additional noisy dimensions
#' ts <- c(rnorm(100), rnorm(200, 5), rnorm(150, 13), runif(200, min = -5, max = 20))
#' 
#' a <- 3
#' w <- 1
#' 
#' points <- cbind(30*cos(w * ts), 30*sin(w * ts), a * ts)
#' 
#' ndim <- 6
#' noise <- cbind(matrix(rnorm(length(ts) * 3, sd = 5), ncol = 3),
#'                matrix(rnorm(length(ts) * (ndim - 3), sd = 10), ncol = ndim - 3))
#' 
#' data <- noise
#' data[, 1:3] <- data[, 1:3] + points
#' 
#' pca <- prcomp(data)
#' 
#' \donttest{#compare Euclidean distance with the real position on the helix
#' sleepwalk(list(pca$x[, 1:2], pca$x[, 1:2]), list(data, as.matrix(ts)), 
#'           compare = "distances", pointSize = 3)}
#' #the same, but with saving the web page to an HTML file
#' sleepwalk(list(pca$x[, 1:2], pca$x[, 1:2]), list(data, as.matrix(ts)), 
#'           compare = "distances", pointSize = 3,
#'           saveToFile = paste0(tempdir(), "/test.html"))
#' 
#'   
#' @importFrom jsonlite toJSON
#' @importFrom stats median
#' @import jrc
#' @export
sleepwalk <- function( embeddings, featureMatrices = NULL, maxdists = NULL, pointSize = 1.5, titles = NULL,
                       distances = NULL, same = c( "objects", "features" ), compare = c("embeddings", "distances"),
                       saveToFile = NULL, ncol = NULL, nrow = NULL, on_selection = NULL, mode = c("canvas", "svg")) {
  same = match.arg( same )
  compare = match.arg( compare )
  mode = match.arg( mode )
  
  if(is.null(featureMatrices)) {
    if(same == "features")
      stop("In the `same features` mode feature matrices must be defined")
    if(is.null(distances))
      stop("One of the two arguments must be defined: 'featureMatrices', 'distances'")
    stopifnot(nrow(distances) == ncol(distances))
  }
  
  stopifnot( is.numeric(pointSize) && length(pointSize) == 1 )
  
  rm(list = ls(envir = .slw), envir = .slw)
  if(is.null(on_selection)) {
    .slw$on_selection <- function(points, emb) {
      if(length(points) > 0){
        message(paste0("You've selected ", length(points), " points from the embedding ", emb, "."))
        message(paste0("The indices of the selected points are now stored in the variable 'selPoints'."))
      }
    }
  } else {
    stopifnot(is.function(on_selection))
    
    n_args <- length(formals(on_selection))
    stopifnot(n_args < 3)
    
    .slw$on_selection <- function(points, emb) {
      if(n_args == 0) on_selection()
      if(n_args == 1) on_selection(points)
      if(n_args == 2) on_selection(points, emb)
    }
      on_selection
  }
  
  #if there is only one embedding
  if(!is.null(dim(embeddings))) 
    embeddings <- list(embeddings)
  stopifnot( is.list(embeddings) )
  
  if(length(embeddings) == 1 && same == "features")
    same = "objects"
  
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
  
  if(!is.null(titles)) {
    stopifnot(length(titles) == length(embeddings))  
  } else {
    if(!is.null(names(embeddings))) {
      titles <- names(embeddings)
    } else {
      titles <- rep("", length(embeddings))
    }
  }
  
  if(!is.null(names(embeddings)))
    embeddings <- unname(embeddings)
      
  #estimate maxdists from the data
  if(is.null(maxdists)) {
    if(!is.null(featureMatrices)) {
      maxdists <- sapply(1:length(featureMatrices), function(i) {
        message(paste0("Estimating 'maxdist' for feature matrix "), i)
        pairs <- cbind(sample(nrow(featureMatrices[[i]]), 1500, TRUE), 
                       sample(nrow(featureMatrices[[i]]), 1500, TRUE))
        median(sqrt(rowSums((featureMatrices[[i]][pairs[, 1], , drop = FALSE] - featureMatrices[[i]][pairs[, 2], , drop = FALSE])^2))) 
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
    jrc::openPage( FALSE, system.file( package="sleepwalk" ), paste0("sleepwalk_", mode, ".html") )
    
    if( same == "objects" ) 
      jrc::sendData( "mode", "A" )
    else
      jrc::sendData( "mode", "B" )
    
    jrc::sendData( "n_charts", length(embeddings) )
    jrc::sendData( "titles", titles, TRUE )
    jrc::sendData( "maxdist", maxdists, TRUE )
    jrc::sendData( "embedding", embeddings, TRUE )
    if(!is.null(featureMatrices)) {
      jrc::sendData( "featureMatrix", featureMatrices, TRUE )
    } else {
      jrc::sendData( "distance", distances, TRUE )
    }
    jrc::sendData( "pointSize", pointSize )
    if(!is.null(ncol))
      jrc::sendData( "ncol", ncol )
    if(!is.null(nrow))
      jrc::sendData( "nrow", nrow )
    jrc::sendData( "compare", compare )
    jrc::sendCommand( "set_up_chart()" )
  } else {
    content <- readLines(paste0(system.file( package="sleepwalk" ), "/",  paste0("sleepwalk_", mode, ".html")), warn = F)
    
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
      paste0("titles = ", toJSON(titles), ";"),
      ifelse(!is.null(nrow), 
        paste0("nrow = ", nrow, ";"), ""),
      ifelse(!is.null(ncol),
        paste0("ncol = ", ncol, ";"), ""),
      paste0("compare = '", compare, "';"),
      "set_up_chart();"
    )
    
    content <- c(content[1:(length(content) - 3)], newLines, "</script>", "</body>", "</html>")
    
    writeLines(content, saveToFile)
  }
}

#' Make a snapshot of the currently running Sleepwalk app
#' 
#' This function produces a static plot that shows a given state of the 
#' currently active Sleepwalk app. Double click on a point in the web browser
#' will generate a command that can reproduce this exact state of the app.
#' 
#' @param point an index of the focus point (i.e. the one, over which the mouse is hovering at the moment).
#' To learn the index of a point double click on it in the web browser.
#' @param emb an index of the embedding of the focus point. To learn the index of an embedding double click 
#' on any of its points in the web browser.
#' @param returnList if \code{TRUE} returns a list of \code{ggplot} objects (one per embedding), that can be easily 
#' modified later. Otherwise returns a single \code{ggplot} objects with all the embeddings.
#' 
#' @return a \code{ggplot} object or a list of \code{ggplot} objects.
#' @examples 
#' \donttest{data("iris")
#' sleepwalk(iris[, c(1,3)], iris[1:4], pointSize = 4)
#' slw_snapshot(10)}
#' 
#' @importFrom httpuv service
#' @import ggplot2
#' @importFrom scales squish
#' @importFrom cowplot plot_grid
#' @export
slw_snapshot <- function(point, emb = 1, returnList = FALSE) {
  stopifnot(is.numeric(point))
  stopifnot(is.numeric(emb))
  
  if(length(point) > 1) {
    warning("More than one focuse point is provided, only the first one will be used.")
    point <- point[1]
  }
  if(length(emb) > 1) {
    warning("More than one focuse embedding is provided, only the first one will be used.")
    emb <- emb[1]
  }
  
  en <- new.env()
  jrc::setEnvironment(en)
  en$finished <- 0
  jrc::sendCommand(paste0("getSnapshotData(", point - 1, ", ", emb - 1, ");"))

  for( i in 1:(10/0.05) ) {
    service(100)
    if( en$finished > 0 ) 
      break
    
    Sys.sleep( .05 )
  }
  
  jrc::setEnvironment(globalenv())
  
  if( en$finished == 0 )
    stop( "Failed to get embedding data from the server" )

  maxdists <- en$maxdists
  colours <- c("#000000", "#1A1935", "#15474E", "#2B6F39", "#767B33", "#C17A6F", "#D490C6", "#C3C0F2")
  if(is.list(en$embs)) {
    n_charts <- length(en$embs)    
  } else {
    n_charts <- dim(en$embs)[1]
  }
  
  if(n_charts == 1) {
    data <- as.data.frame(cbind(en$embs[1, , ], en$dists[1, ]))
    colnames(data) <- c("x1", "x2", "dists")
    ggplot() + geom_point(aes(x = data$x1, y = data$x2, colour = data$dists), size = en$pointSize/2) +
      scale_color_gradientn(colours = colours, limits = c(0, maxdists), oob = squish) +
      ggtitle(en$titles) +
      theme(axis.title = element_blank(), axis.line = element_blank(), panel.grid.major = element_blank(),
            axis.text = element_blank(), axis.ticks = element_blank(), panel.grid.minor = element_blank(),
            legend.position = "bottom", legend.title = element_blank()) + guides(colour = guide_colourbar(barwidth = 15, barheight = 0.5))
  } else {
    plots <- lapply(1:n_charts, function(i) {
      if(is.list(en$embs)) {
        data <- as.data.frame(en$embs[[i]])
      } else {
        data <- as.data.frame(en$embs[i, , ])
      }
      colnames(data) <- c("x1", "x2")
      if(is.list(en$dists)) {
        data$dists <- en$dists[[i]]
        md <- maxdists[i]
      } else {
        if(dim(en$dists)[1] == 1) {
          data$dists <- en$dists[1, ]
          md <- maxdists
        } else {
          data$dists <- en$dists[i, ]
          md <- maxdists[i]
        }
      }
      ggplot() + geom_point(aes(x = data$x1, y = data$x2, colour = data$dists), size = en$pointSize/2) +
        scale_color_gradientn(colours = colours, limits = c(0, md), oob = squish) +
        ggtitle(en$titles[i]) +
        theme(axis.title = element_blank(), axis.line = element_blank(), panel.grid.major = element_blank(),
              axis.text = element_blank(), axis.ticks = element_blank(), panel.grid.minor = element_blank(),
              legend.position = "bottom", legend.title = element_blank()) + guides(colour = guide_colourbar(barwidth = 15, barheight = 0.5))
    })
    if(returnList) {
      plots
    } else {
      plot_grid(plotlist = plots)
    }
  }
}

