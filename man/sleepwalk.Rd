\name{sleepwalk}
\alias{sleepwalk}

\title{ interactively explore a 2D embedding }
\description{ A function to interactively explore a 2D embedding of some higher-dimensional
point cloud, as produced by a dimension reduction method such as MDS, t-SNE, or the like }
\usage{
sleepwalk(embedding, featureMatrix, maxdist)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{embedding}{ a matrix with one row per data point and two columns, for the x and y
     coordinates of the points in the embedding }
  \item{featureMatrix}{
    a matrix with one row per data point and one column for each feature-space dimension.
    The displayed distances will be calculated as Euclidean distances of the rows of this 
    matrix }
  \item{maxdist}{
    the maximum distance (in feature space) that should still be covered by the colour scale;
    higher distances are shown in light gray
  }
}
\details{
The function opens a browser window and displays the embedding as point cloud. When the user
moves the mouse over a point, all data points change colour such that their colour indicates
the feature-space distance to the point under the mouse cursor. This allows to quickly and
intuitively check how tight clusters are, how faithful the embedding is, and how similar
the clusters are. 
}
\value{
None.
}
\references{
}
\author{
Simon Anders
}

\seealso{
sleepwalkMulti
}
\examples{
}