.slw_get_assay_data <- function(object, assay, slot_name) {
  tryCatch(
    Seurat::GetAssayData(object, assay = assay, layer = slot_name),
    error = function(e) Seurat::GetAssayData(object, assay = assay, slot = slot_name)
  )
}

#' Interactively explore a Seurat object's embeddings
#'
#' A convenience wrapper around \code{\link{sleepwalk}} for \code{Seurat}
#' objects: pulls a 2D embedding and a feature matrix straight from the
#' object's dimensional reductions and assay data, instead of having to
#' extract them yourself.
#'
#' @param object A \code{Seurat} object.
#' @param embedding Name of one, or a vector of names of several,
#' dimensional reductions stored in \code{object} (e.g. \code{"umap"},
#' \code{"tsne"}) to use as the 2D embedding(s). Each must have exactly 2
#' dimensions; use \code{\link[Seurat]{Reductions}(object)} to see what's
#' available.
#' @param feature What to use as the feature space for computing the
#' displayed distances: either the name of a dimensional reduction stored
#' in \code{object} (e.g. \code{"pca"}, using all of its dimensions), or
#' the name of an assay data layer (\code{"data"}, \code{"counts"}, or
#' \code{"scale.data"}) to pull from \code{object}'s assay.
#' @param assay Which assay to pull \code{feature} from, if \code{feature}
#' names an assay data layer. Defaults to \code{\link[Seurat]{DefaultAssay}(object)}.
#' @param features Which features (genes) to use, if \code{feature} names
#' an assay data layer. Defaults to \code{\link[Seurat]{VariableFeatures}(object)}
#' if set, otherwise to all features in the assay - which can be slow and
#' memory-heavy for large datasets, so consider running
#' \code{\link[Seurat]{FindVariableFeatures}} first, or passing your own
#' gene list here.
#' @param titles Titles for each embedding; defaults to the \code{embedding}
#' name(s) themselves when more than one is given.
#' @param ... Further arguments passed on to \code{\link{sleepwalk}} (e.g.
#' \code{pointSize}, \code{on_selection}, \code{saveToFile}, \code{metric}).
#'
#' @return None. See \code{\link{sleepwalk}}.
#'
#' @author Simon Anders, Svetlana Ovchinnikova, Claude
#'
#' @examples
#' library(Seurat)
#' data("pbmc_small")
#'
#' # using the defaults (embedding = "umap", feature = "pca"): compares the
#' # UMAP embedding against distances in PCA space, a standard scRNA-seq QC
#' # check - does the embedding respect PCA-space neighborhoods?
#' pbmc_small <- RunUMAP(pbmc_small, dims = 1:15)
#' \dontrun{
#' sleepwalk_seurat(pbmc_small)
#'
#' # same idea, against the tSNE embedding instead
#' sleepwalk_seurat(pbmc_small, embedding = "tsne", feature = "pca", pointSize = 4)
#'
#' # colour by distances in (variable-feature) expression space instead
#' sleepwalk_seurat(pbmc_small, embedding = "tsne", feature = "data")
#' }
#'
#' @export
sleepwalk_seurat <- function(object, embedding = "umap", feature = "pca",
                              assay = NULL, features = NULL, titles = NULL, ...) {
  if (!requireNamespace("Seurat", quietly = TRUE))
    stop("sleepwalk_seurat() requires the 'Seurat' package; install it with install.packages(\"Seurat\").")

  stopifnot(is.character(embedding), length(embedding) >= 1)
  stopifnot(is.character(feature), length(feature) == 1)

  available_reductions <- Seurat::Reductions(object)

  get_embedding <- function(key) {
    if (!(key %in% available_reductions))
      stop(sprintf(
        "'%s' not found among object's reductions: %s",
        key, paste(available_reductions, collapse = ", ")
      ))
    emb <- Seurat::Embeddings(object, reduction = key)
    if (ncol(emb) != 2)
      stop(sprintf(
        paste(
          "Reduction '%s' has %d dimensions; sleepwalk needs a 2D embedding.",
          "Did you mean a different reduction (e.g. 'umap' or 'tsne'),",
          "or to use '%s' as 'feature' instead?"
        ),
        key, ncol(emb), key
      ))
    emb
  }

  embeddings <- lapply(embedding, get_embedding)
  names(embeddings) <- embedding

  if (feature %in% available_reductions) {
    featureMatrix <- Seurat::Embeddings(object, reduction = feature)
  } else if (feature %in% c("data", "counts", "scale.data")) {
    assay <- assay %||% Seurat::DefaultAssay(object)
    mat <- .slw_get_assay_data(object, assay, feature)
    features <- features %||% Seurat::VariableFeatures(object)
    if (length(features) == 0)
      features <- rownames(mat)
    featureMatrix <- t(as.matrix(mat[features, , drop = FALSE]))
  } else {
    stop(sprintf(
      "'%s' is neither one of object's reductions (%s) nor one of 'data'/'counts'/'scale.data'.",
      feature, paste(available_reductions, collapse = ", ")
    ))
  }

  if (is.null(titles) && length(embedding) > 1)
    titles <- embedding

  sleepwalk(
    if (length(embeddings) > 1) embeddings else embeddings[[1]],
    featureMatrices = featureMatrix,
    titles = titles,
    ...
  )
}
