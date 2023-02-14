# Not %in% function
`%ni%` <- Negate('%in%')

# Find matches to any of a vector of character strings (for filtering)
greplany <- function(patterns, v) {
  match <- rep(FALSE, length(v))
  for (pattern in patterns) {
    match <- match | grepl(pattern, v)
  }
  return(match)
}

# Get base_dir
get_base_dir <- function() {
  if (Sys.info()["nodename"]=="Alexs-MacBook-Air-2.local") {
    "/Volumes/TracerX/working/VHL_GERMLINE/tidda/"
  } else {
    "/camp/project/tracerX/working/VHL_GERMLINE/tidda/"
  }
}

# Get cluster centroids from a Seurat object
get_centroids <- function(seu, reduction, ...) {
  dplyr::tibble(
    x = seu@reductions[[reduction]]@cell.embeddings[,1],
    y = seu@reductions[[reduction]]@cell.embeddings[,2],
    seu@meta.data
  ) %>%
    dplyr::group_by(...) %>%
    dplyr::summarise(x = median(x), y = median(y))
}
