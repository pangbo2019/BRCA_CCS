# 2. Jaccard correlation coefficient
#' @param A gene sets
#' @param B gene sets

jaccard_similarity <- function(A, B) {
  intersection = length(intersect(A, B))
  union = length(A) + length(B) - intersection
  return (intersection/union)
}


