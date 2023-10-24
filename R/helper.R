

#' Creates co-occurrence matrix given training trials.
#'
#' Given a training order (list of words and objects per trial), returns a matrix of
#' tallied word-object co-occurrences across the trials.
#'
#' @return A list of trials (with nested words and objects per trial) and perf (P(correct referent | word))
#' @export
create_cooc_matrix <- function(train) {
  Nwords = length(unique(unlist(train$words)))
  Nobjs = length(unique(unlist(train$objs)))
  M = matrix(0, nrow=Nwords, ncol=Nobjs)
  rownames(M) = sort(unique(unlist(train$words)))
  colnames(M) = sort(unique(unlist(train$objs)))
  # iterate over training trials
  for(i in 1:length(train$words)) {
    M[train$words[[i]], train$objs[[i]]] = M[train$words[[i]], train$objs[[i]]] + 1
  }
  return(M)
}

# DELETE, OR NEEDED?
# given a model knowledge matrix, retrieve Luce choice (proportion correct) per item
get_perf <- function(m) {
  perf <- rep(0, nrow(m))
  names(perf) <- rownames(m)
  for (ref in colnames(m)) {
    if (!(ref %in% rownames(m))) {
      next
    }
    correct <- m[ref, ref]
    total <- sum(m[ref,])
    if (total == 0) {
      next
    }
    perf[ref] <- correct / total
  }
  return(perf)
}
