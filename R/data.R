#' Import contact matrix
#'
#' @param filename the filename containing the matrix, csv with first column as row.names
#' @param replace_NA_w_zero replace NAs with zero
#'
#' @return a matrix
#' @export
import_matrix = function(filename, replace_NA_w_zero=FALSE, ...) {
  tmp = read.csv(filename, stringsAsFactors=FALSE, ...)
  n = ncol(tmp)-1
  M = data.matrix(tmp[,1+(1:n)])
  rownames(M) = tmp[,1]
  if (isTRUE(replace_NA_w_zero)) {
    M[is.na(M)] = 0
  }
  M
}
