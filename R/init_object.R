#' Create an object named flamingo for data organization and an object named flamingo_prediction for result organization
#' @keywords flamingo_object
#' @return Objects named flamingo and flamingo_prediction are ready for use.
#' @export
# library(Matrix)
setClass("flamingo", slots=list(IF="sparseMatrix", n_frag='numeric',chr_name='character'))
setClass("flamingo_prediction", slots=list(id="numeric", coordinates='matrix',input_n='numeric'))