#' backbone reconstruction of flamingorlite
#' Core function for 3D genome structure reconstruction using low-rank matrix completion
#' @param temp_folder tempfolder.
#' @param sample_rate Fraction of available entries in Hi-C to be used during the reconstruction. Default = 0.75.
#' @param lambda Weights for all sampled entries. Default = 10.
#' @param r Weights for distance between consecutive points. Default = 1.
#' @param max_dist Maximum allowed distance betwee two consecutive points. Default = 0.1
#' @param error_threshold Error thresholds for reconstruction. Default = 1e-4.
#' @param max_iter Maximum iterations. Default = 500.
#' @param alpha Convertion factor between interaction frequency and pairwise distance. Default = -0.25.
#' @param inf_dist Maximun allowed distance betwee any two points. Default = 3.
#' @keywords flamingo_backbone
#' @return A flamingo_prediction object containing the fragment id and 3D coordinates
#' @export
flamingo_backbone <- function(temp_folder,
                              sample_rate,
                              lambda,
                              r,
                              max_dist,
                              error_threshold,
                              max_iter,
                              alpha,
                              inf_dist)
{
  input_if = read.table(paste0(temp_folder,"/IF_backbone.txt"))
  flamingo_backbone_prediction = flamingo_basic(input_if = input_if,
                                                sample_rate = sample_rate,
                                                lambda = lambda,
                                                r = r,
                                                max_dist = max_dist,
                                                error_threshold = error_threshold,
                                                max_iter = max_iter,
                                                alpha = alpha,
                                                inf_dist = inf_dist)

  return(flamingo_backbone_prediction)

}