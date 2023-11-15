
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