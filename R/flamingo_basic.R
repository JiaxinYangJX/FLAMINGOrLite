#' main function of flamingo
#' Core function for 3D genome structure reconstruction using low-rank matrix completion
#' @param input_if Input interaction frequency matrix.
#' @param sample_rate Fraction of available entries in Hi-C to be used during the reconstruction. Default = 0.75.
#' @param lambda Weights for all sampled entries. Default = 10.
#' @param r Weights for distance between consecutive points. Default = 1.
#' @param max_dist Maximum allowed distance betwee two consecutive points. Default = 0.1
#' @param error_threshold Error thresholds for reconstruction. Default = 1e-4.
#' @param max_iter Maximum iterations. Default = 500.
#' @param alpha Convertion factor between interaction frequency and pairwise distance. Default = -0.25.
#' @param inf_dist Maximun allowed distance betwee any two points. Default = 2.
#' @keywords flamingo_basic
#' @return A flamingo_prediction object containing the fragment id and 3D coordinates
#' @export


flamingo_basic <- function(input_if,
                           sample_rate = 0.75,
                           lambda = 10,
                           r = 1,
                           max_dist,
                           error_threshold,
                           max_iter,
                           alpha = -0.25,
                           inf_dist = 4)
{
  require(Matrix)

  #### generate pairwise distance
  input_if = as.matrix(input_if)
  input_if[which(is.na(input_if))] = 0

  pd = if2pd(input_if,alpha,inf_dist)
  n = nrow(pd)

  #### invalid idx
  rm_id = which(apply(input_if,1,max)==0)


  #### identify the gram matrix of the input data
  M = pd2gram(pd)


  #### define measurement set omega
  omega = get_measurement_set(input_if)
  n_omega = dim(omega)[1]


  #### sub diagonal set
  diag_term <- which(omega[,2]-omega[,1]==1)
  omega_diag = omega[diag_term,]
  omega <<- omega[unique(c(diag_term,sample(1:n_omega,sample_rate*n_omega))),]
  n_omega <- dim(omega)[1]
  

  if(length(diag_term)==1){
    n_omega_diag=1
    omega_diag = matrix(omega_diag,ncol=2)
  }
  else if(length(diag_term)==0){
    return(NULL)
  }else{
    n_omega_diag <- dim(omega_diag)[1]
  }


  #### pre-calculate related data
  # prepare for A*
  precal_sample = get_element_adjoint_linear(omega)
  func_list_sample = precal_sample$func_list
  all_element_sample = precal_sample$all_element

  # prepare for B*
  precal_subdiag = get_element_adjoint_linear(omega_diag)
  func_list_subdiag = precal_subdiag$func_list
  all_element_subdiag = precal_subdiag$all_element


  #### pre-calculate the b and d
  b = linear_proj(omega,M)

  d = linear_proj(omega_diag,M)

  # control the sub-diagonal
  for(i in 1:length(d)){
    d[i] = min(d[i],max_dist)
  }


  #### run flamingo
  P <- flamingo_worker(omega,
                       omega_diag,
                       func_list_sample,
                       func_list_subdiag,
                       all_element_sample,
                       all_element_subdiag,
                       b,d,n,lambda,r,error_threshold,max_iter)


  #### keep the valid samples
  if(length(rm_id)>0){
    frag_id <- (1:n)[-rm_id]
  }else{
    frag_id <- 1:n
  }

  return(new('flamingo_prediction',id = frag_id, coordinates = P,input_n = n))

}
