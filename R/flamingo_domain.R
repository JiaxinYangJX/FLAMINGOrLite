


flamingo_domain <- function(temp_folder,
                            sample_rate,
                            lambda,
                            r,
                            max_dist,
                            error_threshold,
                            max_iter,
                            alpha,
                            inf_dist,
                            nThread)
{

  # get data domain id
  domain_id = sapply(strsplit(list.files(paste0(temp_folder,'/domain_data/')),
                      split='[_.]+'), function(x) x[3])

  # parallel setting
  require(parallel)
  cl <- parallel::makeCluster(nThread)
  parallel::clusterEvalQ(cl, {
    source('model_utils.R')
    source('flamingo_basic.R')
    source('init_object.R')
    # library(FLAMINGOrLite)
  })
  parallel::clusterExport(cl,c('temp_folder','sample_rate','lambda','r','max_dist',
                               'error_threshold','max_iter','alpha','inf_dist'),
                               envir=environment())

  # parallel computing
  res <- parallel::parSapply(cl,domain_id,function(index){
    # read data
    tmp_input_if <- read.table(paste0(temp_folder,"/domain_data/IF_domain_",index,".txt"))
    tmp_res =  flamingo_basic(tmp_input_if,sample_rate,lambda,r,max_dist,error_threshold,max_iter,alpha,inf_dist)
    return(tmp_res)

  })
  parallel::stopCluster(cl)

  # results
  names(res) <- domain_id

  return(res)
}





