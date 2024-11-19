#' Main function of hierarchical-FLAMINGO
#' @param hic_data Input Hi-C data. Now only support .hic
#' @param file_format Foramt of the input hic data. Now only support .hic.
#' @param domain_res Size of the domains in bps, e.g. 1e6. Try strawr::readHicBpResolutions() to see available resolutions.
#' @param frag_res Size of the fragment in bps, e.g. 5e3. Try strawr::readHicBpResolutions() to see available resolutions.
#' @param chr_name Name of the chromosome, e.g. chr1. Try strawr::readHicChroms() to see available chromosomes.
#' @param normalization Normalization method in .hic file. Try strawr::readHicNormTypes() to see available methods. Could be 'NONE'.
#' @param nThread Number of thread avalable for the reconstruction. Default = 1.
#' @param sample_rate Fraction of available entries in Hi-C to be used during the reconstruction. Default = 0.75.
#' @param lambda Weights for all sampled entries. Default = 10.
#' @param r Weights for distance between consecutive points. Default = 1.
#' @param max_dist Maximum allowed distance betwee two consecutive points. Default = 0.01
#' @param alpha Convertion factor between interaction frequency and pairwise distance. Default = -0.25.
#' @param inf_dist Maximun allowed distance betwee any two points. Default = 2.
#' @param error_threshold Error thresholds for reconstruction. Default = 1e-3.
#' @param max_iter Maximum iterations. Default = 500.
#' @keywords flamingo_main
#' @return A data.frame containing the FLAMINGO predicted 3D structure.
#' @export

flamingo_main <- function(hic_data,
                          file_format,
                          domain_res,
                          frag_res,
                          chr_name,
                          normalization,
                          nThread = 20,
                          sample_rate = 0.75,
                          lambda = 10,
                          r = 1,
                          max_dist = 0.01,
                          alpha = -0.25,
                          inf_dist = 2,
                          error_threshold = 1e-3,
                          max_iter=500)
{
  #### check Args
  b = Sys.time()
  bin_size = domain_res/frag_res
  if(bin_size != round(bin_size)){
    stop('The domain resolution is not an integer multiple of the fragment resolution!')
  }

  library(parallel)
  library(Matrix)
  # source('flamingo_domain.R')
  # source('flamingo_backbone.R')
  # source('data_utils.R')
  # source('init_object.R')
  # source('model_utils.R')
  # source('flamingo_basic.R')
  # source('assemble_structure.R')

  #### create temp folder
  temp_folder = "./tempfolder"
  dir.create(temp_folder)
  dir.create(paste0(temp_folder,'/domain_data'))
  # dir.create(paste0(temp_folder,'/genomic_loc'))

  #### generate flamingo object
  print('Preparing data...')
  a = Sys.time()
  if(file_format=='hic'){
    flamingo_low_res_obj <- construct_obj_from_hic(hic_file = hic_data,
                                                   resolution = domain_res,
                                                   chr_name = chr_name,
                                                   normalization = normalization)
    
    flamingo_high_res_obj <- construct_obj_from_hic(hic_file = hic_data,
                                                    resolution = frag_res,
                                                    chr_name = chr_name,
                                                    normalization = normalization)
  }else if(file_format=='mcool'){
    flamingo_low_res_obj <- construct_obj_from_mcool(mcool_file = hic_data,
                                                     resolution = domain_res,
                                                     chr_name = chr_name,
                                                     normalization = normalization)
    
    flamingo_high_res_obj <- construct_obj_from_mcool(mcool_file = hic_data,
                                                      resolution = frag_res,
                                                      chr_name = chr_name,
                                                      normalization = normalization)
  }else{
    stop("file format must be .hic or .mcool")
  }
  print(paste('Finished time: ',round(as.numeric(difftime(Sys.time(),a,units='mins')),digits=2), ' mins'))
  
  
  #### Divide domain dataset
  print('Dividing domains...')
  a = Sys.time()
  divide_domain(flamingo_high_res_obj = flamingo_high_res_obj,
                flamingo_low_res_obj = flamingo_low_res_obj,
                domain_res = domain_res,
                frag_res = frag_res,
                temp_folder = temp_folder)
  print(paste('Finished time: ',round(as.numeric(difftime(Sys.time(),a,units='mins')),digits=2), ' mins'))


  #### Reconstruct backbone
  print('Reconstructing backbones...')
  a = Sys.time()
  flamingo_backbone_prediction = flamingo_backbone(temp_folder,sample_rate,lambda,r,max_dist,error_threshold,max_iter,alpha,inf_dist)
  print(paste('Finished time: ',round(as.numeric(difftime(Sys.time(),a,units='mins')),digits=2), ' mins'))


  #### Reconstruct domain in parallel
  print('Reconstructing intra-domain structures...')
  a = Sys.time()
  flamingo_intra_domain_prediction = flamingo_domain(temp_folder,sample_rate,lambda,r,max_dist,error_threshold,max_iter,alpha,inf_dist,nThread)
  print(paste('Finished time: ',round(as.numeric(difftime(Sys.time(),a,units='mins')),digits=2), ' mins'))


  print('Assembling structures...')
  a = Sys.time()
  res = assemble_structure(flamingo_high_res_obj,flamingo_backbone_prediction,flamingo_intra_domain_prediction,alpha,inf_dist,max_iter)
  print(paste('Finished time: ',round(as.numeric(difftime(Sys.time(),a,units='mins')),digits=2), ' mins'))

  #### Reformat results

  res$chr = chr_name
  res$start = (res$frag_id-1) * frag_res
  res$end = res$frag_id * frag_res
  res = res[,c('chr','start','end','x','y','z')]
  print(paste('Reconstruction successful! Finished time: ',round(as.numeric(difftime(Sys.time(),b,units='mins')),digits=2), ' mins'))

  return(res)
}
