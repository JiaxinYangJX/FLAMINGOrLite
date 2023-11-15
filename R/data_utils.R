




construct_obj_from_hic <- function(hic_file,
                                   resolution,
                                   chr_name,
                                   normalization)
{
  #### get data
  options(scipen = 999)
  chr_number <- gsub("chr","",chr_name)
  normalized_data = strawr::straw(normalization,hic_file,chr_number,chr_number,unit='BP',binsize=resolution)

  #### reformat contact frequency matrix
  n <- max(normalized_data[,2])/resolution + 1
  i_ind <- (normalized_data[,1]/resolution) + 1
  j_ind <- (normalized_data[,2]/resolution) + 1
  input_if = Matrix::sparseMatrix(i=i_ind,j=j_ind,x=normalized_data[,3],dims=c(n,n))

  #### generate object
  res = new('flamingo',IF=input_if,n_frag=n,chr_name=chr_name)
  return(res)

}



# construct_obj_from_mcool <- function(hic_file,
#                                      resolution,
#                                      chr_name,
#                                      normalization)
# {
#   options(scipen = 999)
#   all_dir = rhdf5::h5ls(mcool_file)
#   parent_dir = all_dir[1,2]
#   target_dir = paste(c("",parent_dir,resolution),collapse='/')
#   mcool_dat = rhdf5::h5read(mcool_file,target_dir)
#   available_normalization = setdiff(names(mcool_dat$bins),c('chrom','start','end','weight'))
#   if(!normalization %in% available_normalization){
#     stop(
#       paste('Normalization method not exist in the .mcool data! Must be one of: ',paste(available_normalization,collapse=', '))
#          )
#   }
#   csr_rawcount = data.frame(bin_1 = mcool_dat$pixels$bin1_id,
#                             bin_2 = mcool_dat$pixels$bin2_id,
#                             value = mcool_dat$pixels$count)
#   normalization_file = mcool_dat$bins[[normalization]]
#   chr_id <- which(mcool_dat$bins$chrom == chr_name)
#   offset <- mcool_dat$indexes$chrom_offset[which(mcool_dat$chroms$name==chr_name)]
#   n <- length(chr_id)
#   csr_rawcount <- subset(csr_rawcount,csr_rawcount[,1] %in% chr_id & csr_rawcount[,2] %in% chr_id)
#   csr_rawcount[,1] <- csr_rawcount[,1]-offset+1
#   csr_rawcount[,2] <- csr_rawcount[,2]-offset+1
#   csr_rawcount <- as.matrix(csr_rawcount)
#   normalization_file <- normalization_file[chr_id]
#   for(i in 1:(dim(csr_rawcount)[1])){
#     csr_rawcount[i,3] <- csr_rawcount[i,3]/(normalization_file[csr_rawcount[i,1]]*normalization_file[csr_rawcount[i,2]])
#   }
#   input_if <- Matrix::sparseMatrix(i=csr_rawcount[,1],j=csr_rawcount[,2],x=csr_rawcount[,3],dims=c(n,n))
#   if(n<n_row){
#     input_if <- as.matrix(input_if)
#   }else{
#     input_if <- convert_huge_mat(input_if)
#   }
#   input_if <- input_if + t(input_if)
#   diag(input_if) <- diag(input_if)/2
#   pd <- input_if^(alpha)
#   res = new('flamingo',IF=input_if,PD=pd,n_frag=n,chr_name=chr_name)
#   return(res)

# }


check_data_availability <- function(input_if){
  input_if <- as.matrix(input_if)
  invalid_id <- which(apply(input_if,1,max)==0)
  if(length(invalid_id) == dim(input_if)[1]){
    return(F)
  }else{
    return(T)
  }
}




flamingo_obj2if <- function(flamingo_obj)
{

  input_if = as.matrix(flamingo_obj@IF)
  input_if <- input_if + t(input_if)
  diag(input_if) <- diag(input_if)/2
  return(input_if)

}


flamingo_pred_obj2coord <- function(flamingo_prediction,
                                    genomic_loc)
{

  genomic_loc = genomic_loc[flamingo_prediction@id,]
  coord = flamingo_prediction@coordinates[flamingo_prediction@id,]

  res = cbind(genomic_loc,coord)
  colnames(res) = c('chr','start','end','x','y','z')

  return(res)

}



divide_domain <- function(flamingo_high_res_obj,
                          flamingo_low_res_obj,
                          domain_res,
                          frag_res,
                          temp_folder)
{

  #### create temp folder
  options(scipen = 999)

  res = list()

  bin_size = domain_res/frag_res
  n = flamingo_high_res_obj@n_frag
  chr = flamingo_high_res_obj@chr_name
  input_if = flamingo_high_res_obj@IF

  # write out high resolution domain matrix
  print('Writing out temp files...')
  n_domain = ceiling(n/bin_size)
  for(i in 1:n_domain){

    # contact frequency
    start_id <- (i-1)*bin_size+1
    end_id <- min(n,i*bin_size)
    tmp_input_if <- as.matrix(input_if[start_id:end_id,start_id:end_id])
    tmp_input_if <- tmp_input_if + t(tmp_input_if)
    diag(tmp_input_if) <- diag(tmp_input_if)/2

    # 1D coordinate
    start_loc = ((start_id:end_id)-1)*frag_res+1
    end_loc = (start_id:end_id)*frag_res
    tmp_frag = data.frame(chr=chr,start=start_loc,end=end_loc)

    # availabe data
    if (check_data_availability(tmp_input_if)){
      write.table(tmp_input_if,paste0(temp_folder,"/domain_data/IF_domain_",i,".txt"),col.names = F,row.names = F,sep="\t",quote=F)
      write.table(tmp_frag,paste0(temp_folder,"/genomic_loc/genomic_loc_domain_",i,".txt"),col.names = F,row.names = F,sep="\t",quote=F)
    }
    # write.table(tmp_input_if,paste0("temp/domain_data/IF_domain_",i,".txt"),col.names = F,row.names = F,sep="\t",quote=F)
    # write.table(tmp_frag,paste0("temp/genomic_loc/genomic_loc_domain_",i,".txt"),col.names = F,row.names = F,sep="\t",quote=F)

  }

  # write out low resolution backbone matrix
  input_if = as.matrix(flamingo_low_res_obj@IF)
  input_if <- input_if + t(input_if)
  diag(input_if) <- diag(input_if)/2
  write.table(input_if,paste0(temp_folder,"/IF_backbone.txt"),col.names = F,row.names = F,sep="\t",quote=F)

}

