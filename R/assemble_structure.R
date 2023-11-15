
assemble_structure <- function(flamingo_high_res_obj,
                               flamingo_backbone_prediction,
                               flamingo_intra_domain_prediction,
                               alpha,
                               inf_dist,
                               max_iter=500)
                               
{
  #### read all the infomation
  n = nrow(flamingo_backbone_prediction@coordinates) # numbers of domain
  backbone_id <- flamingo_backbone_prediction@id # domain id
  backbone <- as.matrix(flamingo_backbone_prediction@coordinates) # backbone coord
  n_point_per_domain = flamingo_intra_domain_prediction[[1]]@input_n # domain size

  #### mean distance of sub-diagonal
  scaler <- as.matrix(dist(backbone))
  scaler <- mean(scaler[which(row(scaler)+1 == col(scaler))])

  #### read intra-domain structure
  all_points <- list() # all structures in each domain, with invalid point
  val_id_list <- list()
  val_domain <- c()
  for(counter in backbone_id)
  {
    if(!counter %in% names(flamingo_intra_domain_prediction)) next

    tmp_domain_res = flamingo_intra_domain_prediction[[as.character(counter)]]
    tmp_id <- tmp_domain_res@id # valid id
    p_t <- as.matrix(tmp_domain_res@coordinates) # all points
    p_t_valid = p_t[tmp_id,]
    tmp_center <- backbone[which(backbone_id==counter),]

    # rescale and put to the new center
    old_center <- apply(p_t_valid,2,mean) # based on valid point
    radius <- max(apply(p_t_valid,1,function(x){norm(x-old_center,'2')}))
    new_loc <- t(apply(p_t,1,function(x){
      
      (x-old_center)/radius*scaler+tmp_center
      
    }))

    # save structure
    all_points[[as.character(counter)]] <- as.data.frame(new_loc)
    val_id_list[[as.character(counter)]] <- tmp_id
    val_domain = c(val_domain, counter)

  }

  total_struc <- backbone[match(val_domain,backbone_id),]

  # get the id list
  id_list = list()
  for (i in 1:length(val_domain)){
    tmp_id <- (val_domain[i]-1)*n_point_per_domain + val_id_list[[i]]
    id_list[[as.character(val_domain[i])]] = tmp_id

  }

  # start_id_list
  start_id_list = sapply(val_id_list,function(x)x[1])
  end_id_list = sapply(val_id_list,function(x) tail(x,n=1))

  #### prepare input data
  # domain 2 to domain 1 direction, center
  direction <- total_struc[2,]-total_struc[1,]
  direction <- as.vector(direction/norm(direction,'2'))

  # domain 1 valid tail to domain 1 center direction
  domain_1_tail_coord = all_points[[1]][end_id_list[1],]
  tail_direction <- domain_1_tail_coord - total_struc[1,]
  tail_direction <- as.vector(as.matrix(tail_direction/norm(tail_direction,'2')))

  # rotation matrix to place the first domain
  r_mat <- rotation_matrix(matrix(tail_direction,ncol=1),as.matrix(direction,nrow=1))
  all_points[[1]] <- rotate(all_points[[1]],total_struc[1,],r_mat)

  # whole pd
  input_if = flamingo_high_res_obj@IF
  input_if <- input_if + t(input_if)
  diag(input_if) <- diag(input_if)/2
  pd <- if2pd(input_if, alpha, inf_dist)
  error <- 10
  d = get_dist_vec(all_points,id_list,pd,inf_dist)
  t = 1e-4

  # rotation
  g_p <- matrix(0,length(all_points)-1,3)
  y_s_p <-  matrix(0,length(all_points)-1,3)
  y_e_p <-  matrix(0,length(all_points)-1,3)
  error_p = error_change <- 1
  N <- dim(do.call(rbind,all_points))[1]
  p_p <- matrix(0,N,3)
  iter <- 1
  while(error>1e-4&error_change>1e-5){
    y_s <- get_point(all_points,start_id_list)[2:length(val_domain),]
    y_e <- get_point(all_points,end_id_list)[1:length(val_domain)-1,]
    d_tild <- evaluate_dist(y_s,y_e)
    g <- 4*(d_tild-d)*(y_s-y_e)
    
    y = g-g_p
    s = y_s-y_s_p
    t_k = (sum(diag(t(s)%*%s)))/(sum(diag(t(s)%*%y)))
    
    y_s_prim <- y_s-t_k*g
    g_p <- g
    y_s_p <- y_s

    #update structure
    for(i in 2:length(all_points)){
      tmp_points <- all_points[[i]]
      tmp_id = start_id_list[i]
      tmp_start <- tmp_points[tmp_id,]
      tmp_center <- total_struc[i,]
      old_direction <- as.vector(as.matrix(tmp_start-tmp_center))+1e-3
      new_direction <- y_s_prim[i-1,] - tmp_center
      r_mat <- rotation_matrix(old_direction,new_direction)
      all_points[[i]] <- rotate(all_points[[i]],total_struc[i,],r_mat)
    }
    #update structure
    if(iter == 1){
      y_e_prim <- y_e-t_k*g
      y_e_p <- y_e
      for(i in 1:(length(all_points)-1)){
        tmp_points <- all_points[[i]]
        tmp_id = end_id_list[i]
        tmp_end <- tmp_points[tmp_id,]
        tmp_center <- total_struc[i,]
        old_direction <- as.vector(as.matrix(tmp_end-tmp_center))+1e-3
        new_direction <- y_e_prim[i,] - tmp_center
        r_mat <- rotation_matrix(old_direction,new_direction)
        all_points[[i]] <- rotate(all_points[[i]],total_struc[i,],r_mat)
      }
    }
    iter <- iter + 1
    cur_p <- do.call(rbind,all_points)
    all_points_err <- norm(cur_p-p_p,'F')/norm(cur_p,'F')
    p_p <- cur_p
    error <- norm(d_tild-d,'2')/norm(d,'2')
    error_change <- abs(error-error_p)
    error_p <- error
    # print(error)
    if(all_points_err<0.0001){
      break
    }
    if(iter > max_iter){
      break
    }
    
  }

  all_points_valid = data.frame()
  for(i in 1:length(all_points)){
    all_points_valid = rbind(all_points_valid, all_points[[i]][val_id_list[[i]],])
  }

  # smooth
  # all_points <- smt(all_points)
  res = data.frame(frag_id = unlist(id_list),x=all_points_valid[,1],y=all_points_valid[,2],z=all_points_valid[,3])
  return(res)
}
