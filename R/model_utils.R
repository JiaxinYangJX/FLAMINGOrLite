
############# utils ##################

#### contact frequency to pairwise distance
if2pd <- function(input_if,
                  alpha = -0.25,
                  inf_dist = 2)
{

  pd = input_if ^ (alpha)
  pd[which(pd==Inf|is.na(pd))] = inf_dist
  diag(pd) = 0

  return(pd)

}


####
convert_index <- function(x)
{

  df <- expand.grid(as.vector(x),as.vector(x))

  df[,3] <- (df[,2] == df[,1])*2-1

  return(df)

}


#### 3D coordinate to gram matrix
coord2gram <- function(P)
{

  return(as.matrix(P%*%t(P)))

}


#### pairwise distance to gram matrix
pd2gram <- function(pd)
{
  n = nrow(pd)
  H = Matrix::Diagonal(x=rep(1,n))-1/n*rep(1,n)%*%t(rep(1,n))
  M = -1/2 * H %*% pd^2 %*% H
  M = as.matrix(M)

  return(M)

}


#### get all mearsurement set from contact frequency
get_measurement_set <- function(input_if,
                                if_threshold = 0)
{

  # sparse matrix
  input_if = as(input_if,'sparseMatrix')
  col_j = findInterval(seq(input_if@x)-1,input_if@p[-1])+1
  row_i = input_if@i+1
  x_x = input_if@x
  df = data.frame(row_i,col_j,x_x)

  # contact frequency > 0 (if_threshold) ; off-diagonal
  omega = subset(df[,1:2],df[,3] > if_threshold & df[,1] > df[,2])
  omega = as.matrix(omega)

  return(omega)

}


#### get the a bind of omega set
get_bind_set <- function(omega,
                         bind_size)
{

  bind_term = which(omega[,1]-omega[,2] == bind_size)
  n_omega_bind = length(bind_term)

  omega_bind = omega[bind_term,,drop=FALSE] # avoid vector with only one sample

  return(omega_bind)

}


#### get the downsampled set of omega
get_sample_set <- function(omega,
                           sample_rate)
{
  n_omega = dim(omega)[1]

  omega_sample <<- omega[sample(1:n_omega, sample_rate*n_omega),]

  return(omega_sample)

}


#### get pre-calculated elements for adjoint linear projection
get_element_adjoint_linear <- function(omega)
{

  n_omega = nrow(omega)

  # prepare for A*
  w_list = lapply(1:n_omega,function(x){
    cbind(convert_index(omega[x,]),x)
  })
  w_list = do.call(rbind,w_list)

  all_element = unique(data.frame(w_list[,1:2]))
  loc = prodlim::row.match(data.frame(w_list[,1:2]),all_element)

  func_list <- split(data.frame(w_list[,c(3,4)]),loc)

  res = list("func_list" = func_list, "all_element" = all_element)
  return(res)

}





############# forward ################

linear_proj <- function(omega,x)
{
  # Linear projection: A(X)
  # fi(X) = <X, omega_i>
  # where i = (i_1, i_2) represents the index of DNA fragment pairs
  # omega_i = e_i1i1 + e_i2i2 - e_i1i2 - e_i2i1
  # where e_ab represents a matrix which has 1 at entry (a,b) and 0 otherwise.

  proj <- apply(omega,1,function(y){

    y <- as.numeric(y)

    p <- x[y[1],y[1]]+x[y[2],y[2]]-x[y[1],y[2]]-x[y[2],y[1]]

    return(p)
  })

  return(proj)

}


linear_proj_adj <- function(x,func_list,all_element,n)
{
  # Linear projection: A(X)
  # fi(X) = <X, omega_i>
  # f*(X) = SUM(X_i * omega_i)

  tmp <- sapply(func_list,function(y){

    sum(x[y[,2]] * y[,1])

  })

  mat = Matrix::sparseMatrix(i=all_element[,1],j=all_element[,2],x=tmp,dims=c(n,n))

  return(mat)

}



################# Gradient ##########################

#### trace(PP^T)
# 2P
grad_rank <- function(P)
{
  return(2*P)
}



#### || A(PP^T) + y ||^2
# A^* (A(PP^T) + y) P
grad_linear <- function(P,y,omega,func_list,all_element)
{

  n = nrow(P)
  x = coord2gram(P) # PP^T

  x = linear_proj(omega,x) # A(PP^T)

  x = linear_proj_adj(x + y, func_list,all_element,n)

  x = x %*% P

  return(x)

}



#### || A(PP^T) - b + gamma ||^2
# A^* (A(PP^T) - b + gamma) P
grad_A <- function(P,gamma,b,omega_sample,func_list_sample,all_element_sample){

  return(grad_linear(P,-b+gamma,omega_sample,func_list_sample,all_element_sample))

}


#### || B(PP^T) - d ||^2
# B^* (B(PP^T) - d) P
grad_B <- function(P,d,omega_subdiag,func_list_subdiag,all_element_subdiag)
{

  return(grad_linear(P,-d,omega_subdiag,func_list_subdiag,all_element_subdiag))

}


#### BB decent
BB_decent <- function(P_prev,
                      P_curr,
                      grad_prev,
                      grad_curr)
{

  s = as.matrix(P_curr - P_prev)
  y = as.matrix(grad_curr - grad_prev)

  # step size
  t_k = sum(diag( t(s) %*% s )) / sum(diag( t(s) %*% y ))

  # update
  P_next = as.matrix(P_curr - t_k * grad_curr)

  return(P_next)

}


################ FLAMINGO ######################

#### gradient of FLAMINGO
flamingo_grad <- function(P,
                          omega_sample,
                          omega_subdiag,
                          func_list_sample,
                          func_list_subdiag,
                          all_element_sample,
                          all_element_subdiag,
                          gamma,b,d,lambda,r)
{

  g <- grad_rank(P) +
       lambda * grad_B(P,d,omega_subdiag,func_list_subdiag,all_element_subdiag) +
       r * grad_A(P,gamma,b,omega_sample,func_list_sample,all_element_sample)

  return(g)
}


#### FLAMINGO optimization worker
flamingo_worker <- function(omega_sample,
                            omega_subdiag,
                            func_list_sample,
                            func_list_subdiag,
                            all_element_sample,
                            all_element_subdiag,
                            b,
                            d,
                            n,
                            lambda,
                            r,
                            error_threshold=1e-3,
                            max_iter=300)

{
  # initialization
  q = 3
  P_curr <- matrix(rnorm(n*q),n,q)
  P_prev <- matrix(rnorm(n*q),n,q)
  gamma <- rep(0,nrow(omega_sample))
  error <- 10

  for (iter in 1:max_iter)
  {
    # calculate gradient
    if(iter == 1){
      grad_prev <- flamingo_grad(P_prev,
                                 omega_sample,
                                 omega_subdiag,
                                 func_list_sample,
                                 func_list_subdiag,
                                 all_element_sample,
                                 all_element_subdiag,
                                 gamma,b,d,lambda,r)
    }else{
      grad_prev <- grad_curr
    }

    grad_curr <- flamingo_grad(P_curr,
                               omega_sample,
                               omega_subdiag,
                               func_list_sample,
                               func_list_subdiag,
                               all_element_sample,
                               all_element_subdiag,
                               gamma,b,d,lambda,r)

    # BB decent
    P_next = BB_decent(P_prev, P_curr, grad_prev, grad_curr)

    # update
    P_prev <- P_curr
    P_curr <- P_next
    error <- norm(P_curr - P_prev,"f")

    # gc()

    if (error < error_threshold) break

  }

  return(P_next)

}



################ rotation ######################

#### generate rotation matrix
rotation_matrix = function(x,y)
{

  u=x/sqrt(sum(x^2))

  v=y-sum(u*y)*u
  v=v/sqrt(sum(v^2))

  cost=sum(x*y)/sqrt(sum(x^2))/sqrt(sum(y^2))

  sint=sqrt(max(1-cost^2,0));
  #print(cost)
  diag(length(x)) - u %*% t(u) - v %*% t(v) +
    cbind(u,v) %*% matrix(c(cost,-sint,sint,cost), 2) %*% t(cbind(u,v))
}


#### apply rotation
rotate <- function(x,c,r)
{

  tmp_c <- t(apply(x,1,function(y){y-c}))
  z = as.matrix(tmp_c) %*% r
  z = t(apply(z,1,function(y){y+c}))
  return(z)

}


#### get the start points of all domain
# get_start_point <- function(all_points,val_id_list){
#   start_id = sapply(val_id_list,function(x) x[1])

#   # get start point of each domain, valid
#   start_coord = t(mapply(function(x,y) x[y,],all_points,start_id))

#   start_coord = matrix(unlist(start_coord),ncol=3)
#   start_coord = start_coord[-1,]

#   return(start_coord)
# }


# #### get the end points of all domain
# get_end_point <- function(all_points,val_id_list){
#   end_id = sapply(val_id_list,function(x) tail(x,n=1))

#   # get start point of each domain, valid
#   end_coord = t(mapply(function(x,y) x[y,],all_points,end_id))

#   end_coord = matrix(unlist(end_coord),ncol=3)
#   end_coord = end_coord[1:nrow(end_coord)-1,]

#   return(end_coord)
# }


get_point <- function(all_points,id_list){

  # get start point of each domain, valid
  coord = t(mapply(function(x,y) x[y,],all_points,id_list))

  coord = matrix(unlist(coord),ncol=3)

  return(coord)
}


ave_dist <- function(r,pd){
  res = pd[row(pd)-col(pd)==r | col(pd)-row(pd)==r]
  return(mean(res[res<Inf],na.rm=T))
}


get_dist_vec <- function(all_points,id_list,pd,inf_dist){
  n <- length(all_points)
  start_id <- sapply(id_list[2:(n)],function(y){y[1]})
  end_id <- sapply(id_list[1:(n-1)],function(y){tail(y,n=1)})
  dist_vec <- c()

  # observed distance between the adjoint point in two domain
  for(i in 1:(n-1)){

    tmp_dist = pd[start_id[i],end_id[i]]

    if(is.na(tmp_dist) | tmp_dist == inf_dist){
      # observation not available
      tmp_dist = ave_dist(start_id[i]-end_id[i],pd)

    }else if(tmp_dist == 0){
      tmp_dist = 0.001
    }

    dist_vec = c(dist_vec,tmp_dist)
  }

  return(dist_vec)
}


#start optimization

evaluate_dist <- function(start_point,end_point){

  apply(start_point-end_point,1,function(x){norm(x,'2')})

}

smt <- function(o){
  o_smt <- o
  for(i in 2:c(dim(o_smt)[1]-2)){

    o_smt[i,] <- apply(o_smt[(i-1):(i+1),],2,mean)


  }
  return(o_smt)
}
