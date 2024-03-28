library(data.table)
library(dplyr)
library(stringr)
library(shiny)
library(igraph)
library(ggplot2)


# dp_list is a list of all your data perturbation result
subset_Graph_d1 <- function(dp_list,keep_trait,B=2000,check=TRUE,show=TRUE,maxit=10000){
  set.seed(721)
  keep_trait_id = is.element(dp_list$trait_vec,keep_trait)
  
  warning=warn_len=NULL
  B_list = sample(1:length(dp_list$obs_graph_list),B,replace = FALSE)
  total_B = length(dp_list$obs_graph_list)
  new_obs_graph_list = lapply(dp_list$obs_graph_list[B_list],function(G) {
    G = G[keep_trait_id,keep_trait_id];
    #G[,BW_id] = 0;
    GG = G * t(G); diag(GG) = 0;
    diag(G)=rowSums(GG);
    G})
  
  new_dir_graph_list1 = vector(mode='list', length=length(new_obs_graph_list))
  for(i in 1:length(new_obs_graph_list)){
    G_obs = new_obs_graph_list[[i]]
    iter = 0
    G_dir = G_obs %*% solve(diag(nrow(G_obs))+G_obs);
    G_dir0 = G_dir;
    converge = 0;
    while(any(abs(diag(G_dir)) > 1e-4) & iter < maxit){
      G_dir = G_obs %*% solve(diag(nrow(G_obs))+G_obs);
      GG = G_dir * t(G_obs); diag(GG) = 0;
      diag(G_obs) = rowSums(GG);
      iter = iter + 1;
    }
    if(iter>=maxit){
      G_dir = G_dir0
      converge = 1
    }
    new_dir_graph_list1[[i]]$converge = converge
    new_dir_graph_list1[[i]]$G_dir = G_dir
  }
  
  
  new_dir_graph_list = lapply(new_dir_graph_list1, function(x){x$G_dir})
  new_dir_graph_conv_vec = unlist(lapply(new_dir_graph_list1, function(x){x$converge}))
  conv_len = sum(new_dir_graph_conv_vec)
  new_obs_pval_graph_list = dp_list$obs_graph_pval_list[B_list]
  max_eigen_list = lapply(new_dir_graph_list, function(x) {max(abs(eigen(x)$values))})
  rm_dp_id = which(unlist(max_eigen_list)>1)
  warn_len = length(rm_dp_id)
  if(check & warn_len>0){
    if(warn_len==B & show){
      new_obs_graph_list = new_obs_graph_list
      new_dir_graph_list = new_dir_graph_list
    }else if(warn_len==B & !show){
      new_obs_graph_list = new_obs_graph_list
      new_dir_graph_list = NULL
      warning = TRUE
    }else{
    new_obs_graph_list = new_obs_graph_list[-rm_dp_id]
    new_dir_graph_list = new_dir_graph_list[-rm_dp_id]
    }
  }
  if(warn_len==0){warn_len=NULL}
  obs_graph_mean = apply(simplify2array(new_obs_graph_list), 1:2, mean)
  obs_graph_sd = apply(simplify2array(new_obs_graph_list), 1:2, sd)
  obs_graph_pval = pnorm(-abs(obs_graph_mean/obs_graph_sd))*2
  obs_graph_pval[which(obs_graph_pval==0)] = 5e-8
  colnames(obs_graph_mean)=colnames(obs_graph_sd)=colnames(obs_graph_pval)=dp_list$trait_vec[keep_trait_id]
  rownames(obs_graph_mean)=rownames(obs_graph_sd)=rownames(obs_graph_pval)=dp_list$trait_vec[keep_trait_id]
  
  pval_3dmat = simplify2array(new_obs_pval_graph_list)
  pval_3dmat = pval_3dmat[keep_trait_id,keep_trait_id,]
  dims = dim(pval_3dmat)
  twoDimMat <- matrix(pval_3dmat,prod(dims[1:2]), dims[3])
  twoDimMat[!apply(twoDimMat,1,sum)<1e-100,] -> twoDimMat
  lambda = eigen(cor(t(twoDimMat)))$value
  Me = ceiling(length(lambda) - sum((lambda>1)*(lambda-1)))
  
  if(!is.null(new_dir_graph_list)){
    dir_graph_mean = apply(simplify2array(new_dir_graph_list), 1:2, mean)
    dir_graph_sd = apply(simplify2array(new_dir_graph_list), 1:2, sd)
    dir_graph_pval = pnorm(-abs(dir_graph_mean/dir_graph_sd))*2
    dir_graph_pval[which(dir_graph_pval==0)] = 5e-8
    colnames(dir_graph_mean)=colnames(dir_graph_sd)=colnames(dir_graph_pval)=dp_list$trait_vec[keep_trait_id]
    rownames(dir_graph_mean)=rownames(dir_graph_sd)=rownames(dir_graph_pval)=dp_list$trait_vec[keep_trait_id]
    eigen_dir = eigen(dir_graph_mean)
    if(any(abs(eigen_dir$values)>1)){warning=TRUE}
  }else{
    dir_graph_mean = dir_graph_pval = dir_graph_sd = NULL
  }
  dp_list = list()
  dp_list$obs_graph_list = new_obs_graph_list
  dp_list$dir_graph_list = new_dir_graph_list
  dp_res = list()
  dp_res$obs_graph_mean = obs_graph_mean
  dp_res$obs_graph_sd = obs_graph_sd
  dp_res$obs_graph_pval = obs_graph_pval
  dp_res$dir_graph_mean = dir_graph_mean
  dp_res$dir_graph_sd = dir_graph_sd
  dp_res$dir_graph_pval = dir_graph_pval
  dp_res$Me = Me
  return(list(dp_res=dp_res,dp_list=dp_list,warning=warning,warn_len=warn_len,conv_len=conv_len,total_B=total_B))
}

# G_mean and G_pval are matrices output from subset_Graph_d1(), e.g. obs_graph_mean and obs_graph_pval, or dir_graph_mean and dir_graph_pval
# Me is the number of effective tests output from subset_Graph_d1(), significance threshold is by default 0.05/Me 
# thres1 is the secondary p-value threshold with light-colored edges
plot_graph <- function(G_mean,G_pval,Me,thres1=0.05,Bonferroni=FALSE){
  set.seed(1)
  if(Bonferroni==FALSE & !is.null(Me)){
    thres = 0.05/Me
  }else{
    thres = 0.05/(nrow(G_mean)^2-nrow(G_mean))
  }
  G_pval[is.na(G_pval)] = 1
  if(thres1==0){
    A1 = G_pval<thres
  }else{
    A1 = G_pval<thres1
  }
  A = G_mean
  A[!A1] = 0
  diag(A) = 0
  A2 = G_pval
  A2[!A1] = 0
  diag(A2) = 0
  p2 = graph_from_adjacency_matrix(adjmatrix=A2,mode='directed',weighted=T)
  p = graph_from_adjacency_matrix(adjmatrix=A,mode='directed',weighted=T)
  E(p)$color <- ifelse(E(p)$weight<0,rgb(178,34,34,1,max=255),rgb(79,201,120,1,max=255))
  E(p)$color <- dplyr::case_when(
    E(p)$weight<0 & E(p2)$weight<thres ~ rgb(178,34,34,255,max=255),
    E(p)$weight<0 & E(p2)$weight<thres1 ~ rgb(178,34,34,100,max=255),
    E(p)$weight>0 & E(p2)$weight<thres ~ rgb(79,201,120,255,max=255),
    E(p)$weight>0 & E(p2)$weight<thres1 ~ rgb(79,201,120,100,max=255))
  E(p)$weight <- 1
  edge_width = dplyr::case_when(E(p2)$weight<thres ~ 3,
#                                E(p2)$weight<thres2 ~ 2,
                                E(p2)$weight<thres1 ~ 1)
  V(p)$color <- ifelse(V(p)$name %in% c("CAD","Stroke","T2D","Asthma","AD","AF"),'lightblue','orange') 
  # lightblue color for disease nodes, and orange for risk factor nodes
  return(list(p=p,edge_width=edge_width))
}


