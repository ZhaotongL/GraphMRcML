Graph_Screen <- function(b_mat,se_mat,n_vec,IV_list,R_list,rho_mat,c_vec=rep(1,length(n_vec)),sig.cutoff=5e-08){
  n_trait = length(n_vec)
  N_combination = n_trait * (n_trait - 1) / 2
  if(length(IV_list)!=N_combination){stop("The length of IV_list must be equal to N_combination!")}
  m_block = length(R_list) ## LD blocks
  if(sum(unlist(lapply(R_list,function(x){nrow(x$R)})))!=nrow(b_mat)){stop("LD matrix must contain all SNPs used in the analysis!")}

  IJ_snp_list = vector("list", N_combination)
  k = 1
### Screening for IVs for each pair of traits ###
  ## i:trait 1; j:trait 2
  for(i in 1:(n_trait-1)){
    for(j in (i+1):n_trait){
      SNP_used = IV_list[[k]]
      b_X = b_mat[,i]
      b_Y = b_mat[,j]
      se_X = se_mat[,i]
      se_Y = se_mat[,j]
      n_X = n_vec[i]
      n_Y = n_vec[j]
      keep_ind = intersect(which(!is.na(b_X)),which(!is.na(b_Y)))
      pvalue.X = pnorm(-abs(b_X/se_X))*2
      pvalue.Y = pnorm(-abs(b_Y/se_Y))*2
      cor_X = b_X / sqrt(b_X^2 + (n_X-2)*se_X^2)
      cor_Y = b_Y / sqrt(b_Y^2 + (n_Y-2)*se_Y^2)

      ind_X = which(pvalue.X<(sig.cutoff))
      ind_Y = which(pvalue.Y<(sig.cutoff))
      # With Screening
      intersect.ind.X.Y = intersect(ind_X,ind_Y)
      ind_X_new = setdiff(ind_X,
                          intersect.ind.X.Y[(abs(cor_X)[intersect.ind.X.Y])<
                                              (abs(cor_Y)[intersect.ind.X.Y])])
      ind_Y_new = setdiff(ind_Y,
                          intersect.ind.X.Y[(abs(cor_X)[intersect.ind.X.Y])>
                                              (abs(cor_Y)[intersect.ind.X.Y])])
      ind_X_new = intersect(ind_X_new,keep_ind)
      ind_Y_new = intersect(ind_Y_new,keep_ind)
      ind_X_new_snp = rownames(b_mat)[ind_X_new]
      ind_Y_new_snp = rownames(b_mat)[ind_Y_new]

      ind_X_new1 = ind_X_new[is.element(ind_X_new_snp,SNP_used)]
      ind_Y_new1 = ind_Y_new[is.element(ind_Y_new_snp,SNP_used)]
      IJ_snp_list[[k]]$ind_i_new = ind_X_new1
      IJ_snp_list[[k]]$ind_j_new = ind_Y_new1
      k = k + 1
    }
  }
### End of screening 1 ###
### Screening for LD of b_mat ###
  DP_mat_list = vector("list", m_block)
  for(i in 1:m_block){
      R = R_list[[i]]$R
      match_order = match(rownames(b_mat),rownames(R))
      match_order = match_order[!is.na(match_order)]
      R = R[match_order,match_order,drop=FALSE]
      snp_in_block = R_list[[i]]$snp
      PXR = kronecker(rho_mat,R)
      eigen_decomp = eigen(PXR)
      D = diag(sqrt(zapsmall(eigen_decomp$value,digits=6)))
      DP_mat_list[[i]]$V = eigen_decomp$vector %*% D
      DP_mat_list[[i]]$snp = rownames(R)
  }
### End of screening 2 ###

### Screnning for inflation factor (c) from LDSC ###
  for(i in 1:n_trait){
      se_mat[,i] = max(1,sqrt(c_vec[i])) * se_mat[,i]
  }
### End of screnning 3 ###
  out = list()
  out$IJ_snp_list = IJ_snp_list
  out$DP_mat_list = DP_mat_list
  out$b_mat = b_mat
  out$se_mat = se_mat
  return(out)
}

Generate_Perturb <- function(b_mat,se_mat,n_vec,rho_mat,DP_mat_list){
  n_trait = length(n_vec)
  m_used = nrow(b_mat)
  e_mat_dp = matrix(0,ncol=n_trait,nrow=m_used)
  for(i in 1:length(DP_mat_list)){
      V = DP_mat_list[[i]]$V
      X = rnorm(nrow(V))
      e_vec = V %*% X
      e_mat = matrix(e_vec,ncol=n_trait,byrow=FALSE)
      snp_ind = which(is.element(rownames(b_mat),DP_mat_list[[i]]$snp))
      e_mat_dp[snp_ind,] = se_mat[snp_ind,,drop=FALSE] * e_mat
  }
  b_mat_dp = b_mat + e_mat_dp

  return(b_mat_dp)
}

Graph_Estimate <- function(b_mat,se_mat,n_vec,rho_mat,IJ_snp_list,t,random_start=10){
  n_trait = length(n_vec)
  obs_graph = matrix(1,nrow=n_trait,ncol=n_trait)
  obs_graph_pval = obs_graph_se = matrix(0,nrow=n_trait,ncol=n_trait)
  k = 1
  for(i in 1:(n_trait-1)){
    for(j in (i+1):n_trait){
      ind_i_new = IJ_snp_list[[k]]$ind_i_new
      ind_j_new = IJ_snp_list[[k]]$ind_j_new
### for perturbed data, need to select again! ###
      ind_i_new = ind_i_new[abs(b_mat[ind_i_new,i]/se_mat[ind_i_new,i])>t]
      ind_j_new = ind_j_new[abs(b_mat[ind_j_new,j]/se_mat[ind_j_new,j])>t]
###
      rho_ij = rho_mat[i,j]
      ItoJ_cML_O_res = mr_cML_O(b_exp=b_mat[ind_i_new,i],
                                 b_out=b_mat[ind_i_new,j],
                                 se_exp=se_mat[ind_i_new,i],
                                 se_out=se_mat[ind_i_new,j],
                                 n = min(n_vec[i],n_vec[j]),
                                 rho = rho_ij,
                                 random_start = random_start,t=t)
      JtoI_cML_O_res = mr_cML_O(b_exp=b_mat[ind_j_new,j],
                                     b_out=b_mat[ind_j_new,i],
                                     se_exp=se_mat[ind_j_new,j],
                                     se_out=se_mat[ind_j_new,i],
                                     n = min(n_vec[i],n_vec[j]),
                                     rho = rho_ij,
                                     random_start = random_start,t=t)

     obs_graph[i,j] = ItoJ_cML_O_res$BIC_theta
     obs_graph[j,i] = JtoI_cML_O_res$BIC_theta
     obs_graph_se[i,j] = ItoJ_cML_O_res$BIC_se
     obs_graph_se[j,i] = JtoI_cML_O_res$BIC_se
     obs_graph_pval[i,j] = ItoJ_cML_O_res$BIC_p
     obs_graph_pval[j,i] = JtoI_cML_O_res$BIC_p
      k = k + 1

    }
  }
  dir_graph = obs_graph %*% solve(diag(n_trait)+obs_graph)
#  eigen_og = eigen(dir_graph)
#  if(any(abs(eigen_og$values)>1)){warning("Some eigenvalues' absolute values are greater than 1!")}

  out = list()
  out$obs_graph = obs_graph
  out$obs_graph_se = obs_graph_se
  out$obs_graph_pval = obs_graph_pval
  out$dir_graph = dir_graph

  return(out)

}

Graph_Perturb <- function(b_mat,se_mat,n_vec,rho_mat,IV_list,R_list,c_vec=rep(1,length(n_vec)),
                          sig.cutoff=5e-08,num_pert=100,random_start=10,seed=0,trait_vec=NULL,curse=F){
    set.seed(seed)
    if(curse){
        t = -qnorm(sig.cutoff/2)
    }else{t=0}

  screen_res = Graph_Screen(b_mat=b_mat,se_mat=se_mat,
                            n_vec=n_vec,IV_list=IV_list,R_list=R_list,rho_mat=rho_mat,c_vec=c_vec,
                            sig.cutoff=sig.cutoff)
  if(is.null(trait_vec)){trait_vec=colnames(screen_res$b_mat)}
  out_list = pbmcapply::pbmclapply(1:num_pert,function(i){    b_mat_dp = Generate_Perturb(b_mat=screen_res$b_mat,
                                se_mat=screen_res$se_mat,
                                n_vec=n_vec,rho_mat=rho_mat,DP_mat_list=screen_res$DP_mat_list);
                                Graph_Estimate(b_mat=b_mat_dp,
                                     se_mat=screen_res$se_mat,
                                     n_vec=n_vec,rho_mat=rho_mat,
                                     IJ_snp_list=screen_res$IJ_snp_list,
                                     t=t,
                                     random_start=random_start)},ignore.interactive=TRUE)
  obs_graph_list = lapply(out_list,function(x){x$obs_graph})
  obs_graph_se_list = lapply(out_list,function(x){x$obs_graph_se})
  obs_graph_pval_list = lapply(out_list,function(x){x$obs_graph_pval})
  dir_graph_list = lapply(out_list,function(x){x$dir_graph})

  out = list()
  out$obs_graph_list = obs_graph_list
  out$obs_graph_se_list = obs_graph_se_list
  out$obs_graph_pval_list = obs_graph_pval_list
  out$dir_graph_list = dir_graph_list
  out$trait_vec = trait_vec
  return(out)
}

Graph_summarize <- function(Graph_Perturb_list){
    obs_graph_mean = apply(simplify2array(Graph_Perturb_list$obs_graph_list), 1:2, mean)
    obs_graph_sd = apply(simplify2array(Graph_Perturb_list$obs_graph_list), 1:2, sd)
    obs_graph_pval = pnorm(-abs(obs_graph_mean/obs_graph_sd))*2

    pval_3dmat = simplify2array(Graph_Perturb_list$obs_graph_pval_list)
    dims = dim(pval_3dmat)
    twoDimMat <- matrix(pval_3dmat,prod(dims[1:2]), dims[3])
    twoDimMat[!apply(twoDimMat,1,sum)<1e-100,] -> twoDimMat
    lambda = eigen(cor(t(twoDimMat)))$value
    Me = ceiling(length(lambda) - sum((lambda>1)*(lambda-1)))

    dir_graph_mean = apply(simplify2array(Graph_Perturb_list$dir_graph_list), 1:2, mean)
    dir_graph_sd = apply(simplify2array(Graph_Perturb_list$dir_graph_list), 1:2, sd)
    dir_graph_pval = pnorm(-abs(dir_graph_mean/dir_graph_sd))*2

    colnames(obs_graph_mean)=colnames(obs_graph_sd)=colnames(obs_graph_pval)=Graph_Perturb_list$trait_vec
    colnames(dir_graph_mean)=colnames(dir_graph_sd)=colnames(dir_graph_pval)=Graph_Perturb_list$trait_vec
    rownames(obs_graph_mean)=rownames(obs_graph_sd)=rownames(obs_graph_pval)=Graph_Perturb_list$trait_vec
    rownames(dir_graph_mean)=rownames(dir_graph_sd)=rownames(dir_graph_pval)=Graph_Perturb_list$trait_vec

    eigen_dir = eigen(dir_graph_mean)
    if(any(abs(eigen_dir$values)>1)){warning("Some eigenvalues' absolute values are greater than 1!")}


    out = list()
    out$obs_graph_mean = obs_graph_mean
    out$obs_graph_sd = obs_graph_sd
    out$obs_graph_pval = obs_graph_pval
    out$dir_graph_mean = dir_graph_mean
    out$dir_graph_sd = dir_graph_sd
    out$dir_graph_pval = dir_graph_pval
    out$Me = Me

    return(out)
}

subset_Graph <- function(dp_list,keep_trait,alpha=1,check=TRUE){
    keep_trait_id = is.element(dp_list$trait_vec,keep_trait)

    new_obs_graph_list = lapply(dp_list$obs_graph_list,function(G) {alpha*G[keep_trait_id,keep_trait_id]})
    new_dir_graph_list = lapply(new_obs_graph_list, function(G_obs) {G_obs %*% solve(diag(nrow(G_obs))+G_obs)})
    max_eigen_list = lapply(new_dir_graph_list, function(G_dir) {max(abs(eigen(G_dir)$values))})
    rm_dp_id = which(unlist(max_eigen_list)>1)
    if(check & length(rm_dp_id)>0){
        new_dir_graph_list = new_dir_graph_list[-rm_dp_id]
    }

    obs_graph_mean = apply(simplify2array(new_obs_graph_list), 1:2, mean)
    obs_graph_sd = apply(simplify2array(new_obs_graph_list), 1:2, sd)
    obs_graph_pval = pnorm(-abs(obs_graph_mean/obs_graph_sd))*2

    pval_3dmat = simplify2array(dp_list$obs_graph_pval_list)
    pval_3dmat = pval_3dmat[keep_trait_id,keep_trait_id,]
    dims = dim(pval_3dmat)
    twoDimMat <- matrix(pval_3dmat,prod(dims[1:2]), dims[3])
    twoDimMat[!apply(twoDimMat,1,sum)==0,] -> twoDimMat
    lambda = eigen(cor(t(twoDimMat)))$value
    Me = ceiling(length(lambda) - sum((lambda>1)*(lambda-1)))

    dir_graph_mean = apply(simplify2array(new_dir_graph_list), 1:2, mean)
    dir_graph_sd = apply(simplify2array(new_dir_graph_list), 1:2, sd)
    dir_graph_pval = pnorm(-abs(dir_graph_mean/dir_graph_sd))*2

    colnames(obs_graph_mean)=colnames(obs_graph_sd)=colnames(obs_graph_pval)=dp_list$trait_vec[keep_trait_id]
    colnames(dir_graph_mean)=colnames(dir_graph_sd)=colnames(dir_graph_pval)=dp_list$trait_vec[keep_trait_id]
    rownames(obs_graph_mean)=rownames(obs_graph_sd)=rownames(obs_graph_pval)=dp_list$trait_vec[keep_trait_id]
    rownames(dir_graph_mean)=rownames(dir_graph_sd)=rownames(dir_graph_pval)=dp_list$trait_vec[keep_trait_id]

    eigen_dir = eigen(dir_graph_mean)
    if(any(abs(eigen_dir$values)>1)){warning("Some eigenvalues' absolute values are greater than 1!")}
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
    return(list(dp_res=dp_res,dp_list=dp_list))
}



