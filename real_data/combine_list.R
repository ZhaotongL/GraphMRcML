dp_obs_list = dp_dir_list = dp_obs_pval_list = list()
for(i in 1:200){
    load(paste0('dp_list_seed',i,'.RData'))    
    dp_obs_list = append(dp_obs_list,dp_list$obs_graph_list)
    dp_obs_pval_list = append(dp_obs_pval_list,dp_list$obs_graph_pval_list)
}
dp_obs_list = lapply(dp_obs_list,function(G){diag(G)=0;G})
dp_obs_pval_list = lapply(dp_obs_pval_list,function(G){diag(G)=0;G})
traits = dp_list$trait_vec
dp_list = list()
dp_list$obs_graph_list = dp_obs_list
dp_list$obs_graph_pval_list = dp_obs_pval_list
dp_list$trait_vec = traits 


B = length(dp_list$obs_graph_list)
save(dp_list,file=paste0('dp_list_B',B,'.RData'))
