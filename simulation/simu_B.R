## change sample size n!
n_gwas = 300000
##
source('/home/panwei/lin00374/cML/graph/cMLgraph.R')
source('/home/panwei/lin00374/cML/cML_O_old.R')
load('./data/b_simu_true.RData')
load('./data/6trait.RData')
IV_list = get(load('./data/6trait_IV.RData'))
load('./data/6trait_screen_res.RData')

array_id = Sys.getenv('SLURM_ARRAY_TASK_ID')

set.seed(array_id)
traits=colnames(dat$b_mat)
b_mat = dat$b_mat; pval_mat = dat$pval_mat;
G = matrix(0,nrow=6,ncol=6)
rownames(G)=colnames(G)=traits
G['BMI','FG'] = 0.08;
G['BMI','TG'] = 0.23;
G['BMI','CAD'] = 0.22;
G['BMI','AF'] = 0.32;
G['LDL','CAD'] = 0.42;
G['FG','CAD'] = 0.33;
G['CAD','AF'] = 0.17;
G['AF','CAD'] = 0.10;
b = matrix(0,nrow=nrow(b_mat),ncol=ncol(b_mat))
colnames(b)=traits
BMI_iv = which(pval_mat[,'BMI']<5e-8)
LDL_iv = which(pval_mat[,'LDL']<5e-8)
FG_iv = which(pval_mat[,'FG']<5e-8)
TG_iv = which(pval_mat[,'TG']<5e-8)
CAD_iv = which(pval_mat[,'CAD']<5e-8)
AF_iv = which(pval_mat[,'AF']<5e-8)

b[BMI_iv,'BMI'] = b_mat[BMI_iv,'BMI']
b[BMI_iv,'FG'] = G['BMI','FG']*b[BMI_iv,'BMI']
b[BMI_iv,'TG'] = G['BMI','TG']*b[BMI_iv,'BMI']
b[BMI_iv,'CAD'] = ((0.22+0.33*0.08+0.1*0.32)/(1-0.1*0.17)) *b[BMI_iv,'BMI']
b[BMI_iv,'AF'] = ((0.32+0.17*0.22+0.17*0.33*0.08)/(1-0.1*0.17))*b[BMI_iv,'BMI']
b[LDL_iv,'LDL'] = b_mat[LDL_iv,'LDL'] 
b[LDL_iv,'CAD'] = (0.42/(1-0.1*0.17))*b[LDL_iv,'LDL']
b[LDL_iv,'AF'] = (0.17*0.42)/(1-0.1*0.17)*b[LDL_iv,'LDL']
b[FG_iv,'FG'] = b_mat[FG_iv,'FG']
b[FG_iv,'CAD'] = (0.33/(1-0.1*0.17))*b[FG_iv,'FG']
b[FG_iv,'AF'] = (0.17*0.33)/(1-0.1*0.17)*b[FG_iv,'FG']
b[TG_iv,'TG'] = b_mat[TG_iv,'TG'] 
b[CAD_iv,'CAD'] = (1/(1-0.1*0.17))*b_mat[CAD_iv,'CAD']
b[CAD_iv,'AF'] = 0.17/(1-0.1*0.17)*b[CAD_iv,'CAD']
b[AF_iv,'AF'] = (1/(1-0.1*0.17))*b_mat[AF_iv,'AF']
b[AF_iv,'CAD'] = (0.1/(1-0.1*0.17))*b[AF_iv,'AF']

rownames(b)=rownames(b_mat)

se_mat = dat$se_mat
n_vec = rep(n_gwas,6)
se_mat[which(!is.na(se_mat))] = 1/sqrt(n_gwas)
b_mat_simu = Generate_Perturb(b_mat=b,
                            se_mat=se_mat,
                            n_vec=n_vec,rho_mat=dat$rho_mat,DP_mat_list=screen_res$DP_mat_list)

Graph_Perturb_simu <- function(b_mat,se_mat,screen_res,
                          n_vec,rho_mat,
                          num_pert=100,random_start=10,seed=0,trait_vec=NULL,curse=F){
    set.seed(seed)
    if(curse){
        t = -qnorm(sig.cutoff/2)
    }else{t=0}
  if(is.null(trait_vec)){trait_vec=colnames(screen_res$b_mat)}
  out_list = pbmcapply::pbmclapply(1:num_pert,function(i){    b_mat_dp = Generate_Perturb(b_mat=b_mat,
                                se_mat=se_mat,
                                n_vec=n_vec,rho_mat=rho_mat,DP_mat_list=screen_res$DP_mat_list);
                                Graph_Estimate(b_mat=b_mat_dp,
                                     se_mat=se_mat,
                                     n_vec=n_vec,rho_mat=rho_mat,
                                     IJ_snp_list=screen_res$IJ_snp_list,
                                     t=0,
                                     random_start=random_start)},ignore.interactive=TRUE)
  obs_graph_list = lapply(out_list,function(x){x$obs_graph})
  obs_graph_se_list = lapply(out_list,function(x){x$obs_graph_se})
  obs_graph_pval_list = lapply(out_list,function(x){x$obs_graph_pval})

  out = list()
  out$obs_graph_list = obs_graph_list
  out$obs_graph_se_list = obs_graph_se_list
  out$obs_graph_pval_list = obs_graph_pval_list
  out$trait_vec = trait_vec
  return(out)
}

dp_list = Graph_Perturb_simu(b_mat = b_mat_simu,
                        se_mat = se_mat,
                        screen_res = screen_res,
                        n_vec = n_vec,
                        rho_mat = dat$P,
                        random_start = 5,
                        num_pert = 200, seed=array_id,curse=FALSE)
save(dp_list,file=paste0('./result2/simu',array_id,'.RData'))




