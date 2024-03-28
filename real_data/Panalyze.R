source('/home/panwei/lin00374/cML/graph/cMLgraph.R')
#source('/home/panwei/lin00374/cML/cML_O_old.R')
library(dplyr)

array_id = Sys.getenv('SLURM_ARRAY_TASK_ID')

dat = get(load('/home/panwei/lin00374/cML/graph/17trait/data/AD+AF+11r4d.RData'))
IV_list = get(load('/home/panwei/lin00374/cML/graph/17trait/data/AD+AF+11r4d_IV.RData'))

dp_list = Graph_Perturb(b_mat = dat$b_mat,
                        se_mat = dat$se_mat,
                        n_vec = dat$n_vec,
                        rho_mat = dat$P,
#                        rho_mat = diag(17),
                        IV_list = IV_list,
                        R_list = dat$R_list,
                        random_start = 5,
                        num_pert = 10, seed=array_id,curse=FALSE)


save(dp_list,file=paste0('dp_list_seed',array_id,'.RData'))






