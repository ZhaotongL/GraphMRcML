library(dplyr)
library(tidyverse)
# TG 24097068: ebi-a-GCST002216
# LDL 24097068: ebi-a-GCST002222
# HDL 24097068: ebi-a-GCST002223
# Height 25282103: ieu-a-89
# BMI 25673413: ieu-a-835
# BW 27680694: ieu-a-1083
# DBP 30224653: ieu-b-39
# SBP 30224653: ieu-b-38
# FG 20081858: ebi-a-GCST000568
# Smoke 30643251: ieu-b-25
# Alcohol 30643251: ieu-b-73
# CAD 29212778: ebi-a-GCST005195
# Stoke 29531354: ebi-a-GCST005838
# T2D 22885922: ieu-a-26
# Asthma 29273806: ebi-a-GCST006862
# AF: ebi-a-GCST006414
# AD: /home/panwei/lin00374/TWAS_XSquare/JansenGWAS/AD_sumstats_Jansenetal_2019sept.txt.gz

traits = c('TG','LDL','HDL','Height','BMI','BW','DBP','SBP','FG','Smoke','Alcohol','CAD','Stroke','T2D','Asthma','AF','AD')

TSMR.ID.Vec =
  c("ebi-a-GCST002216","ebi-a-GCST002222","ebi-a-GCST002223",
    "ieu-a-89","ieu-a-835","ieu-a-1083","ieu-b-39",
    "ieu-b-38","ebi-a-GCST000568","ieu-b-25","ieu-b-73","ebi-a-GCST005195",
    "ebi-a-GCST005838","ieu-a-26","ebi-a-GCST006862","ebi-a-GCST006414","Jansen2019")

n_vec = c()
for(t in 1:length(traits)){
    t0 = paste0(traits[t],"=get(load('./data/",TSMR.ID.Vec[t],".RData'))")
    eval(parse(text=t0))
    t1 = paste0("n_t=mean(",traits[t],"$samplesize.outcome)")
    eval(parse(text=t1))
    n_vec = c(n_vec,n_t)
}
n_vec = ceiling(n_vec)

for(t in traits){
    t0 = paste0(t,"=",t,"%>% select(SNP,chr,pos,A1=effect_allele.outcome,A2=other_allele.outcome,b=beta.outcome,se=se.outcome)")
    eval(parse(text=t0))
}

all_trait = paste0(traits,collapse=',')
t1 = paste0("list(",all_trait,")%>% reduce(full_join, by = c('SNP','chr','pos','A1','A2')) -> all_disease")
eval(parse(text=t1))
all_disease$chr = as.numeric(all_disease$chr)
all_disease$pos = as.numeric(all_disease$pos)
all_disease %>% arrange(chr,pos) -> all_disease
print(nrow(all_disease))

b_mat = as.matrix(all_disease[,seq(from=6,to=ncol(all_disease),by=2)])
rownames(b_mat) = all_disease$SNP
colnames(b_mat) = traits
se_mat = as.matrix(all_disease[,seq(from=7,to=ncol(all_disease),by=2)])
rownames(se_mat) = all_disease$SNP
colnames(se_mat) = traits

all_SNP = all_disease %>% select(SNP,chr,pos,A1,A2)
ldBlocks = read.table('/home/panwei/shared/zhaotong_share/ldblocks-eur-hg19.txt',header=TRUE)
matList = list()
matListInd = 1
for(i in 1:nrow(ldBlocks)){
  curChr = substring(ldBlocks[i,1],4)
  curMin = ldBlocks[i,2]
  curMax = ldBlocks[i,3]
  curSnps = all_SNP[all_SNP$pos > curMin & all_SNP$pos <= curMax & all_SNP$chr == curChr,]
  if(nrow(curSnps) > 0){
    mat_res = list()
    if(nrow(curSnps)==1){
        R_flipped = matrix(1,nrow=1)
        colnames(R_flipped) = rownames(R_flipped) = curSnps$SNP
    }else{
        R = ieugwasr::ld_matrix(variants=curSnps$SNP)
        R_snp_df = do.call(rbind.data.frame,strsplit(colnames(R),split='_'))
        colnames(R_snp_df) = c('SNP','A1','A2')
        R_mydat_df = merge(R_snp_df,curSnps,by='SNP',sort=F)
        R_toflip_ind = which(with(R_mydat_df,A1.x==A2.y & A2.x==A1.y))
        R_flipped = R
        for(j in R_toflip_ind){
            R_flipped[j,] = -R_flipped[j,]
            R_flipped[,j] = -R_flipped[,j]
        }
        colnames(R_flipped) = rownames(R_flipped) = as.character(R_snp_df$SNP)
    }
    mat_res$snp = curSnps$SNP
    mat_res$R = R_flipped
    matList[[matListInd]] = mat_res
    matListInd = matListInd+1
  }
}
print(matListInd)

## rho_mat
full_rho_mat = get(load('/home/panwei/lin00374/cML/graph/vcffiles/ldsc_res/rho_thres_0.1.RData'))
rho_mat = full_rho_mat[traits,traits]

dat = list()
dat$b_mat = b_mat
dat$se_mat = se_mat
dat$n_vec = n_vec
dat$P = rho_mat
dat$R_list = matList

save(dat,file='./data/AD+AF+11r4d.RData')
