library(TwoSampleMR)
library(data.table)
library(dplyr)
source('/home/panwei/lin00374/TWAS_XSquare/Feb4_2021/allele_qc.R')
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

TSMR.ID.Vec =
  c("ebi-a-GCST002216","ebi-a-GCST002222","ebi-a-GCST002223",
    "ieu-a-89","ieu-a-835","ieu-a-1083","ieu-b-39",
    "ieu-b-38","ebi-a-GCST000568","ieu-b-25","ieu-b-73","ebi-a-GCST005195",
    "ebi-a-GCST005838","ieu-a-26","ebi-a-GCST006862","ebi-a-GCST006414","Jansen2019")



m = length(TSMR.ID.Vec) 
traits = c('TG','LDL','HDL','Height','BMI','BW','DBP','SBP','FG','Smoke','Alcohol','CAD','Stroke','T2D','Asthma','AF','AD')

for(i in 1:(m)){
    t0 = paste0(traits[i],' <- extract_instruments("',TSMR.ID.Vec[i],'",clump=F)')
    eval(parse(text=t0))
}

AD_all = fread('/home/panwei/lin00374/TWAS_XSquare/JansenGWAS/AD_sumstats_Jansenetal_2019sept.txt.gz')
AD_all = AD_all %>% select(chr=CHR,pos=BP,SNP,effect_allele.outcome=A1,other_allele.outcome=A2,beta.outcome=BETA,se.outcome=SE,pval.outcome=P,eaf.outcome=EAF,samplesize.outcome=Nsum) 
AD = AD_all %>% filter(pval.outcome<5e-08)
all.IV = list()

k = 1

for(exp.index in 1:(m-1)){
    exp_id = TSMR.ID.Vec[exp.index]
  for(out.index in (exp.index+1):length(traits)){
      if(out.index!=length(traits)){
    out_id = TSMR.ID.Vec[out.index]
    
    
    # preprocess data ---------------------------------------------------------
    set.seed(1)
    t0 = paste0("union.IV = union(",traits[exp.index],"$SNP,",traits[out.index],"$SNP)")
    eval(parse(text=t0))
    
    exposure_dat <- extract_outcome_data(union.IV, exp_id, proxies = 0)
    exposure_dat = exposure_dat[which(!duplicated(exposure_dat$SNP)),]
    exposure_dat = exposure_dat[order(exposure_dat$chr,exposure_dat$pos),]
    
    outcome_dat <- extract_outcome_data(union.IV, out_id, proxies = 0)
    outcome_dat = outcome_dat[which(!duplicated(outcome_dat$SNP)),]
    outcome_dat = outcome_dat[order(outcome_dat$chr,outcome_dat$pos),]
      }else{
          t0 = paste0("union.IV = union(",traits[exp.index],"$SNP,","AD$SNP)")
          eval(parse(text=t0))
          exposure_dat <- extract_outcome_data(union.IV, exp_id, proxies = 0)
          exposure_dat = exposure_dat[which(!duplicated(exposure_dat$SNP)),]
          exposure_dat = exposure_dat[order(exposure_dat$chr,exposure_dat$pos),]
          outcome_dat <- AD_all %>% filter(SNP %in% union.IV)
          outcome_dat = outcome_dat[which(!duplicated(outcome_dat$SNP)),]
          outcome_dat = outcome_dat[order(outcome_dat$chr,outcome_dat$pos),]
          outcome_dat$id.outcome = TSMR.ID.Vec[out.index]
      }
    
    common.IV = intersect(exposure_dat$SNP,outcome_dat$SNP)
    
    exposure_dat_new = 
      exposure_dat[which(is.element(exposure_dat$SNP,common.IV)),]
    outcome_dat_new = 
      outcome_dat[which(is.element(outcome_dat$SNP,common.IV)),]
    
    exp.out.flip = 
      allele.qc(a1 = exposure_dat_new$effect_allele.outcome,
                a2 = exposure_dat_new$other_allele.outcome,
                ref1 = outcome_dat_new$effect_allele.outcome,
                ref2 = outcome_dat_new$other_allele.outcome)
    if(sum(exp.out.flip$flip)>0)
    {
      outcome_dat_new$beta.outcome[exp.out.flip$flip] = 
        -outcome_dat_new$beta.outcome[exp.out.flip$flip]
    }
    exposure_dat_new = exposure_dat_new[exp.out.flip$keep,]
    outcome_dat_new = outcome_dat_new[exp.out.flip$keep,]
    
    ### do the clumping
    pval.clump = -log10(exposure_dat_new$pval.outcome) - log10(outcome_dat_new$pval.outcome)
    pval.clump = 10^(-pval.clump)
    
    data.clump = 
      data.frame(rsid = exposure_dat_new$SNP,
                 pval = pval.clump)
    
    clump.result = 
      ieugwasr::ld_clump(dat = data.clump)

    clump.result = merge(clump.result,exposure_dat_new,by.x='rsid',by.y='SNP') %>% select(rsid,chr,pval)

    #all.IV = union(all.IV,clump.result$rsid)
    all.IV[[k]] = as.character(clump.result$rsid)
    k = k+1
  }
}

save(all.IV,file='./data/AD+AF+11r4d_IV.RData')

all.IV = unique(unlist(all.IV))
#for(i in 1:m){
#        exposure_dat <- extract_outcome_data(all.IV, TSMR.ID.Vec[i], proxies = 0)
#        save(exposure_dat,file=paste0('./data/',TSMR.ID.Vec[i],'.RData'))
#}

AD_dat = AD_all %>% filter(SNP %in% all.IV)
AD_dat$id.outcome = 'Jansen2019'
#save(AD_dat,file=paste0('./data/AD.RData'))

all_exposure_dat <- extract_outcome_data(all.IV, TSMR.ID.Vec, proxies = 0)

#all_trait_pivot = tidyr::pivot_wider(data = all_exposure_dat, 
#            id_cols = SNP, 
#            names_from = id.outcome, 
#            values_from = c("effect_allele.outcome","other_allele.outcome"))
#
#check = apply(all_trait_pivot,1,function(a1){any(as.character(a1[3:12])!=as.character(a1[2]),na.rm=T)})
#print(sum(check))

all_exposure_dat <- plyr::rbind.fill(all_exposure_dat,AD_dat)


all_dup_df = all_exposure_dat[duplicated(all_exposure_dat %>% select(SNP,id.outcome)) |
                    duplicated(all_exposure_dat %>% select(SNP,id.outcome),fromLast=TRUE),]
dup_df = all_exposure_dat[duplicated(all_exposure_dat %>% select(SNP,id.outcome)),]

for(i in 1:nrow(dup_df)){
    dup_temp_df = all_exposure_dat %>% filter(SNP==dup_df$SNP[i] & id.outcome==dup_df$id.outcome[i])
    compare_dup_df = all_exposure_dat %>% filter(SNP==dup_df$SNP[i] & id.outcome!=dup_df$id.outcome[i])
    keep_dup = which(dup_temp_df$effect_allele.outcome==compare_dup_df$effect_allele.outcome[1] &
          dup_temp_df$other_allele.outcome==compare_dup_df$other_allele.outcome[1] )
    all_exposure_dat = all_exposure_dat %>% filter(!(SNP==dup_df$SNP[i] & id.outcome==dup_df$id.outcome[i] &
                                effect_allele.outcome == dup_temp_df[-keep_dup,]$effect_allele.outcome &
                                other_allele.outcome == dup_temp_df[-keep_dup,]$other_allele.outcome))
}

all_trait_pivot = tidyr::pivot_wider(data = all_exposure_dat,
            id_cols = SNP,
            names_from = id.outcome,
            values_from = c("effect_allele.outcome","other_allele.outcome"))

m = length(unique(all_exposure_dat$id.outcome))
trait_id = 2:(m+1)
not_check_SNP = all_trait_pivot$SNP
all_checked_SNP = setdiff(all_trait_pivot$SNP,not_check_SNP)
j = 1

while(length(not_check_SNP)>0){
    ref_trait_id = trait_id[j]
    alt_trait_id = trait_id[-j]
    ### ref_trait_id NA
    to_check_SNP_ind = which(all_trait_pivot$SNP %in% not_check_SNP) 
    not_check_SNP = intersect(all_trait_pivot$SNP[which(is.na(all_trait_pivot[,ref_trait_id]))],not_check_SNP)
    check = apply(all_trait_pivot[to_check_SNP_ind,],1,function(a1){any(as.character(a1[alt_trait_id])!=as.character(a1[ref_trait_id]),na.rm=T)})
    flip_SNP = all_trait_pivot[to_check_SNP_ind[check],]$SNP
    t_name = names(all_trait_pivot)[alt_trait_id]
    for(i in flip_SNP){
        t = as.character(all_trait_pivot %>% filter(SNP==i))
        flip_trait = unlist(lapply(strsplit(t_name[which(t[alt_trait_id]!=t[ref_trait_id])],split='outcome_'),function(a){a[2]}))
        all_exposure_dat[which(all_exposure_dat$SNP == i & all_exposure_dat$id.outcome%in% flip_trait),]$beta.outcome = -all_exposure_dat[which(all_exposure_dat$SNP == i & all_exposure_dat$id.outcome%in% flip_trait),]$beta.outcome
        seqinr::swap(all_exposure_dat[which(all_exposure_dat$SNP == i & all_exposure_dat$id.outcome%in% flip_trait),]$effect_allele.outcome,all_exposure_dat[which(all_exposure_dat$SNP == i & all_exposure_dat$id.outcome%in% flip_trait),]$other_allele.outcome)
    }
    all_trait_pivot = tidyr::pivot_wider(data = all_exposure_dat,
            id_cols = SNP,
            names_from = id.outcome,
            values_from = c("effect_allele.outcome","other_allele.outcome"))
    check = apply(all_trait_pivot[to_check_SNP_ind,],1,function(a1){any(as.character(a1[alt_trait_id])!=as.character(a1[ref_trait_id]),na.rm=T)})
    print(sum(check))
    all_checked_SNP = setdiff(all_trait_pivot$SNP,not_check_SNP) 
    print(length(not_check_SNP))
    j = j+1
}




for(i in 1:m){
        exposure_dat <- all_exposure_dat %>% dplyr::filter(id.outcome==TSMR.ID.Vec[i])
        save(exposure_dat,file=paste0('./data/',TSMR.ID.Vec[i],'.RData'))
}



