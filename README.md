# Brief introduction to using GraphMRcML

## Step 1. Data preparation
First, to analyze the relationship among $N$ traits, you need to prepare your pre-processed data, including 
* `dat`: a list of 5 elements:
  * The first two elements are `b_mat` and `se_mat`, which are both $m\times n$ matrices storing the GWAS estimates (and corresponding standard errors) of $M$ SNPs on $N$ traits.
  * The third element `n_vec` is a vector of length $N$, storing the sample size for each trait.
  * The fourth element `P` is an $N\times N$ correlation matrix, storing the correlations between the GWAS estimates among the $N$ traits, which can be obtained by, for example, bivariate LDSC regression.
  * The fifth element `R_list` is a *list* of LD information among the $M$ SNPs, which is splitted into LD blocks. Each element in `R_list` has two elements, `snp` is a vector of SNP names (rsid), and `R`
 is the LD correlation matrix of SNPs in `snp`.

* `IV_list`: a list of $N\times (N-1)/2$ elements, which correspond to a total of $N\times (N-1)/2$ pairwise MR analyses among $N$ traits. Specifically, each element in `IV_list` stores the *union* of
SNPs (rsid) to be used in the *bivariate* MR analysis for each pair of traits.

The sample codes for preprocessing the data are given in `real_data/real_data.R` and `real_data/process.R`.

Briefly, in `real_data.R`, we obtain IVs for each pair of traits by searching GWAS significant SNPs with trait A and GWAS significant SNPs with trait B, and do LD-clumping on the union.
Then this LD-clumped union of SNPs is stored in `IV_list`.
We then extract all $M$ SNPs in `IV_list` for all $N$ traits and perform data harmonizing.

In `process.R`, we construct the LD correlation matrix among the $M$ SNPs by first splitting them into LD blocks, and call `ieugwasr::ld_matrix` to obtain LD matrix.

**Note**
1. Some elements in `b_mat` and `se_mat` can be NA, as long as those SNPs will not be used in the bivariate analyses involving that trait, which is guaranteed by the use of `IV_list`.
2. Based on the preprocessing code, we expect the order in `IV_list` should be pair 1-2, pair 1-3, ..., pair 1-N, pair 2-3, pair 2-4, ..., pair (N-1)-N.

## Step 2. Run GraphMRcML
After preparing the above two datasets, we can implement GraphMRcML. As GraphMRcML relies on data perturbation (bootstrap) for statistical inference, we can run the analysis in parallel. 
`real_data/Panalyze.R` shows an example to run the analysis with 10 perturbations (`num_pert = 10`) in each job, then you can submit, e.g. 200 jobs, to obtain a total of 2000 perturbations.

## Step 3. Summarize and visualize results
After Step 2, you first need to combine results from all 200 jobs (see `real_data/combine_list.R`), then you can summarize the result using `subset_Graph_d1` and visualize it with `plot_graph` in `R/helper.R`.
In particular, `subset_Graph_d1` returns a list storing the final total-effect graph (`obs_graph_mean`,`obs_graph_sd`,`obs_graph_pval`) 
and the final direct-effect graph after network deconvolution (`dir_graph_mean`,`dir_graph_sd`,`dir_graph_pval`). 

It also stores `warn_len` which is the number of direct-effect graphs with a spectral value greater than 1 (which violates the key assumption in ND), 
and `conv_len` which is the number of failed convergence of the algorithm to ensure the diagonal elements of the resulting direct-effect network are zeros.
In general, if the two numbers are large, the resulting direct-effect graph may not be reliable. 
One strategy is to reduce the number of traits and/or highly correlated traits in the model.

## Reference
Lin, Z., Xue, H., & Pan, W. (2023). 
[**Combining Mendelian randomization and network deconvolution for inference of causal networks with GWAS summary data**](https://doi.org/10.1371/journal.pgen.1010762). PLoS genetics, 19(5), e1010762. 



