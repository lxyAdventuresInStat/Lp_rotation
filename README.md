# Lp_rotation
This is a repository for conducting experiments in Rotation to Sparse Loadings using Lp Losses and Related Inference Problems http://arxiv.org/abs/2206.02263.

cl510_oblique_rotation_function.R contains all the auxiliary functions, such as IRGP algorithm, proximal gradient descent algorithm, etc. It is the source for all experiments. The file is the same in all three folders.

#Study 1
Table 1 is produced mainly by cl513obl_experiment.15S.R, cl513obl_experiment.30S.R, where the 15S, 30S indicating the setting and the two files only differ in parameters, not in experiments. There is no need to adjust sample size N for it.
Table 2 is produced by cl505A_experiment.15S.R, cl505A_experiment.30S.R, where the 15S, 30S indicating the setting. We need to adjust N to 400, 800, 1600 for the full table.
The results of the above experiments are saved in Rdata files, and we use res.outputfinal.Rmd to produce the tables and Figure 3.Rmd is a rmarkdown file.

#Study 2
Table 3 is produced by cl505Counter_experiment.15S.R. There is no need to adjust sample size N for it.

#Big-Five Personality Test
Table 4-7 is produced by big5_r_code_final.Rmd. Rmd is a rmarkdown file.
