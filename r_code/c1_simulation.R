args = commandArgs(trailingOnly=TRUE)
path_out <- paste0("~/maxstable/complete/",args[1])

source(paste0("~/maxstable/complete/parameters.R"))
source(paste0("~/maxstable/complete/functions.R"))

a_sim_exp_ns <- sim_expt_2d_nonstat(X, Y, df_true, alpha_true, a_matrix,
                                    b_matrix, g_matrix, n_sim=n_sim)

saveRDS(a_sim_exp_ns, paste0(path_out,"/res/simulation.rds"))