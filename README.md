# debias-covariate-shift

The functions for running low dimensional simulations are found in debiased linear.R, plugin.R and the functions for high dimensional simulations are found in debiased lasso.R and plugin high dim.R in their respective folders.

These functions return a dataframe with 4 columns. Each observation records the result of one trial. The four columns V1,V2,V3,V4 record the test statistic, the standard deviation, whether the test rejected the null and the normalized test statistic respectively.
