what does this do?
did we have it in the R code? will rho ever be -1?

    cov = max(cov, -1.0 + 1e-12)  # guard against division by zero at rho=-1

covariance.py:132

======

does the code apply climate extreme value theory, what are the paraemters in Jusutus code or in our code that decide the number of clusters, k and i want the plots to be coloured and not dotted pointed like in the [text](docs/figures/edc_clusters_map.pdf) and [text](docs/figures/edc_risk_es_map.pdf) [text](docs/figures/lec_clusters_map.pdf) 

======

I want to reproduce the same graphs as in the extremes paper figure nine page 729 application to climate data and use the data used there so the figure nine says Results of the clustering algorithms EDC (panel a) and LEC (panel b) applied to precipitation
data from a historical run of the global climate model AWI-ESM-1-1LR. Parameters used are ν = 5,
α = 1.0, ε = 5. The number of clusters is 104 for the EDC algorithm and 24 for the LEC algorithm.
It was determined using as cut-off threshold the empirical 30%-quantile of the calculated pairwise
dissimilarities... but this time we focus on precipitations instead of temperatuee is this good, how can i link my risk metrics to the precidipations between the clusters generated, 


=====

why do we print plots to pdf files? isn't it better to store as image files and then refer to them in the latex code?