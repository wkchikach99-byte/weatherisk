############################### Functions #####################################
dist_x <- function(x1,x2) {x1-x2}
dist_y <- function(y1,y2) {y1-y2}
# conversion between grid values (x,y), matrix coordinates [i,j] and
# R-internal enumeration of values
grid_number <- function(i,j) {
  if(i<1 | j<1 | i>nrow_grid | j>ncol_grid) {stop(paste0("Index out of Bounds"))}
  return((j-1)*nrow_grid+i)
}
number_grid <- function(n) {
  if(n<1 | n>nrow_grid*ncol_grid) {stop(paste0("Index out of Bounds"))}
  return(c((n-1)%%nrow_grid+1, floor((n-1)/nrow_grid+1)))
}
koord_num <- function(x, y) {
  return(which.min((X-x)^2+(Y-y)^2))
}
number_koord <- function(n) {
  return(c(X[n],Y[n]))
}
array_apply <- function(array, fun) {
  if(dim(array)[1]!=nrow_grid | dim(array)[2] != ncol_grid) {
    stop("Wrong array dimensions")
  }
  return(matrix(sapply(1:nrow_grid, function(i) {
    sapply(1:ncol_grid, function(j) {
      fun(array[i,j,])
    })}), nrow=nrow_grid, ncol=ncol_grid))
}
rad <- function(x) {x*pi/180}
deg <- function(x) {180*x/pi}

plot_map <- function(matrix, von=min(matrix,na.rm=T), bis=max(matrix,na.rm=T),
                     rev=T, main="",print=T, legend=T, margin=0, truncate=F,
                     cr=colorRampPalette(brewer.pal(11,"RdBu")[6:11])(21),
                     contours=0*X) {
  n_col=length(cr)
  if(truncate) {
    matrix[which(matrix<von)]<-von
    matrix[which(matrix>bis)]<-bis
  }
  if(rev){cr = rev(cr)}
  att <- seq(von, bis, length.out = n_col)
  if(margin==0) {
    gr = exp.grid
  } else {
    gr <- expand.grid(y=y_ax[(margin+1):(length(y_ax)-margin)],
                      x=x_ax[(margin+1):(length(x_ax)-margin)])
  }
  if(legend==FALSE) {
    colorkey = FALSE
  }
  else {
    leg <- rep(F,floor(n_col/4)); leg[1] <- T
    colorkey = list(at=c(von,att),labels=list(at=c(von,att[leg]),
                     labels=round(c(1000,att[leg]),2),
                     col=c("white",rep("black",length(att[leg])))))
    # invisible 1000-label to ensure the legend has the same width in all plots
  }
  if(max(contours)>0) panel="draw_contours" else panel=panel.levelplot
  p <- levelplot(matrix~x*y, data=gr,panel=panel,at=seq(von, bis,
                                            length.out=n_col),
                 col.regions = cr,main=main, colorkey=colorkey,
                 contours=contours,
                 par.settings =  list(layout.heights = list(
                   top.padding = 0,
                   main.key.padding = 0,
                   key.axis.padding = 0,
                   axis.xlab.padding = 0,
                   xlab.key.padding = 0,
                   key.sub.padding = 0),
                   layout.widths = list(
                     left.padding = 0,
                     key.ylab.padding = 0,
                     ylab.axis.padding = 0,
                     axis.key.padding = 0,
                     right.padding = -1*(legend==F))))
  if(print==T) {print(p)} else {return(p)}
}
draw_contours <- function(..., contours) {
  gr <- expand.grid(y=y_ax[(0+1):(length(y_ax)-0)],
                    x=x_ax[(0+1):(length(x_ax)-0)])
  panel.levelplot(...)
  for(i in min(contours):max(contours))
    panel.contourplot(gr[,2],gr[,1],1*(matrix(contours,nrow=resolution)==i)
                    ,1:length(contours),contour=T,region=F,at=c(0.5,1.5))
}

dtdiff <- function(x, df) {
  return( gamma((df+1)/2)/gamma(df/2)/sqrt(df*pi)*
            (-df-1)/2*(1+x^2/df)^(-(df+1)/2-1)*2*x/df )
}
pairwise_density_summand <- function(z1,z2, x,y ,df, alpha,a,b,g) {
  cv <- cov_fkt_2d(x,y, alpha,a,b,g)
  c <- sqrt(1-cv*cv)/sqrt(df+1)
  m1 <- ((z2/z1)^(1/df)-cv)/c
  m2 <- ((z1/z2)^(1/df)-cv)/c
  dm1 = dnorm(m1); dm2 = dnorm(m2) ; pm1 = pnorm(m1) ; pm2 = pnorm(m2)
  # log-density:  log(dz1(V)dz2(V)-dz1z2(V)) - V
  return( log(
    (-pt(m1,df+1)/z1/z1 - dt(m1,df+1)*z2^(1/df)*z1^(-1/df-2)/c/df +
       dt(m2,df+1)*z1^(1/df-1)*z2^(-1/df-1)/c/df)*
      ( -pt(m2,df+1)/z2/z2 - dt(m2,df+1)*z1^(1/df)*z2^(-1/df-2)/c/df +
          dt(m1,df+1)*z2^(1/df-1)*z1^(-1/df-1)/c/df)
    +
      (dt(m1,df+1)*z1^(-1/df-2)*z2^(1/df-1) +
         dt(m2,df+1)*z2^(-1/df-2)*z1^(1/df-1) +
         dt(m1,df+1)*z1^(-1/df-2)*z2^(1/df-1)/df +
         dt(m2,df+1)*z2^(-1/df-2)*z1^(1/df-1)/df +
         dtdiff(m1,df+1)*z1^(-2/df-2)*z2^(2/df-1)/c/df +
         dtdiff(m2,df+1)*z2^(-2/df-2)*z1^(2/df-1)/c/df)/c/df 
  )  -  (pt(m1,df+1)/z1 + pt(m2,df+1)/z2)
  )
}

# set the values at margin distance to 0
crop_matrix <- function(matrix, margin) {
  matrix <- 0*X + matrix
  matrix <- matrix[(1+margin):(length(y_ax)-margin),
                   (1+margin):(length(x_ax)-margin)]
  return(matrix)
}

cov_fkt_2d <- function(x,y, alpha=1, a=1, b=1, g=0) {
  # a>0: kleinere Diagonale
  # b>0: Differenz zw. größerer und kleinerer Diagonale
  # g in [-pi/2, pi/2]: Richtung der gr. Diagonale
  # the version in the comments is incorrect: a and (a+b) switched
  # return(exp(-sqrt(t(c(x,y))%*%
  #                    (t(matrix(c(cos(g)/(a+b),
  # -sin(g)/a, sin(g)/(a+b), cos(g)/a),nrow=2))%*%
  #                   matrix(c(cos(g)/(a+b),
  # -sin(g)/a, sin(g)/(a+b), cos(g)/a),nrow=2))%*%
  #                  (c(x,y)))^alpha))
  # return(exp(-sqrt(t(c(x,y))%*%
  #                    matrix(c(sin(g)*sin(g)/(a+b)/(a+b)+cos(g)*cos(g)/a/a,
  #                             sin(g)*cos(g)*(1/(a+b)/(a+b)-1/a/a),
  #                             sin(g)*cos(g)*(1/(a+b)/(a+b)-1/a/a),
  #                             cos(g)*cos(g)/(a+b)/(a+b)+sin(g)*sin(g)/a/a
  #                    ),nrow=2)%*%
  #                    c(x,y))^alpha))
  return(exp(-sqrt(x*x*(sin(g)*sin(g)/a/a+cos(g)*cos(g)/(a+b)/(a+b)) + 
                     2*x*y*sin(g)*cos(g)*(-1/a/a+1/(a+b)/(a+b)) +
                     y*y*(cos(g)*cos(g)/a/a+sin(g)*sin(g)/(a+b)/(a+b)))^alpha))
}
plot_cov <- function(alpha, a,b,g) {
  plot_map(cov_fkt_2d(X,Y,alpha=alpha,a,b,g), von=0, bis=1, rev=F)
}

cov_to_ec <- function(df, cov) {
  return(ifelse(cov==1, 1, 
                2*pt((1-cov)/((df+1)^(-1/2)*sqrt(1-cov*cov) ), df+1)))
}
ec_to_cov <- function(df, ec) {
  ec <- min(ec, cov_to_ec(df, 0))
  uniroot((function (x) cov_to_ec(df, x) - ec), lower = 0, upper = 1)[1]$root
}
cov_fkt_2d_nonstat2 <- function(x,y, alpha=1,
                                a1=1,b1=1,g1=0,
                                a2=1,b2=1,g2=0) {
  # a>0: kleinere Diagonale
  # b>0: Differenz zw. größerer und kleinerer Diagonale
  # g in [0°, 180°]: Richtung der gr. Diagonale
  m1 <- matrix(c(sin(g1)*sin(g1)/a1/a1+cos(g1)*cos(g1)/(a1+b1)/(a1+b1),
                 sin(g1)*cos(g1)*(-1/a1/a1+1/(a1+b1)/(a1+b1)),
                 sin(g1)*cos(g1)*(-1/a1/a1+1/(a1+b1)/(a1+b1)),
                 cos(g1)*cos(g1)/a1/a1+sin(g1)*sin(g1)/(a1+b1)/(a1+b1)
  ),nrow=2)
  m2 <- matrix(c(sin(g2)*sin(g2)/a2/a2+cos(g2)*cos(g2)/(a2+b2)/(a2+b2),
                 sin(g2)*cos(g2)*(-1/a2/a2+1/(a2+b2)/(a2+b2)),
                 sin(g2)*cos(g2)*(-1/a2/a2+1/(a2+b2)/(a2+b2)),
                 cos(g2)*cos(g2)/a2/a2+sin(g2)*sin(g2)/(a2+b2)/(a2+b2)
  ),nrow=2)
  om <- solve((solve(m1)+solve(m2))/2)
  # 1/sqrt(det(m1)) = a1*(a1+b1)
  return(min(1,sqrt(det(om)*a1*(a1+b1)*a2*(a2+b2))*
               exp(-sqrt(t(c(x,y))%*%om%*%c(x,y))^alpha)))
}
plot_ec_nonstat <- function(x,y, df, alpha, a_matrix,b_matrix,g_matrix,
                             local_approx=F) {
  wm <- koord_num(x,y)
  if(local_approx==F) {
    plot_map(sapply(1:n_grid, function(i) {
      cov_to_ec(df, cov_fkt_2d_nonstat2(X[i]-x,Y[i]-y,alpha=alpha,
                          a1=a_matrix[wm], b1=b_matrix[wm], g1=g_matrix[wm],
                          a2=a_matrix[i],b2=b_matrix[i],g2=g_matrix[i]))})
      , von=1, bis=2, rev=T)
  } 
  else {
    plot_map(sapply(1:n_grid, function(i) {
      cov_to_ec(df, cov_fkt_2d_nonstat2(X[i]-x,Y[i]-y,alpha=alpha,
                          a1=a_matrix[wm], b1=b_matrix[wm], g1=g_matrix[wm],
                          a2=a_matrix[wm],b2=b_matrix[wm],g2=g_matrix[wm]))})
      , von=1, bis=2, rev=T)
  }
}

sim_expt_2d <- function(X,Y, df, alpha, a, b, g, n_sim=1) {
  colmax <- function(Z) sapply(1:ncol(Z), function(t) {max(Z[,t])})
  c_df <- 2^(1-df/2)*sqrt(pi)*1/(gamma((df+1)/2)) 
  Cmax <- c_df*(qnorm(0.999)^df)
  lenX <- length(X)
  # covariance matrix for fBm
  lenX <- length(X)
  cov_matrix <- matrix(sapply(1:n_grid, function(i) { sapply(1:n_grid, function(j) {
    cov_fkt_2d(X[i]-X[j],Y[i]-Y[j],alpha,a,b,g)
  })}), nrow=n_grid)
  cov_chol <- chol(cov_matrix)
  if(n_sim > 1) {sim_all <- array(dim=c(nrow_grid, ncol_grid, n_sim))}
  for(nn in 1:n_sim) {
    e_i <- poiss <- numeric(0)
    cnt<-1
    Z <- sim <- numeric(0)
    # fbm <- matrix( t(cov_chol)%*%rnorm(length(X),0,1), nrow=length(ax))
    while(cnt>0) {
      e_i <- c(e_i, rexp(1,1))
      poiss <- c(poiss, 1/sum(e_i))
      # simulate gaussian process
      gauss_pr <- t(cov_chol)%*%rnorm(n_grid,0,1)
      # calculate spectral process
      Z <- rbind(Z, c_df*pmax(0, gauss_pr)^df)
      sim <- rbind(sim, colmax(Z*poiss))
      if (Cmax*poiss[cnt] <= min(colmax(sim))) {break}
      cnt <- cnt+1
    }
    if(n_sim > 1) {
      sim_all[,,nn] <- matrix(colmax(sim), nrow=nrow_grid)
    }
  }
  if(n_sim > 1) {return(sim_all) }
  else {return(matrix(colmax(sim), nrow=nrow_grid)) }
}
sim_expt_2d_nonstat <- function(X,Y, df, alpha,
                                a_matrix, b_matrix, g_matrix, n_sim=1) {
  colmax <- function(Z) sapply(1:ncol(Z), function(t) {max(Z[,t])})
  c_df <- 2^(1-df/2)*sqrt(pi)*1/(gamma((df+1)/2)) 
  Cmax <- c_df*(qnorm(0.999)^df)
  print("Berechnung Cov-Matrix")
  cov_matrix <- matrix(
    sapply(1:n_grid, function(i) { sapply(1:n_grid, function(j) {
      cov_fkt_2d_nonstat2(X[i]-X[j],Y[i]-Y[j], alpha,
                          a_matrix[i],b_matrix[i],g_matrix[i],
                          a_matrix[j],b_matrix[j],g_matrix[j])
    })}), nrow=n_grid)
  cov_chol <- chol(cov_matrix)
  if(n_sim > 1) {sim_all <- array(dim=c(nrow_grid,ncol_grid,n_sim))}
  print(paste0("Simulating ",n_sim," max-stable processes"))
  for(nn in 1:n_sim) {
    if(nn%%10==0) {
      print(nn)
    }
    e_i <- poiss <- numeric(0)
    cnt<-1
    Z <- sim <- numeric(0)
    # fbm <- matrix( t(cov_chol)%*%rnorm(length(X),0,1), nrow=length(ax))
    while(cnt>0) {
      e_i <- c(e_i, rexp(1,1))
      poiss <- c(poiss, 1/sum(e_i))
      # simulate gaussian process
      gauss_pr <- t(cov_chol)%*%rnorm(length(X),0,1)
      # calculate spectral process
      Z <- rbind(Z, c_df*pmax(0, gauss_pr)^df)
      sim <- rbind(sim, colmax(Z*poiss))
      if (Cmax*poiss[cnt] <= min(colmax(sim))) {break}
      cnt <- cnt+1
    }
    if(n_sim > 1) {
      sim_all[,,nn] <- matrix(colmax(sim), nrow=nrow_grid)
    }
  }
  if(n_sim > 1) {return(sim_all) }
  else {return(matrix(colmax(sim), nrow=nrow_grid)) }
}

c_extrcoeff_matrix <- function(a_sim, madogram=FALSE) {
  # first calculate the ranks of each gev distributed time series
  rank_matrix <- matrix(nrow=n_grid, ncol=n_sim)
  print("Calculating Rank Matrix")
  for(s in 1:n_grid) {
    i <- number_grid(s)[1]
    j <- number_grid(s)[2]
    for (k in 1:n_sim) {
      rank_matrix[s,k] <- sum(a_sim[i,j,k]<=a_sim[i,j,1:n_sim])
    }
  }
  print("Calculating Extremal Coefficients")
  extrcoeff_matrix <- matrix(0, ncol=n_grid, nrow=n_grid)
  # calculate upper triangular matrix of extremal coefficients (minus 1)
  for(i in 1:(n_grid-1)) { for(j in (i+1):n_grid) {
    v = mean(abs(rank_matrix[i,1:n_sim]-rank_matrix[j,1:n_sim]))/(2*(n_sim+1))
    if(madogram)
    {
      extrcoeff_matrix[i,j] = v
    }
    else
    {
      extrcoeff_matrix[i,j] = min(1, (1+2*v)/(1-2*v) -1 )
    }
  }}
  return(extrcoeff_matrix+t(extrcoeff_matrix))
}
c_real_extrcoeff_matrix <- function(df,alpha,a_matrix,b_matrix,g_matrix) {
  extrcoeff_matrix <- matrix(0, nrow=n_grid, ncol=n_grid)
  for(i in 1:(n_grid-1)) { for(j in (i+1):n_grid) {
    extrcoeff_matrix[i,j] = cov_to_ec(df, cov_fkt_2d_nonstat2(X[i]-X[j],
                                                              Y[i]-Y[j],alpha,
                                                              a_matrix[i],b_matrix[i],g_matrix[i],
                                                              a_matrix[j],b_matrix[j],g_matrix[j])) - 1
  }}
  return(extrcoeff_matrix+t(extrcoeff_matrix))
}

pairwise_density_optim <- function(z, df, alpha, X ,Y , print=T,max_dist=0,
                                   lower_bounds=c(0.01,0.01), upper_bounds=c(15,15),ensemble=3) {
  n_grid <- length(X)
  if(length(dim(z))==3) {
    z <- sapply(1:(dim(z)[3]), function(i) {as.vector(z[,,i])})
  }
  if(nrow(z) != n_grid) stop("nrow(z) must equal the nrow(distance_matrix)")
  ilist <- rep(1,n_grid)%*%t((1:n_grid)) ; ilist <- ilist[lower.tri(ilist)]
  jlist <- (1:n_grid)%*%t(rep(1,n_grid)) ; jlist <- jlist[lower.tri(jlist)]
  Xlist <- rep(X[ilist]-X[jlist],each=dim(z)[2])
  Ylist <- rep(Y[ilist]-Y[jlist],each=dim(z)[2])
  
  zilist <- as.vector(t(z[ilist,]))
  zjlist <- as.vector(t(z[jlist,]))     
  
  if(max_dist>0) {
    sel <- which(Xlist*Xlist+Ylist*Ylist<=max_dist*max_dist)
    zilist <- zilist[sel]
    zjlist <- zjlist[sel]
    Xlist  <-  Xlist[sel]
    Ylist  <-  Ylist[sel]
  }
  if(length(zilist) > 0)
  {
    llh <- function(par) {
      return(sum(pairwise_density_summand(zilist,zjlist,Xlist,Ylist, df, alpha,
                                          par[1],par[2],par[3])))
    }
    
    parameters_lower_bound <- c( lower_bounds[1], lower_bounds[2], -pi/2)
    parameters_upper_bound <- c( upper_bounds[1], upper_bounds[2], pi/2  )
    parnames <- c("a", "b", "g")
    parscale <- (parameters_upper_bound - parameters_lower_bound)/100
    
    ensemble <- ensemble
    
    set <- matrix(parameters_lower_bound,ensemble,
                  length(parameters_lower_bound),byrow=T) + 
      maximinLHS(ensemble,length(parameters_lower_bound))%*%diag(
        parameters_upper_bound-parameters_lower_bound)
    
    for(i in 1:ensemble) {
      if(print) {print(paste("Lauf",i))}
      fit_<-optim(set[i,], fn=llh, method=c("L-BFGS-B"),
                  lower=parameters_lower_bound, upper=parameters_upper_bound,
                  control=list(fnscale=-1, parscale=parscale, maxit=10000))
      # if +- pi/2 reached, continue from the other side
      if(abs(fit_$par[3])==pi/2) {
        fit_<-optim(c(fit_$par[1],fit_$par[2],-fit_$par[3]), fn=llh, method=c("L-BFGS-B"),
                    lower=parameters_lower_bound, upper=parameters_upper_bound,
                    control=list(fnscale=-1, parscale=parscale, maxit=10000))
      }
      if(print) {
        print("MLE Estimates:")
        names(fit_$par) <- parnames
        print(round(fit_$par, 2))
        print(round(fit_$value, 2))
      }
      if(i == 1) {
        fit <- fit_
      } else {
        if(-fit_$value < -fit$value) {
          fit <- fit_
        }
      }
    }
    return(as.vector(fit$par))
  }
  else
  {
    return(c(0,0,0))
  }
}
pairwise_density_optim_local <- function(a_sim, df,alpha, x,y ,abstand=3, print=T, ensemble=1,
                                         lower_bounds = c(0.01,0.01), upper_bounds = c(15,15)) {
  xy_pos <- which.min(dist_x(x,X)^2+dist_y(y,Y)^2)[1]
  xs <- rep(-abstand:abstand, (2*abstand+1))
  ys <- rep(-abstand:abstand, each=(2*abstand+1))
  which.dist <- which(sqrt(xs*xs+ys*ys)<=abstand & (xs*xs+ys*ys)>0)
  xs <- xs[which.dist] ; ys <- ys[which.dist]
  sel.x <- number_grid(xy_pos)[2] + xs
  sel.y <- number_grid(xy_pos)[1] + ys
  in.bounds <- which(sel.y>0 & sel.y<=nrow_grid & sel.x>0 & sel.x<=ncol_grid)
  sel.x <- sel.x[in.bounds]
  sel.y <- sel.y[in.bounds]
  sel.grid <- sapply(1:length(in.bounds), function(i) {grid_number(sel.y[i], sel.x[i])})
  
  n_sel <- length(sel.x)
  
  # compare all observations in the surrounding (in zlist1) with the observations
  # at the given grid points (in zlist2)
  zlist1 <- as.vector(sapply(1:n_sel, function(i) {a_sim[number_grid(sel.grid[i])[1],number_grid(sel.grid[i])[2],] }))
  zlist2 <- rep(a_sim[number_grid(xy_pos)[1],number_grid(xy_pos)[2],], n_sel)
  
  Xlist <- rep(X[sel.grid]-X[xy_pos],each=n_sim)
  Ylist <- rep(Y[sel.grid]-Y[xy_pos],each=n_sim)
  
  llh <- function(par) {
    return(sum(pairwise_density_summand(zlist1,zlist2,Xlist,Ylist, df, alpha,
                                        par[1],par[2],par[3])))
  }
  
  max_num_additional_calculations = 5
  
  parameters_lower_bound <- c(lower_bounds, -pi/2)
  parameters_upper_bound <- c(upper_bounds, pi/2  )
  parnames <- c( "a", "b", "g")
  parscale <- (parameters_upper_bound - parameters_lower_bound)/100
  
  set <- matrix(parameters_lower_bound,ensemble+max_num_additional_calculations,
                length(parameters_lower_bound),byrow=T) + 
    maximinLHS(ensemble+max_num_additional_calculations,length(parameters_lower_bound))%*%
    diag(parameters_upper_bound-parameters_lower_bound)
  
  num_calc = 1
  num_more_calc_done = 0
  while(num_calc<=ensemble) {
    fit_<-optim(set[num_calc+num_more_calc_done,], fn=llh, method=c("L-BFGS-B"),
                lower=parameters_lower_bound, upper=parameters_upper_bound,
                control=list(fnscale=-1, parscale=parscale, maxit=10000))
    # if +- pi/2 reached, continue from the other side
    if(abs(fit_$par[3])==pi/2) {
      fit_<-optim(c(fit_$par[1],fit_$par[2],-fit_$par[3]), fn=llh, method=c("L-BFGS-B"),
                  lower=parameters_lower_bound, upper=parameters_upper_bound,
                  control=list(fnscale=-1, parscale=parscale, maxit=10000))
    }
    if(print) {
      toprint <- c(num_calc, round(fit_$par,4), round(fit_$value,2))
      names(toprint) <- c("Run Nr.", parnames, "Value")
      print(toprint)
    }
    if(num_calc+num_more_calc_done == 1) {
      fit <- fit_
    } else {
      if(-fit_$value < -fit$value) {
        fit <- fit_
      }
    }
    if(min(abs(fit_$par-parameters_lower_bound))<0.01 ||
       min(abs(fit_$par-parameters_upper_bound))<0.01) {
      # rerun analysis if upper or lower bound reached, but only max. 3 times
      if(num_more_calc_done < max_num_additional_calculations) {
        num_more_calc_done=num_more_calc_done+1
        if(print) {print("Re-run because par limits reached")}
        next;
      }
    }
    num_calc=num_calc+1
  }
  names(fit$par) <- parnames
  return(fit$par)
}

crop_local_estimates <- function(local_estimates_matrix, margin) {
  return(cbind(as.vector(crop_matrix(local_estimates_matrix[,1], margin)),
               as.vector(crop_matrix(local_estimates_matrix[,2], margin)),
               as.vector(crop_matrix(local_estimates_matrix[,3], margin))))
}
plot_local_estimates<- function(x,y, extrcoeff_matrix, local_estimates_matrix)
{
  num <- koord_num(x,y)
  toprint <- c(local_estimates_matrix[num,2]+local_estimates_matrix[num,1],
               local_estimates_matrix[num,1],deg(local_estimates_matrix[num,3]))
  names(toprint) = c("a+b","a","g")
  print(round(toprint, 2))
  p1 <- plot_map(extrcoeff_matrix[,num]+1, von=1, bis=2, legend=F,
                 main="Estimation from Data",print=F)
  p2 <- plot_map(sapply(1:n_grid, function(i) {cov_to_ec(df_true,cov_fkt_2d(
    dist_x(X[i],x),dist_y(Y[i],y),alpha=alpha_true,
    a=local_estimates_matrix[num,1],
    b=local_estimates_matrix[num,2],
    g=local_estimates_matrix[num,3]))})
    , von=1, bis=2, legend=F, main="Estimation from local Estimates",print=F)
  print(p1, split=c(1,1,2,1), more=TRUE)
  print(p2, split=c(2,1,2,1))
}
calc_distance_ellipses <- function(local_estimates_matrix, res=21) {
  xs <- rep(seq(-1,1,length.out=res), res)
  ys <- rep(seq(-1,1,length.out=res), each=res)
  wh <- which(((xs^2+ys^2)<=res^2) & (ys > 0) | (xs > 0))
  # only half-circle needs to be considered bco symmetry 
  xs <- xs[wh]; ys <- ys[wh]
  ell1 <- ell2 <- 0*xs
  dist_matrix <- matrix(0, nrow=nrow(local_estimates_matrix),
                        ncol=nrow(local_estimates_matrix))
  cnt=0
  for(i in 1:(nrow(local_estimates_matrix)-1)) {
    for(j in (i+1):(nrow(local_estimates_matrix))) {
      cnt = cnt + 1
      if((cnt%%100000==0)) print(cnt/1000000)
      mx <- max(local_estimates_matrix[i,1]+local_estimates_matrix[i,2],
                local_estimates_matrix[j,1]+local_estimates_matrix[j,2])
      ell1 <- 1*(cov_fkt_2d(xs,ys,alpha=1,
                            a=local_estimates_matrix[i,1]/mx,
                            b=local_estimates_matrix[i,2]/mx,
                            g=local_estimates_matrix[i,3])>exp(-1))
      ell2 <- 1*(cov_fkt_2d(xs,ys,alpha=1,
                            a=local_estimates_matrix[j,1]/mx,
                            b=local_estimates_matrix[j,2]/mx,
                            g=local_estimates_matrix[j,3])>exp(-1))
      # the value of alpha doesn't matter here as we're only interested
      # in the shape of the ellipses
      if(max(sum(ell1),sum(ell2))==0) {
        dist_matrix[j]=1
      }
      else {
        dist_matrix[i,j] <- 1-(sum(ell1*ell2)+1/2)/(sum(ell1+ell2-ell1*ell2)+1/2)
      }
      # +1/2 because of the point (0,0) which is always included 
    }
  }
  return(100*(dist_matrix + t(dist_matrix)))
}
plot_distance_ellipses <- function(local_estimates_matrix, num,
                                   resolution=21) {
  xs <- rep(seq(-1,1,length.out=resolution), resolution)
  ys <- rep(seq(-1,1,length.out=resolution), each=resolution)
  wh <- which(((xs^2+ys^2)<=1) & (ys > 0) | (xs > 0))
  # only half-circle needs to be considered bco symmetry 
  xs <- xs[wh]; ys <- ys[wh]
  ell1 <- ell2 <- 0*xs
  dist_matrix <- numeric(n_grid)
  for(j in (1:n_grid)) {
    mx <- max(local_estimates_matrix[num,1]+local_estimates_matrix[num,2],
              local_estimates_matrix[j,1]+local_estimates_matrix[j,2])
    ell1 <- 1*(cov_fkt_2d(xs,ys,alpha=1,
                          a=local_estimates_matrix[num,1]/mx,
                          b=local_estimates_matrix[num,2]/mx,
                          g=local_estimates_matrix[num,3])>exp(-1))
    ell2 <- 1*(cov_fkt_2d(xs,ys,alpha=1,
                          a=local_estimates_matrix[j,1]/mx,
                          b=local_estimates_matrix[j,2]/mx,
                          g=local_estimates_matrix[j,3])>exp(-1))
    # the value of alpha doesn't matter here as we're only interested
    # in the shape of the ellipses
    if(max(sum(ell1),sum(ell2))==0) {
      dist_matrix[j]=100
    }
    else {
      dist_matrix[j] <- 1-(sum(ell1*ell2)+1/2)/(sum(ell1+ell2-ell1*ell2)+1/2)
      # +1/2 because of the point (0,0) which is always included 
    }
  }
  plot_map(100*dist_matrix, von=0, bis=100)
}

smooth_local_estimates <- function(local_estimates, smoothing_dist) {
  matrix_smoothing <- function(matr, smoothing_dist, angle=F) {
    center_angle <- function(angle_list, center) {
      # convert elements of angle_list to be in [center-pi/2,center+pi/2]
      angle_list[which(angle_list < center-pi/2)] <- 
        angle_list[which(angle_list < center-pi/2)] + pi
      angle_list[which(angle_list > center+pi/2)] <-
        angle_list[which(angle_list > center+pi/2)] - pi
      return(angle_list)
    }
    xs <- rep(-smoothing_dist:smoothing_dist, (2*smoothing_dist+1))
    ys <- rep(-smoothing_dist:smoothing_dist, each=(2*smoothing_dist+1))
    which.dist <- which(sqrt(xs*xs+ys*ys)<=smoothing_dist)
    xs <- xs[which.dist] ; ys <- ys[which.dist]
    matrix_neu <- 0*matr
    for(i in 1:ncol(matr)) for(j in 1:nrow(matr)) {
      sel.x <- i + xs ; sel.y <- j + ys
      in.bounds <- which(sel.x>0 & sel.y>0 & sel.x<=ncol_grid & sel.y<=nrow_grid)
      sel.x <- sel.x[in.bounds] ; sel.y <- sel.y[in.bounds]
      sel.numbers <- sapply(1:length(in.bounds), function(i) {
        grid_number(sel.x[i], sel.y[i])})
      if(angle) {
        matrix_neu[i,j] <- center_angle(mean(center_angle(matr[sel.numbers],
                                      matr[grid_number(i,j)])), 0)
      } else {
        matrix_neu[i,j] <- mean(matr[sel.numbers])
      }
    }
    return(matrix_neu)
  }
  if(smoothing_dist == 0)
  {
    return(local_estimates)
  }
  else
  {
  return(cbind(as.vector(matrix_smoothing(0*X+local_estimates[,1],
                                          smoothing_dist=smoothing_dist)),
               as.vector(matrix_smoothing(0*X+local_estimates[,2],
                                          smoothing_dist=smoothing_dist)),
               as.vector(matrix_smoothing(0*X+local_estimates[,3],
                                          smoothing_dist=smoothing_dist, angle=T)))
  )
  }
}
clustering <- function(dist_matrix, method="average") {
  return(hclust(as.dist(t(dist_matrix), diag=TRUE),method=method))
}

calculate_and_save_clusters <- function (dissimilarity_matrix, save_name, k_clust=5,width=1000) {
  hc <- clustering(dissimilarity_matrix)
  png(paste0(path,"/clusterings/",save_name,".png"),
      width=width, height=1000, unit="px", res=300)
  clusters <- cutree(hc, k_clust)
  plot_cluster_map(clusters)
  dev.off()
  saveRDS(clusters, paste0(path, "/clusterings/",save_name,".rds"))
}
plot_cluster_map <- function(cluster_data, seed=sample(1:100,1), main="",
                             k_clust=max(cluster_data, na.rm=T), margin=0, print=T) {
  # as "plot_sel_map", but with a different colour for each cluster and without
  # a legend
  set.seed(seed) # to determine color ordering
  if(margin==0) {
    gr = exp.grid
  } else {
    gr <- expand.grid(y=y_ax[(margin+1):(length(y_ax)-margin)],
                      x=x_ax[(margin+1):(length(x_ax)-margin)])
  }
  # colors are randomized to prevent adjacent clusters from having
  # similar colors
  l <- levelplot(cluster_data~x*y, data=gr, col.regions = 
                    brewer.pal(k_clust, "Dark2"),
                  colorkey=F, at=1:(k_clust+1)-0.5,xlab="",ylab="",
                  par.settings =  list(layout.heights = list(
                    top.padding = 0,
                    main.key.padding = 0,
                    key.axis.padding = 0,
                    axis.xlab.padding = 0,
                    xlab.key.padding = 0,
                    key.sub.padding = 0),
                    layout.widths = list(
                      left.padding = 0,
                      key.ylab.padding = 0,
                      ylab.axis.padding = 0,
                      axis.key.padding = 0,
                      right.padding = -1)))
  if(print==T) {print(l)} else {return(l)}
}
cluster_number_threshold_method <- function(hc, threshold) {
  return(length(which(hc$height>threshold)))
}
llh_in_cluster <- function(z,df,alpha,X,Y,
                           locest,
                           max_dist=0, average=F) {
  n_grid <- length(X)
  if(length(dim(z))==3) {
    z <- sapply(1:(dim(z)[3]), function(i) {as.vector(z[,,i])})
  }
  if(nrow(z) != n_grid) stop("nrow(z) must equal the nrow(distance_matrix)")
  ilist <- rep(1,n_grid)%*%t((1:n_grid)) ; ilist <- ilist[lower.tri(ilist)]
  jlist <- (1:n_grid)%*%t(rep(1,n_grid)) ; jlist <- jlist[lower.tri(jlist)]
  Xlist <- rep(X[ilist]-X[jlist],each=dim(z)[2])
  Ylist <- rep(Y[ilist]-Y[jlist],each=dim(z)[2])
  
  zilist <- as.vector(t(z[ilist,]))
  zjlist <- as.vector(t(z[jlist,]))     
  
  if(max_dist>0) {
    sel <- which(Xlist*Xlist+Ylist*Ylist<=max_dist*max_dist)
    zilist <- zilist[sel]
    zjlist <- zjlist[sel]
    Xlist  <-  Xlist[sel]
    Ylist  <-  Ylist[sel]
  }
  
  if(length(zilist)>0)
  {
    llh <- function(par) {
      return(sum(pairwise_density_summand(zilist,zjlist,Xlist,Ylist, df, alpha,
                                          par[1],par[2],par[3])))
    }
    lh <- llh(locest)
    if(average) lh<-lh/length(zilist)
    return(lh)
  }
  else return(0)
}


calc_estimates_in_clusters <- function(clusters, df, alpha, upperbounds,
                                       clusternum=1:5) {
  # needs also: resolution, a_sim_exp_ns
  estimates_in_clusters <- matrix(-Inf, nrow=max(clusters), ncol=5)
  for(i in clusternum) {
    print(i)
    which_cl <- which(clusters==i)
    if(length(which_cl) >= 5)
    {
    estimates_in_clusters[i,1:3] <- pairwise_density_optim(
      t(sapply(which_cl, function(j) {
        a_sim_exp_ns[number_grid(j)[1],number_grid(j)[2],]})),
      df, alpha, X[which_cl], Y[which_cl] , upper_bounds = upperbounds,
      max_dist=4*((max(X)-min(X))/(resolution-1)), print=T, ensemble=3)
    estimates_in_clusters[i,4] <- length(which(clusters==i))
    estimates_in_clusters[i,5] <- llh_in_cluster(
      t(sapply(which_cl, function(j) {
        a_sim_exp_ns[number_grid(j)[1],number_grid(j)[2],]})),
      df, alpha, X[which_cl], Y[which_cl] , 
      estimates_in_clusters[i,1:3],
      max_dist=4*((max(X)-min(X))/(resolution-1)),average=T)
    }
  }
  return(estimates_in_clusters)
}

plot_llh <- function(clusters, estimates, name) {
  llh <- 0*X + sapply(clusters,
                            function(i) {estimates[i,4]})
  png(paste0(path,"estimates_",name,"_b.png"),
      width=450, height=400, unit="px", res=100)
  plot_map(llh, von=-4.5, bis=-3.5, cr=brewer.pal(10,"RdBu"))
  dev.off()
}
