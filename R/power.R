## power for univariate and bivariate latent change score model
## Johnny Zhang 
## Created on Sep 26, 2016

powerLCS<-function(N=100, T=5, R=1000,
	betay=0, my0=0, mys=0, varey=1, vary0=1, varys=1, vary0ys=0, alpha=0.05, ...){
	#if (sum(N < 2*T)>0) stop("The sample size has to be at least 2 times of the number of occasions")
	
	pop.model <- function(T){
	## latent y
	## Intercept
		model<-"y0 =~ 1*y1\n"
		
		## path from y(t-1) to y(t) with path 1
		for (i in 2:T){
			model<-paste(model, "y",i,"~1*y",(i-1),"\n", sep="")
		}
		
		## loading from dy(t) to y(t) with path 1
		for (i in 2:T){
			model<-paste(model, "dy",i,"=~1*y",i,"\n", sep="")
		}
		
		## path from y(t) to dy(t+1) with path betay
		for (i in 2:T){
			model<-paste(model, "dy",i,"~", betay, "*y", (i-1), "\n", sep="")
		}
		
		## latent slope ys factor model
		for (i in 2:T){
			model<-paste(model, "ys=~1*dy", i, "\n", sep="")
		}
		
		## variance for dy constraints to 0
		for (i in 2:T){
			model<-paste(model, "dy",i,"~~0*dy",i,"\n", sep="")
		}
		
		## variance for y constraints to 0
		for (i in 1:T){
			model<-paste(model, "y",i,"~~0*y",i,"\n", sep="")
		}
		
		## variance and covariance for intercept and slope
		model<-paste(model, "ys~~", vary0ys, "*y0\n", sep="")
		model<-paste(model, "y0~~", vary0, "*y0\n", sep="")
		model<-paste(model, "ys~~", varys, "*ys\n", sep="")
		
		model<-paste(model, "ys~", mys, "*1\n", sep="")	
		model<-paste(model, "y0~", my0, "*1\n", sep="")

		## constrain means of y and dy to be zero
		for (i in 1:T){
			model<-paste(model, "y",i,"~0*1\n", sep="")
		}
		for (i in 2:T){
			model<-paste(model, "dy",i,"~0*1\n", sep="")
		}
		
		
		## for observed data part
		## y(t) to Y(t)
		for (i in 1:T){
			model<-paste(model, "y",i,"=~1*", "Y",i, "\n", sep="")				
		}
		
		## set means of Y to be zero
		for (i in 1:T){
			model<-paste(model, "Y",i, "~0*1\n", sep="")		
		}
		
		## set the variance for Y
		for (i in 1:T){
			model<-paste(model, "Y",i, "~~", varey, "*", "Y",i, "\n", sep="")		
		}
		model
	}
	
	fit.model <- function(T){
	## latent y
	## Intercept
		model<-"y0 =~ 1*y1\n"
		
		## path from y(t-1) to y(t) with path 1
		for (i in 2:T){
			model<-paste(model, "y",i,"~1*y",(i-1),"\n", sep="")
		}
		
		## loading from dy(t) to y(t) with path 1
		for (i in 2:T){
			model<-paste(model, "dy",i,"=~1*y",i,"\n", sep="")
		}
		
		## path from y(t) to dy(t+1) with path betay
		for (i in 2:T){
			model<-paste(model, "dy",i,"~start(", betay, ")*y", (i-1)," + betay*y", (i-1), "\n", sep="")
		}
		
		## latent slope ys factor model
		for (i in 2:T){
			model<-paste(model, "ys=~1*dy", i, "\n", sep="")
		}
		
		## variance for dy constraints to 0
		for (i in 2:T){
			model<-paste(model, "dy",i,"~~0*dy",i,"\n", sep="")
		}
		
		## variance for y constraints to 0
		for (i in 1:T){
			model<-paste(model, "y",i,"~~0*y",i,"\n", sep="")
		}
		
		## variance and covariance for intercept and slope
		model<-paste(model, "ys~~start(", vary0ys, ")*y0 + vary0ys*y0\n", sep="")
		model<-paste(model, "y0~~start(", vary0, ")*y0 + vary0*y0\n", sep="")
		model<-paste(model, "ys~~start(", varys, ")*ys + varys*ys\n", sep="")
		
		model<-paste(model, "ys~start(", mys, ")*1 + label('mys')*1\n", sep="")	
		model<-paste(model, "y0~start(", my0, ")*1 + label('my0')*1\n", sep="")

		## constrain means of y and dy to be zero
		for (i in 1:T){
			model<-paste(model, "y",i,"~0*1\n", sep="")
		}
		for (i in 2:T){
			model<-paste(model, "dy",i,"~0*1\n", sep="")
		}
		
		
		## for observed data part
		## y(t) to Y(t)
		for (i in 1:T){
			model<-paste(model, "y",i,"=~1*", "Y",i, "\n", sep="")				
		}
		
		## set means of Y to be zero
		for (i in 1:T){
			model<-paste(model, "Y",i, "~0*1\n", sep="")		
		}
		
		## set the variance for Y
		for (i in 1:T){
			model<-paste(model, "Y",i, "~~start(", varey, ")*", "Y", i, " + varey*Y",i, "\n", sep="")		
		}
		model
	}
	
	sem.est <- function(model, data){
		temp.res <- sem(model=model, data=data)
		label <- temp.res@ParTable$label
		c(temp.res@ParTable$est[label!=""], temp.res@ParTable$se[label!=""])
	}
	## do it once for a given N and T
	fit.once <- function(N, T){
		## generate data
		pop.model.T <- pop.model(T)
		pop.model.T.res <- sem(pop.model.T, do.fit=FALSE)
		pop.model.T.cov <- inspect(pop.model.T.res, "cov.ov")
		pop.model.T.mean <- inspect(pop.model.T.res, "mean.ov")
		ynames <- row.names(pop.model.T.cov)
		gen.data <- lapply(1:R, mvrnorm, n=N, mu=pop.model.T.mean, Sigma=pop.model.T.cov)
		
		## conduct the analysis
		fit.model.T <- fit.model(T)
		fit.res <- lapply(gen.data, sem.est, model=fit.model.T)
		
		## run once to get the model information
		model.info.res <- sem(fit.model.T, gen.data[[1]])
		label <- model.info.res@ParTable$label
		label <- label[label!=""]
		label.unique <- !duplicated(label)
		label <- label[label.unique]
		npar <- length(label)
		## get the parameter estimates, sd, se, power, CI of power
		all.res <- do.call(rbind, fit.res)
		all.res <- all.res[, c(label.unique, label.unique)]
		mc.est <- colMeans(all.res[, 1:npar])
		mc.se <- apply(all.res[, (npar+1):(2*npar)], 2, mean)
		mc.sd <- apply(all.res[, 1:npar], 2, sd)
		
		mc.z.score <- all.res[, 1:npar]/all.res[, (npar+1):(2*npar)]
		mc.z.score.check <- abs(mc.z.score) >= qnorm(1-alpha/2)
		
		mc.power <- colMeans(mc.z.score.check)
		
		pop.par <- unlist(lapply(label, function(x){eval(parse(text=x))}))
		mc.output <- cbind(pop.par, mc.est, mc.sd, mc.se, mc.power, N, T)
		row.names(mc.output) <- label
		label.sort <- sort(label)
		mc.output[label.sort, ]
	}
	
	if (length(N)>1 | length(T)>1){
		all.output <- list()
		for (i in N){
			for (j in T){
				all.output [[paste('N',i,'-T',j, sep="")]]<- fit.once(i,j)
			}
		}
	}else{
		all.output <- fit.once(N,T)
	}
	class(all.output) <- "lcs.power"
	all.output 	
}



plot.lcs.power <- function(x, parameter, ...){
	## x is the output from power analysis
	power.mat <- do.call('rbind', x)
	power.par <- power.mat[rownames(power.mat)==parameter, ]
	unique.N <- unique(power.par[ ,6])
	unique.T <- unique(power.par[ ,7])
	
	if (length(unique.N)==1 & length(unique.T)==1) stop("Multiple N or T is needed for power plot.")
	
	if (length(unique.N)==1){
		## plot the power along T
		plot(power.par[, 7], power.par[, 5], type='l', xlab='Number of Occasions', ylab='Power', ylim=c(0,1))
		points(power.par[, 7], power.par[, 5])
	}
	
	if (length(unique.T)==1){
		plot(power.par[, 6], power.par[, 5], type='l', xlab='Sample size', ylab='Power', ylim=c(0,1))
		points(power.par[, 6], power.par[, 5])
	}
	
	if (length(unique.N)>1 & length(unique.T)>1){
		for (N in unique.N){
			## plot power with time for a given sample size
			temp.power <- power.par[power.par[, 6]==N, ]
			plot(temp.power[, 7], temp.power[, 5], type='l', xlab='Number of Occasions', ylab='Power', ylim=c(0,1))
			points(temp.power[, 7], temp.power[, 5])
			cat ("Press [enter] to continue")
			line <- readline()
		}
		
		for (T in unique.T){
			## plot power with time for a given sample size
			temp.power <- power.par[power.par[, 7]==T, ]
			plot(temp.power[, 6], temp.power[, 5], type='l', xlab='Sample size', ylab='Power', ylim=c(0,1))
			points(temp.power[, 6], temp.power[, 5])
			cat ("Press [enter] to continue")
			line <- readline()
		}
	}
}

powerBLCS<-function(N=100, T=5, R=1000,
	betay=0, my0=0, mys=0, varey=1, vary0=1, varys=1, vary0ys=0, alpha=0.05,
	betax=0, mx0=0, mxs=0, varex=1, varx0=1, varxs=1, varx0xs=0, varx0y0=0,  
	varx0ys=0, vary0xs=0, varxsys=0, gammax=0, gammay=0, ...){
	
	pop.model <- function(T){
		## for y
		## latent y
		## Intercept
		model<-"y0 =~ 1*y1\n"
		
		## path from y(t-1) to y(t) with path 1
		for (i in 2:T){
			model<-paste(model, "y",i,"~1*y",(i-1),"\n", sep="")
		}
		
		## loading from dy(t) to y(t) with path 1
		for (i in 2:T){
			model<-paste(model, "dy",i,"=~1*y",i,"\n", sep="")
		}
		
		## path from y(t) to dy(t+1) with path betay
		for (i in 2:T){
			model<-paste(model, "dy",i,"~", betay, "*y", (i-1), "\n", sep="")
		}
		
		## latent slope ys factor model
		for (i in 2:T){
			model<-paste(model, "ys=~1*dy", i, "\n", sep="")
		}
		
		## variance for dy constraints to 0
		for (i in 2:T){
			model<-paste(model, "dy",i,"~~0*dy",i,"\n", sep="")
		}
		
		## variance for y constraints to 0
		for (i in 1:T){
			model<-paste(model, "y",i,"~~0*y",i,"\n", sep="")
		}
		
		## variance and covariance for intercept and slope
		model<-paste(model, "ys~~", vary0ys, "*y0\n", sep="")
		model<-paste(model, "y0~~", vary0, "*y0\n", sep="")
		model<-paste(model, "ys~~", varys, "*ys\n", sep="")
		
		model<-paste(model, "ys~", mys, "*1\n", sep="")	
		model<-paste(model, "y0~", my0, "*1\n", sep="")

		## constrain means of y and dy to be zero
		for (i in 1:T){
			model<-paste(model, "y",i,"~0*1\n", sep="")
		}
		for (i in 2:T){
			model<-paste(model, "dy",i,"~0*1\n", sep="")
		}
		
		
		## for observed data part
		## y(t) to Y(t)
		for (i in 1:T){
			model<-paste(model, "y",i,"=~1*", "Y",i, "\n", sep="")				
		}
		
		## set means of Y to be zero
		for (i in 1:T){
			model<-paste(model, "Y",i, "~0*1\n", sep="")		
		}
		
		## set the variance for Y
		for (i in 1:T){
			model<-paste(model, "Y",i, "~~", varey, "*", "Y",i, "\n", sep="")		
		}
		
		
		## for x
		## latent x
		## Intercept
		model<-paste(model, "x0 =~ 1*x1\n")
		
		## path from x(t-1) to x(t) with path 1
		for (i in 2:T){
			model<-paste(model, "x",i,"~1*x",(i-1),"\n", sep="")
		}
		
		## loading from dx(t) to x(t) with path 1
		for (i in 2:T){
			model<-paste(model, "dx",i,"=~1*x",i,"\n", sep="")
		}
		
		## path from x(t) to dx(t+1) with path betax
		for (i in 2:T){
			model<-paste(model, "dx",i,"~", betax, "*x", (i-1), "\n", sep="")
		}
		
		## latent slope xs factor model
		for (i in 2:T){
			model<-paste(model, "xs=~1*dx", i, "\n", sep="")
		}
		
		## variance for dx constraints to 0
		for (i in 2:T){
			model<-paste(model, "dx",i,"~~0*dx",i,"\n", sep="")
		}
		
		## variance for x constraints to 0
		for (i in 1:T){
			model<-paste(model, "x",i,"~~0*x",i,"\n", sep="")
		}
		
		## variance and covariance for intercept and slope
		model<-paste(model, "xs~~", varx0xs, "*x0\n", sep="")
		model<-paste(model, "x0~~", varx0, "*x0\n", sep="")
		model<-paste(model, "xs~~", varxs, "*xs\n", sep="")
		
		model<-paste(model, "xs~", mxs, "*1\n", sep="")	
		model<-paste(model, "x0~", mx0, "*1\n", sep="")

		## constrain means of x and dx to be zero
		for (i in 1:T){
			model<-paste(model, "x",i,"~0*1\n", sep="")
		}
		for (i in 2:T){
			model<-paste(model, "dx",i,"~0*1\n", sep="")
		}
		
		
		## for observed data part
		## x(t) to X(t)
		for (i in 1:T){
			model<-paste(model, "x",i,"=~1*", "X",i, "\n", sep="")				
		}
		
		## set means of X to be zero
		for (i in 1:T){
			model<-paste(model, "X",i, "~0*1\n", sep="")		
		}
		
		## set the variance for X
		for (i in 1:T){
			model<-paste(model, "X",i, "~~", varex, "*", "X",i, "\n", sep="")		
		}
	
	
	
		## coupling effects
		for (i in 2:T){
			model<-paste(model, "dy",i,"~", gammax, "*x",i-1, "\n", sep="")
		}
	
		for (i in 2:T){
			model<-paste(model, "dx",i,"~", gammay, "*y",i-1, "\n", sep="")
		}
	
		model<-paste(model, "x0~~", varx0y0, "*y0\n", sep="")	
		model<-paste(model, "x0~~", varx0ys, "*ys\n", sep="")	
		model<-paste(model, "y0~~", vary0xs, "*xs\n", sep="")
		model<-paste(model, "xs~~", varxsys, "*ys\n", sep="")
		model
	}
	
	fit.model <- function(T){
		## latent y
		## Intercept
		model<-"y0 =~ 1*y1\n"
		
		## path from y(t-1) to y(t) with path 1
		for (i in 2:T){
			model<-paste(model, "y",i,"~1*y",(i-1),"\n", sep="")
		}
		
		## loading from dy(t) to y(t) with path 1
		for (i in 2:T){
			model<-paste(model, "dy",i,"=~1*y",i,"\n", sep="")
		}
		
		## path from y(t) to dy(t+1) with path betay
		for (i in 2:T){
			model<-paste(model, "dy",i,"~start(", betay, ")*y", (i-1)," + betay*y", (i-1), "\n", sep="")
		}
		
		## latent slope ys factor model
		for (i in 2:T){
			model<-paste(model, "ys=~1*dy", i, "\n", sep="")
		}
		
		## variance for dy constraints to 0
		for (i in 2:T){
			model<-paste(model, "dy",i,"~~0*dy",i,"\n", sep="")
		}
		
		## variance for y constraints to 0
		for (i in 1:T){
			model<-paste(model, "y",i,"~~0*y",i,"\n", sep="")
		}
		
		## variance and covariance for intercept and slope
		model<-paste(model, "ys~~start(", vary0ys, ")*y0 + vary0ys*y0\n", sep="")
		model<-paste(model, "y0~~start(", vary0, ")*y0 + vary0*y0\n", sep="")
		model<-paste(model, "ys~~start(", varys, ")*ys + varys*ys\n", sep="")
		
		model<-paste(model, "ys~start(", mys, ")*1 + label('mys')*1\n", sep="")	
		model<-paste(model, "y0~start(", my0, ")*1 + label('my0')*1\n", sep="")

		## constrain means of y and dy to be zero
		for (i in 1:T){
			model<-paste(model, "y",i,"~0*1\n", sep="")
		}
		for (i in 2:T){
			model<-paste(model, "dy",i,"~0*1\n", sep="")
		}
		
		
		## for observed data part
		## y(t) to Y(t)
		for (i in 1:T){
			model<-paste(model, "y",i,"=~1*", "Y",i, "\n", sep="")				
		}
		
		## set means of Y to be zero
		for (i in 1:T){
			model<-paste(model, "Y",i, "~0*1\n", sep="")		
		}
		
		## set the variance for Y
		for (i in 1:T){
			model<-paste(model, "Y",i, "~~start(", varey, ")*", "Y", i, " + varey*Y",i, "\n", sep="")		
		}
		
		## latent x
		## Intercept
		model<-paste(model, "x0 =~ 1*x1\n")
		
		## path from x(t-1) to x(t) with path 1
		for (i in 2:T){
			model<-paste(model, "x",i,"~1*x",(i-1),"\n", sep="")
		}
		
		## loading from dx(t) to x(t) with path 1
		for (i in 2:T){
			model<-paste(model, "dx",i,"=~1*x",i,"\n", sep="")
		}
		
		## path from x(t) to dx(t+1) with path betax
		for (i in 2:T){
			model<-paste(model, "dx",i,"~start(", betax, ")*x", (i-1)," + betax*x", (i-1), "\n", sep="")
		}
		
		## latent slope xs factor model
		for (i in 2:T){
			model<-paste(model, "xs=~1*dx", i, "\n", sep="")
		}
		
		## variance for dx constraints to 0
		for (i in 2:T){
			model<-paste(model, "dx",i,"~~0*dx",i,"\n", sep="")
		}
		
		## variance for x constraints to 0
		for (i in 1:T){
			model<-paste(model, "x",i,"~~0*x",i,"\n", sep="")
		}
		
		## variance and covariance for intercept and slope
		model<-paste(model, "xs~~start(", varx0xs, ")*x0 + varx0xs*x0\n", sep="")
		model<-paste(model, "x0~~start(", varx0, ")*x0 + varx0*x0\n", sep="")
		model<-paste(model, "xs~~start(", varxs, ")*xs + varxs*xs\n", sep="")
		
		model<-paste(model, "xs~start(", mxs, ")*1 + label('mxs')*1\n", sep="")	
		model<-paste(model, "x0~start(", mx0, ")*1 + label('mx0')*1\n", sep="")

		## constrain means of x and dx to be zero
		for (i in 1:T){
			model<-paste(model, "x",i,"~0*1\n", sep="")
		}
		for (i in 2:T){
			model<-paste(model, "dx",i,"~0*1\n", sep="")
		}
		
		
		## for observed data part
		## x(t) to X(t)
		for (i in 1:T){
			model<-paste(model, "x",i,"=~1*", "X",i, "\n", sep="")				
		}
		
		## set means of X to be zero
		for (i in 1:T){
			model<-paste(model, "X",i, "~0*1\n", sep="")		
		}
		
		## set the variance for X
		for (i in 1:T){
			model<-paste(model, "X",i, "~~start(", varex, ")*", "X", i, " + varex*X",i, "\n", sep="")		
		}
		
		## coupling effects
		for (i in 2:T){
			model<-paste(model, "dy",i,"~start(", gammax, ")*x", i-1, " + gammax*x", i-1, "\n", sep="")
		}
	
		for (i in 2:T){
			model<-paste(model, "dx",i,"~start(", gammay, ")*y", i-1, " + gammay*y", i-1, "\n", sep="")
		}
		
		model<-paste(model, "x0~~start(", varx0y0, ")*y0 + varx0y0*y0\n", sep="")	
		model<-paste(model, "x0~~start(", varx0ys, ")*ys + varx0ys*ys\n", sep="")	
		model<-paste(model, "y0~~start(", vary0xs, ")*xs + vary0xs*xs\n", sep="")
		model<-paste(model, "xs~~start(", varxsys, ")*ys + varxsys*ys\n", sep="")
		
		model
	}
	
	sem.est <- function(model, data){
		temp.res <- sem(model=model, data=data)
		label <- temp.res@ParTable$label
		c(temp.res@ParTable$est[label!=""], temp.res@ParTable$se[label!=""])
	}
	## do it once for a given N and T
	fit.once <- function(N, T){
		## generate data
		pop.model.T <- pop.model(T)
		pop.model.T.res <- sem(pop.model.T, do.fit=FALSE)
		pop.model.T.cov <- inspect(pop.model.T.res, "cov.ov")
		pop.model.T.mean <- inspect(pop.model.T.res, "mean.ov")
		ynames <- row.names(pop.model.T.cov)
		gen.data <- lapply(1:R, mvrnorm, n=N, mu=pop.model.T.mean, Sigma=pop.model.T.cov)
		
		## conduct the analysis
		fit.model.T <- fit.model(T)
		fit.res <- lapply(gen.data, sem.est, model=fit.model.T)
		
		## run once to get the model information
		model.info.res <- sem(fit.model.T, gen.data[[1]])
		label <- model.info.res@ParTable$label
		label <- label[label!=""]
		label.unique <- !duplicated(label)
		label <- label[label.unique]
		npar <- length(label)
		## get the parameter estimates, sd, se, power, CI of power
		all.res <- do.call(rbind, fit.res)
		all.res <- all.res[, c(label.unique, label.unique)]
		mc.est <- colMeans(all.res[, 1:npar])
		mc.se <- apply(all.res[, (npar+1):(2*npar)], 2, mean)
		mc.sd <- apply(all.res[, 1:npar], 2, sd)
		
		mc.z.score <- all.res[, 1:npar]/all.res[, (npar+1):(2*npar)]
		mc.z.score.check <- abs(mc.z.score) >= qnorm(1-alpha/2)
		
		mc.power <- colMeans(mc.z.score.check)
		
		pop.par <- unlist(lapply(label, function(x){eval(parse(text=x))}))
		mc.output <- cbind(pop.par, mc.est, mc.sd, mc.se, mc.power, N, T)
		row.names(mc.output) <- label
		label.sort <- sort(label)
		mc.output[label.sort, ]
	}
	
	if (length(N)>1 | length(T)>1){
		all.output <- list()
		for (i in N){
			for (j in T){
				all.output [[paste('N',i,'-T',j, sep="")]]<- fit.once(i,j)
			}
		}
	}else{
		all.output <- fit.once(N,T)
	}
	class(all.output) <- "lcs.power"
	all.output 	
}



