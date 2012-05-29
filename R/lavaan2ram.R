## Lavaan to ram
lavaan2ram<-function(fitModel, digits=2, zero.print="0", ram.out=TRUE, fit=FALSE){
	parTable<-fitModel@ParTable
	parEst<-fitModel@Fit@est
	parSE<-fitModel@Fit@se
	fitInd<-fitMeasures(fitModel)
	ngroup<-fitModel@Data@ngroups
	if (ngroup>1){
		A<-S<-Ase<-Sse<-M<-Mse<-list()
		for(g in 1:ngroup) {
        	obsVar <- lavaan:::vnames(parTable, "ov", group=g)
        	latVar <- lavaan:::vnames(parTable, "lv", group=g)
        	varName<-c(obsVar, latVar)
			manifest<-length(obsVar)
			latent<-length(latVar)
			nrow<-length(varName)
			A[[g]]<-S[[g]]<-Ase[[g]]<-Sse[[g]]<-matrix(0,nrow,nrow,dimnames=list(varName, varName))
			M[[g]]<-Mse[[g]]<-matrix(0,nrow,1,dimnames=list(varName, 'M'))
			
			for (j in parTable$id){
				if (parTable$group[j]==g){
					if (parTable$op[j]=="~"){
						A[[g]][parTable$lhs[j], parTable$rhs[j]]<-parEst[j]
						Ase[[g]][parTable$lhs[j], parTable$rhs[j]]<-parSE[j]
					}
					if (parTable$op[j]=="=~"){
						A[[g]][parTable$rhs[j], parTable$lhs[j]]<-parEst[j]
						Ase[[g]][parTable$rhs[j], parTable$lhs[j]]<-parSE[j]
					}
					if (parTable$op[j]=="~~"){
						S[[g]][parTable$lhs[j], parTable$rhs[j]]<-parEst[j]
						Sse[[g]][parTable$lhs[j], parTable$rhs[j]]<-parSE[j]
						S[[g]][parTable$rhs[j], parTable$lhs[j]]<-parEst[j]
						Sse[[g]][parTable$rhs[j], parTable$lhs[j]]<-parSE[j]
					}
					if (parTable$op[j]=="~1"){
						M[[g]][parTable$lhs[j], 'M']<-parEst[j]
						Mse[[g]][parTable$lhs[j], 'M']<-parSE[j]
					}
				}
			}
			
		}
		
		## Print some results
		## Print the Model fit
		if (fit){
			cat("Model fit statistics and indices\n")
			print(fitInd)
		}
		if (ram.out){
			for (g in 1:ngroup){
				## Print the parameter estimates
				A.na<-A[[g]]
				A.na[A[[g]]==0]<-NA
				S.na<-S[[g]]
				S.na[S[[g]]==0]<-NA
				Ase.na<-Ase[[g]]
				Ase.na[Ase[[g]]==0]<-NA
				Sse.na<-Sse[[g]]
				Sse.na[Sse[[g]]==0]<-NA
				cat("\n--------------------\n")
	 			cat("Group ",g,"\n")
	  			cat("--------------------\n")
			
	  			cat("\n--------------------\n")
	  			cat("Parameter estimates:\n")
	  			cat("--------------------\n")
	  			cat("\nMatrix A\n\n")
	  			print(A.na, digits=digits,na.print = zero.print)
	  			cat("\nMatrix S\n\n")
	  			print(S.na,digits=digits,na.print = zero.print)
	  
	  			cat("\n----------------------------------------\n")
	  			cat("Standard errors for parameter estimates:\n")
	  			cat("----------------------------------------\n")
	  			cat("\nMatrix A\n\n")
	  			print(Ase.na,digits=digits,na.print = zero.print)
	  			cat("\nMatrix S\n\n")
	  			print(Sse.na,digits=digits,na.print = zero.print)
	  			cat("\n\n")
			}
		}
	  	lname<-NULL
  		if (latent>0) lname = varName[(manifest+1):nrow]
		invisible(list(A=A, S=S, Ase=Ase, Sse=Sse, M=M, Mse=Mse, fit=fitInd, lavaan=fitModel, nvar=nrow, manifest=manifest,latent=latent,lname=lname,varname=varName))        	
	}else{
		obsVar <- lavaan:::vnames(parTable, "ov")
        latVar <- lavaan:::vnames(parTable, "lv")
		varName<-c(obsVar, latVar)
		manifest<-length(obsVar)
		latent<-length(latVar)
	
		nrow<-length(varName)
		A<-S<-Ase<-Sse<-matrix(0,nrow,nrow,dimnames=list(varName, varName))
		M<-Mse<-matrix(0,nrow,1,dimnames=list(varName, 'M'))
		for (j in parTable$id){
			if (parTable$op[j]=="~"){
				A[parTable$lhs[j], parTable$rhs[j]]<-parEst[j]
				Ase[parTable$lhs[j], parTable$rhs[j]]<-parSE[j]
			}
			if (parTable$op[j]=="=~"){
				A[parTable$rhs[j], parTable$lhs[j]]<-parEst[j]
				Ase[parTable$rhs[j], parTable$lhs[j]]<-parSE[j]
			}
			if (parTable$op[j]=="~~"){
				S[parTable$lhs[j], parTable$rhs[j]]<-parEst[j]
				Sse[parTable$lhs[j], parTable$rhs[j]]<-parSE[j]
				S[parTable$rhs[j], parTable$lhs[j]]<-parEst[j]
				Sse[parTable$rhs[j], parTable$lhs[j]]<-parSE[j]
			}
			if (parTable$op[j]=="~1"){
				M[parTable$lhs[j], 'M']<-parEst[j]
				Mse[parTable$lhs[j], 'M']<-parSE[j]
			}
		}
	
	## Print some results
	## Print the Model fit
	if (fit){
	cat("Model fit statistics and indices\n")
	print(fitInd)
	}
	## Print the parameter estimates
	A.na<-A
	A.na[A==0]<-NA
	S.na<-S
	S.na[S==0]<-NA
	Ase.na<-Ase
	Ase.na[Ase==0]<-NA
	Sse.na<-Sse
	Sse.na[Sse==0]<-NA
	
	if (ram.out){
	  cat("\n--------------------\n")
	  cat("Parameter estimates:\n")
	  cat("--------------------\n")
	  cat("\nMatrix A\n\n")
	  print(A.na, digits=digits,na.print = zero.print)
	  cat("\nMatrix S\n\n")
	  print(S.na,digits=digits,na.print = zero.print)
	  
	  cat("\n----------------------------------------\n")
	  cat("Standard errors for parameter estimates:\n")
	  cat("----------------------------------------\n")
	  cat("\nMatrix A\n\n")
	  print(Ase.na,digits=digits,na.print = zero.print)
	  cat("\nMatrix S\n\n")
	  print(Sse.na,digits=digits,na.print = zero.print)
		cat("\n\n")
	}
  lname<-NULL
  if (latent>0) lname = varName[(manifest+1):nrow]
	invisible(list(A=A, S=S, Ase=Ase, Sse=Sse,  M=M, Mse=Mse, fit=fitInd, lavaan=fitModel, nvar=nrow, manifest=manifest,latent=latent,lname=lname,varname=varName))
	}
}

ramVF<-function(ramout, ylim, xlim, ninterval=10, scale=.1, length=.25, scatter=TRUE, n=20, alpha=.95, ...){
	ind<-which(ramout$lavaan@ParTable$label=="betax")[1]
	betax<-ramout$lavaan@Fit@est[ind]
	
	ind<-which(ramout$lavaan@ParTable$label=="gammax")[1]
	gammax<-ramout$lavaan@Fit@est[ind]
	
	ind<-which(ramout$lavaan@ParTable$label=="betay")[1]
	betay<-ramout$lavaan@Fit@est[ind]
	
	ind<-which(ramout$lavaan@ParTable$label=="gammay")[1]
	gammay<-ramout$lavaan@Fit@est[ind]
	
	ind<-which(ramout$lavaan@ParTable$label=="mxs")[1]
	mux<-ramout$lavaan@Fit@est[ind]
	
	ind<-which(ramout$lavaan@ParTable$label=="mys")[1]
	muy<-ramout$lavaan@Fit@est[ind]
	
	x<-seq(xlim[1],xlim[2], length=ninterval)
	y<-seq(ylim[1],ylim[2], length=ninterval)
	xy<-expand.grid(x,y)
	
	x<-xy[,1]
	y<-xy[,2]
	
	dx<-mux + betax*x + gammay*y
	dy<-muy + betay*y + gammax*x
	
	x1<-x+scale*dx
	y1<-y+scale*dy
	
	plot(x,y,type='n', ...)
	arrows(x,y,x1,y1,length=length,...)	
	
	if (scatter){
		alldata<-ramout$lavaan@Data@X[[1]]
		xdata<-c(alldata[1:n, ramout$info$x])
		ydata<-c(alldata[1:n, ramout$info$y])
		points(xdata, ydata, col='grey', ...)
		
		## add confidence interval
		##
		xall <- c(alldata[,ramout$info$x])
		yall <- c(alldata[,ramout$info$y])
		
		require(ellipse)		       
		       lines(ellipse(cor(xall, yall), level=alpha, scale=c(sd(xall),sd(yall)), centre=c(mean(xall),mean(yall))), lwd=1.5, col="green")
	}
	
	invisible(cbind(x,y,x1,y1))
}
