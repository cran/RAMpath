## Lavaan to ram
lavaan2ram<-function(fitModel, digits=2, zero.print="0", ram.out=TRUE, fit=FALSE){
	parTable<-fitModel@ParTable
	parEst<-fitModel@Fit@est
	parSE<-fitModel@Fit@se
	fitInd<-fitMeasures(fitModel)
	
	varALL<-unique(c(parTable$rhs, parTable$lhs))
	varALL<-varALL[varALL!=""]
	latVar<-unique(parTable$lhs[parTable$op=="=~"])
	
	obsVar<-NULL
	for (i in 1:length(varALL)){
		if (!(varALL[i] %in% latVar)){ 
			obsVar<-c(obsVar, varALL[i])
		}
	}
	varName<-c(obsVar, latVar)
	manifest<-length(obsVar)
	latent<-length(latVar)
	
	nrow<-length(varName)
	A<-S<-Ase<-Sse<-matrix(0,nrow,nrow,dimnames=list(varName, varName))
	
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
	invisible(return(list(A=A, S=S, Ase=Ase, Sse=Sse, fit=fitInd, lavaan=fitModel, nvar=nrow, manifest=manifest,latent=latent,lname=lname,varname=varName)))
}

ramVF<-function(ramout, ylim, xlim, ninterval=10, scale=.1, length=.25, ...){
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
	
	invisible(cbind(x,y,x1,y1))
}
