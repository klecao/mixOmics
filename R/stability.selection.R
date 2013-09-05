stability.function.spls=function(X, Y = NULL, 
						 DA=FALSE, # discriminant analysis
                         mode = c("regression", "canonical"),
                         ncomp=1, #integer, which component are we looking at?
						  keepX = ncol(X),
						  keepY = ncol(Y),
						  near.zero.var=TRUE,
                        indice.X, #list of length ncomp, names of the variables to keep on each comp
                        indice.Y,
                         method = c("Mfold","bootstrap","LOO"),
                         num.method=10, #Mfold, number bootstrap
                         repeat.CV=1,...)
{

#-- validation des arguments --#
    if (length(dim(X)) != 2 || !is.numeric(X)) 
        stop("'X' must be a numeric matrix.")
     
    X = as.matrix(X)
    n = nrow(X)
	p = ncol(X)

	if((missing(indice.X))&(ncomp>1))         
	stop("indice.X is missing")
 

	if(DA==TRUE)
	{
		# / Testing the input Y
    if (is.null(dim(Y))) {
        Y = as.factor(Y)	
        ind.mat = unmap(as.numeric(Y))					
    	}else{stop("'Y' should be a factor or a class vector.")}	
	
	    if ((n != nrow(ind.mat))) 
        stop("unequal number of rows in 'X' and 'Y'.")
        
		Y=ind.mat
		mode="regression"
		keepY=ncol(Y) #no selection in spls-da
		if(ncomp>1)
		{indice.Y=list();for(i in 1:(ncomp-1)) indice.Y[[i]]=1:ncol(Y)}
	}else{
		Y = as.matrix(Y)
		if ((n != nrow(Y))) stop("unequal number of rows in 'X' and 'Y'.")
    	mode = match.arg(mode)
		if((missing(indice.Y))&(ncomp>1)) {indice.Y=list(ncomp-1);for(i in 1:(ncomp-1)) indice.Y[[i]]=1:ncol(Y)}
	}

    q = ncol(Y)      
    if (is.null(ncomp) || !is.numeric(ncomp) || ncomp <= 0)
        stop("invalid number of variates, 'ncomp'.")
    
    ncomp = round(ncomp)
    if(ncomp > p) {
        warning("Reset maximum number of variates 'ncomp' to ncol(X) = ", p, ".")
        ncomp = p}
	
    if (length(keepX) != 1) 
        stop("length of 'keepX' must be equal to ", 1, ".")
     
    if (length(keepY) != 1) 
        stop("length of 'keepY' must be equal to ", 1, ".")
     
    if (any(keepX > p)) 
        stop("each component of 'keepX' must be lower or equal than ", p, ".")
     
    if (any(keepY > q)) 
        stop("each component of 'keepY' must be lower or equal than ", q, ".")
        
    if(ncomp>1){
    	if (any(sapply(indice.X,length))>p)
        	stop("each component of 'indice.X' must be lower or equal than ", p, ".")
        
        if(length(indice.X)!=(ncomp-1))
        	stop("indice.X must be of length ",ncomp-1,".")	
    
    	if (any(sapply(indice.Y,length))>q)
        	stop("each component of 'indice.Y' must be lower or equal than ", q, ".")

        if(length(indice.Y)!=(ncomp-1))
        	stop("indice.Y must be of length ",ncomp-1,".")	
     	}
     
#-- initialisation des matrices --#
    X.names = dimnames(X)[[2]]
    if (is.null(X.names)) X.names = paste("X", 1:p, sep = "")
     
    if (dim(Y)[2] == 1) Y.names = "Y"
    else {
        Y.names = dimnames(Y)[[2]]
        if (is.null(Y.names)) Y.names = paste("Y", 1:q, sep = "")
    }
     
    ind.names = dimnames(X)[[1]]
    if (is.null(ind.names)) {
        ind.names = dimnames(Y)[[1]]
        rownames(X) = ind.names
    }
     	
    if (is.null(ind.names)) {
        ind.names = 1:n
        rownames(X) = rownames(Y) = ind.names
    }
     

	X=scale(X)
    Y=scale(Y)
    means.X=attr(X,"scaled:center")
    sigma.X=attr(X,"scaled:scale")
    means.Y=attr(Y,"scaled:center")
    sigma.Y=attr(Y,"scaled:scale")


	if(near.zero.var == TRUE){ 
		nzv=nearZeroVar(X)
		if(length(nzv$Position)>0)
		{
    	  	warning("Zero- or near-zero variance predictors. 
  			Reset predictors matrix to not near-zero variance predictors.
  			See $nzv for problematic predictors.")
			X2=X[,-nzv$Position]
			if(ncomp>1)
			{
				a=match(unlist(indice.X),colnames(X))
				indice.X2.unlist=match(colnames(X)[a],colnames(X2))
				indice.X2=relist(indice.X2.unlist,indice.X) 
				#peut avoir des Na si un indice de indice.X est dans nzv$Position
				for(i in 1:(ncomp-1))
				{
					a=indice.X2[[i]]
					if(sum(is.na(a))>0) indice.X2[[i]]=indice.X2[[i]][-which(is.na(a))]
					}
				}
		}else{X2=X;if(ncomp>1){indice.X2=indice.X}}
	}

#---------------- choice of the method
	TOT=matrix(0,nrow=ncol(X),ncol=1)
	if(method=="LOO"){method="Mfold";num.method=n;repeat.CV=1}
	if(method=="Mfold")
	{
		nk=n/num.method
		indice.boot=matrix(0,nrow=num.method*repeat.CV,ncol=n-nk)
		compt=0
		for(j in 1:repeat.CV)
		{
			folds = split(sample(1:n), rep(1:num.method, length = n)) 
			for(i in 1:num.method) {compt=compt+1;indice.boot[compt,]=(1:n)[-folds[[i]]]}
			}
		m=repeat.CV*num.method
		}
	
	if(method=="bootstrap")
	{
		n=nrow(X)
		indice.boot=matrix(0,nrow=num.method,ncol=n)
		for(i in 1:num.method){indice.boot[i,]=sample(1:n,n,replace=TRUE)}
		m=num.method
		}
	m=nrow(indice.boot) #should be the same as the m before; to check
print(m)
	Xh=X2
	Yh=Y
	print(indice.X2)
	print(indice.X)
	print(nzv$Position)
	print(dim(Xh))
	
#----------------- deflation of X and Y, if ncomp>1 
	if(ncomp>1) 
	{
		uload=matrix(0,nrow=p,ncol=ncomp-1)
		vload=matrix(0,nrow=q,ncol=ncomp-1)
		colnames(uload)=colnames(vload)=paste("comp",1:(ncomp-1))
		for(i in 1:(ncomp-1))
		{
	res=plsalgo1dim(Xh=Xh[,indice.X2[[i]],drop=FALSE],Yh=Yh[,indice.Y[[i]],drop=FALSE])

		signature.u=matrix(0,nrow=p,ncol=1)
		signature.v=matrix(0,nrow=q,ncol=1)
		signature.u[indice.X2[[i]]]=res$uload
		signature.v[indice.Y[[i]]]=res$vload
		
		uload[,i]=signature.u
		vload[,i]=signature.v
		
		res2=deflation(Xh,Yh,uload=signature.u,vload=signature.v,variate.X=res$variate.X,variate.Y=res$variate.Y,mode=mode)
        Yh=res2$Yh
        Xh=res2$Xh
		#print(uload)
		}
		
	}
#print(indice.boot[1,])
non.conv=0
	for(i in 1:m)
	{	
		X.train=Xh[indice.boot[i,],]
		Y.train=Yh[indice.boot[i,],]
		sp.i=try(spls(X.train,Y.train,keepX=keepX,keepY=keepY,ncomp=1,near.zero.var=FALSE,mode=mode),TRUE)
		#print(sp.i)
		if(any(class(sp.i)!="try-error")){
		indice=which(sp.i$loadings$X[,1]!=0)
		indice=names(indice)
		a=match(indice,colnames(X))
		TOT[a,1]=TOT[a,1]+1 #record the selected variables			
		}else{non.conv=non.conv+1}
	
	}
	TOT=TOT/m	

out=list(frequency=TOT,data=list(X=X,Y=Y),keepX=keepX,keepY=keepY,method=method,num.method=num.method,repeat.CV=repeat.CV,mode=mode,ncomp=ncomp,nzv=nzv,non.conv=non.conv,call=match.call())	
if(ncomp>1)
{
out$indice.X=indice.X
out$indice.Y=indice.Y
out$loadings.X=uload
out$loadings.Y=vload
}
#out=list(frequency=TOT)

if(DA) {class(out)="splsda.stab"}else{class(out)="spls.stab"}

return(out)
}