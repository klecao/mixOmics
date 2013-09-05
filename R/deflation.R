deflation=function(Xh,Yh,uload,vload,mode,variate.X,variate.Y)
{
	X.temp=Xh
 	Y.temp=Yh
	t=variate.X
	u=variate.Y
	n=nrow(X.temp)
	p=ncol(X.temp)
	q=ncol(Y.temp)
	
	
	na.X = FALSE
    na.Y = FALSE
    is.na.X = is.na(X.temp)
    is.na.Y = is.na(Y.temp)
	if (any(is.na.X)) na.X = TRUE
    if (any(is.na.Y)) na.Y = TRUE
    n.ones = rep(1, n)
	p.ones = rep(1, p)
  	q.ones = rep(1, q)

       #-- deflation des matrices --#
        if (na.X) {
            X.aux = X.temp
            X.aux[is.na.X] = 0
            c = crossprod(X.aux, t)				
            T = drop(t) %o% p.ones
            T[is.na.X] = 0
            t.norm = crossprod(T)				
            c = c / diag(t.norm)
        }
        else {
            c = crossprod(X.temp, t) / drop(crossprod(t))
        }	
		
        X.temp = X.temp - t %*% t(c)   
         
        #-- mode canonique --#
        if (mode == "canonical") {
            if (na.Y) {
                Y.aux = Y.temp
                Y.aux[is.na.Y] = 0
                e = crossprod(Y.aux, u)
                U = drop(u) %o% q.ones
                U[is.na.Y] = 0
                u.norm = crossprod(U)				
                e = e / diag(u.norm)					
            }
            else {
                e = crossprod(Y.temp, u) / drop(crossprod(u))
            }
			
            Y.temp = Y.temp - u %*% t(e)
        }
         
        #-- mode regression --#
        if(mode == "regression") {
            if (TRUE) {#na.Y
                Y.aux = Y.temp
                Y.aux[is.na.Y] = 0
                d = crossprod(Y.aux, t)
                T = drop(t) %o% q.ones
                T[is.na.Y] = 0
                t.norm = crossprod(T)				
                d = d / diag(t.norm)
            }
            else {				
                d = crossprod(Y.temp, t) / drop(crossprod(t))
            }
             
            Y.temp = Y.temp - t %*% t(d)
        }

out=list(Xh=X.temp,Yh=Y.temp,variate.X=t,variate.Y=u)
}