 plsalgo1dim=function(Xh,Yh)
 { 
 	X.temp=Xh
 	Y.temp=Yh
  
	n=nrow(X.temp)
	p=ncol(X.temp)
	q=ncol(Y.temp)

    tol = 1e-06

  	na.X = FALSE
    na.Y = FALSE
    is.na.X = is.na(X.temp)
    is.na.Y = is.na(Y.temp)
	if (any(is.na.X)) na.X = TRUE
    if (any(is.na.Y)) na.Y = TRUE
    n.ones = rep(1, n)
	p.ones = rep(1, p)
  	q.ones = rep(1, q)
  
   #nx = p - keepX[h]
    #    ny = q - keepY[h]
         
        #-- svd de M = t(X)*Y --#
        X.aux = X.temp		
        if (na.X) X.aux[is.na.X] = 0
         
        Y.aux = Y.temp       	
        if (na.Y) Y.aux[is.na.Y] = 0
         
        M = crossprod(X.aux, Y.aux)
        svd.M = svd(M, nu = 1, nv = 1)
        a.old = svd.M$u
        b.old = svd.M$v
         
        #-- latent variables --#
        if (na.X) {
            t = X.aux %*% a.old
            A = drop(a.old) %o% n.ones
            A[t(is.na.X)] = 0
            a.norm = crossprod(A)
            t = t / diag(a.norm)
            t = t / drop(sqrt(crossprod(t)))			
        }
        else {
            t = X.temp %*% a.old ##/ drop(crossprod(a.old))
            t = t / drop(sqrt(crossprod(t)))
        }
         
        if (na.Y) {
            u = Y.aux %*% b.old
            B = drop(b.old) %o% n.ones
            B[t(is.na.Y)] = 0
            b.norm = crossprod(B)
            u = u / diag(b.norm)
            u = u / drop(sqrt(crossprod(u)))			
        }
        else {
            u = Y.temp %*% b.old ##/ drop(crossprod(b.old))
            u = u / drop(sqrt(crossprod(u)))
        }
         
        iter = 1
         
        #-- boucle jusqu'? convergence de a et de b --#
        repeat {
            if (na.X) a = t(X.aux) %*% u
            else a = t(X.temp) %*% u
			
			if (na.Y) b = t(Y.aux) %*% t
            else b = t(Y.temp) %*% t
             
            a = a / drop(sqrt(crossprod(a)))
            b = b / drop(crossprod(t))
			 
            if (na.X) {
                t = X.aux %*% a
                A = drop(a) %o% n.ones
                A[t(is.na.X)] = 0
                a.norm = crossprod(A)
                t = t / diag(a.norm)
                t = t / drop(sqrt(crossprod(t)))			
            }
            else {
                t = X.temp %*% a ##/ drop(crossprod(a))
                t = t / drop(sqrt(crossprod(t)))
            }
             
            if (na.Y) {
                u = Y.aux %*% b
                B = drop(b) %o% n.ones
                B[t(is.na.Y)] = 0
                b.norm = crossprod(B)
                u = u / diag(b.norm)
                u = u / drop(sqrt(crossprod(u)))			
            }
            else {
                u = Y.temp %*% b ##/ drop(crossprod(b))
                u = u / drop(sqrt(crossprod(u)))
            }
           
            if (crossprod(a - a.old) < tol) break
             
            if (iter == max.iter) {
                warning(paste("Maximum number of iterations reached for the component", h),
                        call. = FALSE)
                break
            }
             
            a.old = a
            b.old = b
            iter = iter + 1
        }
         
   out=list(uload=a,vload=b,variate.X=t,variate.Y=u)

         
}#end function