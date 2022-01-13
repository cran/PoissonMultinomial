#' @title Probability Mass Function of Poisson-Multinomial Distribution
#' 
#' @description Computes the pmf of Poisson-Multinomial distribution (PMD), 
#' specified by the success probability matrix, using various methods. This function 
#' is capable of computing all probability mass points as well as
#' of pmf at certain point(s). 
#' 
#' @references 
#' Lin, Z., Wang, Y., and Hong, Y. (2022). The Poisson Multinomial Distribution and Its Applications in Voting Theory, Ecological Inference, and Machine Learning, arXiv:2201.04237. 
#'  
#'
#' @param pmat      An \eqn{\rm n \times m} success probability matrix. Here, \eqn{\rm n} is the number of independent trials, and
#'                  \eqn{\rm m} is the number of categories.
#'                  Each row of pmat describes the success probability for the corresponding
#'                  trial and it should add up to 1.
#' @param method     Character string stands for the method selected by users to 
#'                   compute the cdf. The method can only be one of 
#'                   the following three: 
#'                   \code{"DFT-CF"},
#'                   \code{"NA"},
#'                   \code{"SIM"}. 
#' @param xmat       A matrix with \eqn{\rm m} columns that specifies where the pmf is to be computed. Each row of the matrix should has the form
#'                   \eqn{\rm x = (x_{1}, \ldots, x_{m})} which is used for computing 
#'                   \eqn{\rm P(X_{1}=x_{1}, \ldots, X_{m} = x_{m})}, the values of \eqn{\rm x} should sum up to \eqn{\rm n}. It can be a vector of length \eqn{\rm m}. If \code{xmat} is \code{NULL}, the pmf at all probability mass points will be computed.
#' @param B          Number of repeats used in the simulation method. It is ignored for methods other than
#'                   the \code{"SIM"} method.
#'                   
#' @details
#' Consider \eqn{\rm n} independent trials and each trial leads to a success outcome for exactly one of the \eqn{\rm m} categories. 
#' Each category has varying success probabilities from different trials. The Poisson multinomial distribution (PMD) gives the probability 
#' of any particular combination of numbers of successes for the \eqn{\rm m} categories. 
#' The success probabilities form an \eqn{\rm n \times m} matrix, which is called the success probability matrix and denoted by \code{pmat}.  
#' For the methods we applied in \code{dpmd}, \code{"DFT-CF"} is an exact method that computes all probability mass points of the distribution,
#' using multi-dimensional FFT algorithm. When the dimension of \code{pmat} increases, the computation burden of \code{"DFT-CF"} may challenge the capability 
#' of a computer because the method automatically computes all probability mass points regardless of the input of \code{xmat}.
#' 
#' \code{"SIM"} is a simulation method that generates random samples from the distribution, and uses relative frequency to estimate the pmf. Note that the accuracy and running time will be affected by user choice of \code{B}. 
#' Usually \code{B}=1e5 or 1e6 will be accurate enough. Increasing \code{B} to larger than 1e8 will heavily increase the
#' computational burden of the computer. 
#' 
#' \code{"NA"} is an approximation method that uses a multivariate normal distribution to approximate 
#' the pmf at the points specified in \code{xmat}. This method requires an input of \code{xmat}.
#'
#' Notice if \code{xmat} is not specified then it will be set as \code{NULL}. In this case, \code{dpmd} will 
#' compute the entire pmf if the chosen method is \code{"DFT-CF"} or \code{"SIM"}. 
#' If \code{xmat} is provided, only the pmf at the points specified 
#' by \code{xmat} will be outputted.
#' 
#' @return           
#' For a given \code{xmat}, \code{dpmd} returns the pmf at points specified by \code{xmat}. 
#'
#' If \code{xmat} is \code{NULL}, all probability mass points for the distribution specified by the success probability matrix \code{pmat} will be computed, and the results are
#' stored and outputted in a multi-dimensional array, denoted by \code{res}. Note the dimension of 
#' \code{pmat} is \eqn{\rm n \times m}, thus \code{res} will be an \eqn{\rm (n+1)^{(m-1)}} array. Then 
#' the value of the pmf \eqn{\rm P(X_{1}=x_{1}, \ldots, X_{m} = x_{m})} can be extracted as \eqn{\rm res[x_{1}+1, \ldots, x_{m-1}+1]}.
#'
#' For example, for the \code{pmat} matrix in the example section, the array element \code{res[1,2,1]=0.90} gives 
#' the value of the pmf \eqn{\rm P(X_{1}=0, X_{2}=1, X_{3}=0, X_{4}=2)=0.90}. 
#'                    
#' @examples
#' pp <- matrix(c(.1, .1, .1, .7, .1, .3, .3, .3, .5, .2, .1, .2), nrow = 3, byrow = TRUE)
#' x <- c(0,0,1,2) 
#' x1 <- matrix(c(0,0,1,2,2,1,0,0),nrow=2,byrow=TRUE)
#'
#' dpmd(pmat = pp)
#' dpmd(pmat = pp, xmat = x1)
#' dpmd(pmat = pp, xmat = x)
#'
#' dpmd(pmat = pp, xmat = x, method = "NA" )
#' dpmd(pmat = pp, xmat = x1, method = "NA" )
#'
#' dpmd(pmat = pp, method = "SIM", B = 1e3)
#' dpmd(pmat = pp, xmat = x, method = "SIM", B = 1e3)
#' dpmd(pmat = pp, xmat = x1, method = "SIM", B = 1e3)
#' 
#' @export
dpmd <-function(pmat, xmat = NULL, method="DFT-CF", B=1e3)
{
  if(is.vector(xmat))
  {
    xmat=as.matrix(t(xmat))
  }
  
  chck = pmat.check(pmat,xmat)
  if(chck!=1){stop(chck)}
  # xmat should be matrix
  switch(method,
         "DFT-CF"={
           mm=ncol(pmat) # ncol of pmat
           nn=nrow(pmat) # nrow of pmat
           
           nn.vec=rep(nn+1, mm-1)
           l.vec=rep(0, mm-1)
           cn.vec=cumprod(nn.vec)
           cn.vec=c(1, cn.vec[-(mm-1)])
           cn.vec=cn.vec[length(cn.vec):1]
           cn.vec=as.integer(cn.vec) #((n+1)^(m-2),...,(n+1)^2,(n+1),1)
           
           nnt=prod(nn.vec) # (n+1)^(m-1) probability mass points
           
           #browser()
           
           res0=pmn_mdfft_arma(nnt, pmat, nn.vec, l.vec, cn.vec)
           
           #example an_array[k + 27 * (j + 12 * i)]
           #print(round(res0, 9))
           
           res=array(0, nn.vec)
           
           res.expr="res[idx[1]"
           if(mm>=3)
           {
             for(i in 2:(mm-1))
             {
               res.expr=paste0(res.expr, ", idx[", i, "]")
             }
           }
           res.expr=paste0(res.expr, "]=res0[i]")
           
           #browser()
           
           #print(nnt)
           
           for(i in 1:nnt)
           {
             idx=l.vec.compute(k=i, cn.vec=cn.vec, m=mm)
             #print(idx)
             eval(parse(text=res.expr))
           }
           
           res=round(res, 10)

           if(!is.null(xmat))
           {
            nrow.x=nrow(xmat)
            res.x=rep(NA,nrow.x)
            if(nrow.x!=1)
            {
              for(j in 1:nrow.x)
              {
                idx.x=xmat[j,1:(mm-1)]+1
                res.x.expr="res.x[j]=res[idx.x[1]"
                if(mm>=3)
                {
                  for(i in 2:(mm-1))
                  {
                    res.x.expr=paste0(res.x.expr, ", idx.x[", i, "]")
                  }
                }
                res.x.expr=paste0(res.x.expr, "]")
                eval(parse(text=res.x.expr))
              }
              res=matrix(res.x,ncol=1)
            }
            else
            {
              res.x.expr="res.x=res[idx.x[1]"
              idx.x=xmat[1:(mm-1)]+1
                if(mm>=3)
                {
                  for(i in 2:(mm-1))
                  {
                    res.x.expr=paste0(res.x.expr, ", idx.x[", i, "]")
                  }
                }
                res.x.expr=paste0(res.x.expr, "]")
                eval(parse(text=res.x.expr))
                res=res.x
            }
           }
           
           
         },
         "SIM"={
            if(is.null(xmat))
            {
              mm=ncol(pmat) # ncol of pmat
              nn=nrow(pmat) # nrow of pmat
              nn.vec=rep(nn+1, mm-1)
              l.vec=rep(0, mm-1)
              cn.vec=cumprod(nn.vec)
              cn.vec=c(1, cn.vec[-(mm-1)])
              cn.vec=cn.vec[length(cn.vec):1]
              cn.vec=as.integer(cn.vec) #((n+1)^(m-2),...,(n+1)^2,(n+1),1)
              nnt=prod(nn.vec) # (n+1)^(m-1) probability mass points
             
              res0 = pmd_simulation_allpoints(pmat, nnt, l.vec, cn.vec, B)
             
              res=array(0, nn.vec)
             
              res.expr="res[idx[1]"
              if(mm>=3)
              {
                for(i in 2:(mm-1))
                {
                  res.expr=paste0(res.expr, ", idx[", i, "]")
                }
              }
              res.expr=paste0(res.expr, "]=res0[i]")
             
              #browser()
             
              #print(nnt)
             
              for(i in 1:nnt)
              {
                idx=l.vec.compute(k=i, cn.vec=cn.vec, m=mm)
                #print(idx)
                eval(parse(text=res.expr))
              }
             
              res=round(res, 10)
            }
            else
            {
              nrow.x = nrow(xmat)
              res = rep(NA,nrow.x)
              if(nrow.x!=1)
              {
                for (i in 1:nrow.x) {
                  res[i] = pmd.by.demands(xmat[i,],pmat,B)
                }
                res=matrix(res,ncol=1)
              }
              else
              {
                res = pmd.by.demands(xmat,pmat,B)
              }
            }
         },
         "NA"=   {
           mm=ncol(pmat) # m categories
           nn=nrow(pmat) # n people
           
           if(is.null(xmat))
           {
            stop("Value of xmat is not assigned.")
           }
           
           nrow.x = nrow(xmat)
           
           for (i in 1:nrow.x) 
           {
            if(sum(xmat[i,])>nn)
            {
              stop("Sum of a row of xmat greater than n.")
            }
           }
           
           mm = mm - 1
           
           # asymptotic sigma
           n = nn
           m = mm
           if(n==1) pmat = t(as.matrix(pmat[,1:m])) else pmat = as.matrix(pmat[,1:m])
           sig = matrix(0,m,m)
           for (i in 1:n) {
             sig = sig + diag(pmat[i,],nrow = m) - pmat[i,]%*%t(pmat[i,])
           }
           
           # asymptotic mu
           mu = matrix(0, nrow = 1, ncol = ncol(pmat))
           for (i in 1:n) {
             mu = mu + pmat[i,]
           }
           mu = as.vector(mu)


           res = rep(NA,nrow.x)


           if(nrow.x!=1)
           {
            for (i in 1:nrow.x) 
            {
              x_vec = xmat[i,1:mm]
              lb = as.numeric(x_vec - 0.5)
              ub = as.numeric(x_vec+0.5)
              res0 = 0
           
              res0 = mvtnorm::pmvnorm(lower=lb,upper = ub, mean = mu, sigma = sig)
              res[i] = res0[[1]]
            }
              res=matrix(res,ncol=1)
            }
            else
            {
              x_vec = xmat[1:mm]
              lb = as.numeric(x_vec - 0.5)
              ub = as.numeric(x_vec+0.5)
              res0 = 0
              res0 = mvtnorm::pmvnorm(lower=lb,upper = ub, mean = mu, sigma = sig)
              res = res0[[1]]
            }
         },
         
  )
  
  return(res)
}



########################################################################################
#' @title Cumulative Distribution Function of Poisson-Multinomial Distribution
#' 
#' @description Computes the cdf of 
#' Poisson-Multinomial distribution that is specified by the success probability matrix, 
#' using various methods.
#'  
#' @param pmat      An \eqn{\rm n \times m} success probability matrix. Here, \eqn{\rm n} is the number of independent trials, and
#'                  \eqn{\rm m} is the number of categories.
#'                  Each row of pmat describes the success probability for the corresponding
#'                  trial and it should add up to 1.
#' @param method     Character string stands for the method selected by users to 
#'                   compute the cdf. The method can only be one of 
#'                   the following three: 
#'                   \code{"DFT-CF"},
#'                   \code{"NA"},
#'                   \code{"SIM"}. 
#' @param B          Number of repeats used in the simulation method. It is ignored for methods other than
#'                   the \code{"SIM"} method.
#' @param xmat       A matrix with \eqn{\rm m} columns.  Each row has the form \eqn{\rm x = (x_{1},\ldots,x_{m})} for computing the cdf at \eqn{\rm x},
#'                   \eqn{\rm P(X_{1} \leq x_{1},\ldots, X_{m} \leq x_{m})}. It can also be a vector with length \eqn{\rm m}.
#' 
#' @details 
#' See Details in \code{dpmd} for the definition of the PMD, the introduction of notation, and the description of the three methods (\code{"DFT-CF"}, \code{"NA"}, and \code{"SIM"}).
#' \code{ppmd} computes the cdf by adding all probability 
#' mass points within hyper-dimensional space bounded by \code{x} as in the cdf. 
#' 
#' @return 
#' The value of cdf \eqn{\rm P(X_{1} \leq x_{1},\ldots, X_{m} \leq x_{m})} at 
#' \eqn{\rm x = (x_{1},\ldots, x_{m})}.
#' 
#' @examples
#' pp <- matrix(c(.1, .1, .1, .7, .1, .3, .3, .3, .5, .2, .1, .2), nrow = 3, byrow = TRUE)
#' x <- c(3,2,1,3)
#' x1 <- matrix(c(0,0,1,2,2,1,0,0),nrow=2,byrow=TRUE)
#'
#' ppmd(pmat = pp, xmat = x)
#' ppmd(pmat = pp, xmat = x1)
#'
#' ppmd(pmat = pp, xmat = x, method = "NA")
#' ppmd(pmat = pp, xmat = x1, method = "NA")
#'
#' ppmd(pmat = pp, xmat = x, method = "SIM", B = 1e3)
#' ppmd(pmat = pp, xmat = x1, method = "SIM", B = 1e3)
#'
#' @export
ppmd = function(pmat,xmat,method="DFT-CF",B=1e3){
  
  if(is.vector(xmat))
  {
    xmat=as.matrix(t(xmat))
  }
  
  chck = pmat.check(pmat,xmat)
  if(chck!=1){stop(chck)}
  x = xmat
  nrow.x = nrow(xmat)
  
  nn = nrow(pmat)
  mm = ncol(pmat)
  #idx formed
  nn.vec=rep(nn+1, mm-1)
  l.vec=rep(0, mm-1)
  cn.vec=cumprod(nn.vec)
  cn.vec=c(1, cn.vec[-(mm-1)])
  cn.vec=cn.vec[length(cn.vec):1]
  cn.vec=as.integer(cn.vec)
  
  nnt=prod(nn.vec)
  idx = as.data.frame(matrix(0,nrow = nnt,ncol=mm))
  idx0 = as.data.frame(matrix(0,nrow = nnt,ncol=(mm-1)))
  #transfer l.vec to result vector(idx)
  for(i in 1:nnt)
  {
    idx[i,1:(mm-1)]=l.vec.compute(k=i, cn.vec=cn.vec, m=mm)
    idx0[i,] = l.vec.compute(k=i, cn.vec=cn.vec, m=mm)
    idx[i,1:(mm-1)] = idx[i,1:(mm-1)]-1
    idx[i,mm] = nn - sum(idx[i,1:(mm-1)])
    #print(idx)
  }
  
  # filter probability mass points
  conditions = c()
  expr0 = 'which('
  expr1 = ')'
  for(i in 1:mm){
    conditions[i] = paste0('idx$V',i,'<=',x[i])
  }
  cond = conditions[1]
  for(i in 1:(mm-1)){
  cond = paste0(cond,'&',conditions[i+1])
  }
  expr = paste0(expr0,cond,expr1)
  index = eval(parse(text=expr))
  points = idx[index,]
  switch(method,
         "DFT-CF" = {
           res = dpmd(pmat)
           if(nrow.x==1){
            # filter probability mass points
            conditions = c()
            expr0 = 'which('
            expr1 = ')'
            for(i in 1:mm){
              conditions[i] = paste0('idx$V',i,'<=',x[i])
            }
            cond = conditions[1]
            for(i in 1:(mm-1)){
              cond = paste0(cond,'&',conditions[i+1])
            }
            expr = paste0(expr0,cond,expr1)
            index = eval(parse(text=expr))
            points = idx[index,]


            temp.index = idx0[index,]
            prob  = 0
            res.expr="prob = prob + res[temp[1]"
            if(mm>=3)
            {
              for(i in 2:(mm-1))
              {
                res.expr=paste0(res.expr, ", temp[", i, "]")
              }
            }
            res.expr=paste0(res.expr, "]")
            
            for(i in 1:nrow(temp.index)){
              temp=as.numeric(temp.index[i,])
              eval(parse(text=res.expr))
            }
           }
           prob.mat = matrix(rep(NA,nrow.x),ncol=1)
           for(j in 1:nrow.x){
            # filter probability mass points
            conditions = c()
            expr0 = 'which('
            expr1 = ')'
            for(i in 1:mm){
              conditions[i] = paste0('idx$V',i,'<=',x[j,i])
            }
            cond = conditions[1]
            for(i in 1:(mm-1)){
              cond = paste0(cond,'&',conditions[i+1])
            }
            expr = paste0(expr0,cond,expr1)
            index = eval(parse(text=expr))
            points = idx[index,]


            temp.index = idx0[index,]
            prob  = 0
            res.expr="prob = prob + res[temp[1]"
            if(mm>=3)
            {
              for(i in 2:(mm-1))
              {
                res.expr=paste0(res.expr, ", temp[", i, "]")
              }
            }
            res.expr=paste0(res.expr, "]")
            
            for(i in 1:nrow(temp.index)){
              temp =  as.numeric(temp.index[i,])
              eval(parse(text=res.expr))
            }
            prob.mat[j] = prob
           }
           prob <- prob.mat
         },
         "SIM" = {
             if(nrow.x==1){
              # filter probability mass points
              conditions = c()
              expr0 = 'which('
              expr1 = ')'
              for(i in 1:mm){
                conditions[i] = paste0('idx$V',i,'<=',x[i])
              }
              cond = conditions[1]
              for(i in 1:(mm-1)){
                cond = paste0(cond,'&',conditions[i+1])
              }
              expr = paste0(expr0,cond,expr1)
              index = eval(parse(text=expr))
              points = idx[index,]
              points.pos = points[which(points[,mm]>=0),]
              prob = 0
              for(i in 1:nrow(points.pos)){
                prob = prob + pmd.by.demands(as.numeric(points.pos[i,]),pmat,B)
              }
             }
             #nrow.x>1
             prob.mat = matrix(rep(NA,nrow.x),ncol=1)
             for (j in 1:nrow.x) {
              # filter probability mass points
              conditions = c()
              expr0 = 'which('
              expr1 = ')'
              for(i in 1:mm){
                conditions[i] = paste0('idx$V',i,'<=',x[j,i])
              }
              cond = conditions[1]
              for(i in 1:(mm-1)){
                cond = paste0(cond,'&',conditions[i+1])
              }
              expr = paste0(expr0,cond,expr1)
              index = eval(parse(text=expr))
              points = idx[index,]
              points.pos = points[which(points[,mm]>=0),]
              prob = 0
              for(i in 1:nrow(points.pos)){
                prob = prob + pmd.by.demands(as.numeric(points.pos[i,]),pmat,B)
              }
              prob.mat[j]=prob
             }
             prob=prob.mat
         },
         "NA" = {
             if(nrow.x==1){
              # filter probability mass points
              conditions = c()
              expr0 = 'which('
              expr1 = ')'
              for(i in 1:mm){
                conditions[i] = paste0('idx$V',i,'<=',x[i])
              }
              cond = conditions[1]
              for(i in 1:(mm-1)){
                cond = paste0(cond,'&',conditions[i+1])
              }
              expr = paste0(expr0,cond,expr1)
              index = eval(parse(text=expr))
              points = idx[index,]

              prob = 0
              points.pos = points[which(points[,mm]>=0),]
              if(nrow(points.pos)!=0){
                for(i in 1:nrow(points.pos)){
                  prob = prob + dpmd(pmat, xmat = as.matrix(points.pos[i,]), method="NA")
                }
              }
             }
             #nrow.x>1
             prob.mat = matrix(rep(NA,nrow.x),ncol=1)
             for (j in 1:nrow.x) {
              # filter probability mass points
              conditions = c()
              expr0 = 'which('
              expr1 = ')'
              for(i in 1:mm){
                conditions[i] = paste0('idx$V',i,'<=',x[j,i])
              }
              cond = conditions[1]
              for(i in 1:(mm-1)){
                cond = paste0(cond,'&',conditions[i+1])
              }
              expr = paste0(expr0,cond,expr1)
              index = eval(parse(text=expr))
              points = idx[index,]

              prob = 0
              points.pos = points[which(points[,mm]>=0),]
              if(nrow(points.pos)!=0){
                for(i in 1:nrow(points.pos)){
                  prob = prob + dpmd(pmat, xmat = as.matrix(points.pos[i,]), method="NA")
                }
              }
              prob.mat[j]=prob
             }
             prob=prob.mat

         })
  return(prob)
}
########################################################################################
#' @title Poisson-Multinomial Distribution Random Number Generator
#' @description Generates random samples from the PMD specified by the success probability matrix.
#'  
#' @param pmat      An \eqn{\rm n \times m} success probability matrix, where \eqn{\rm n} is the number of independent trials and \eqn{\rm m} is the number of categories.
#'                  Each row of pmat contains the success probabilities for the corresponding
#'                  trial, and each row adds up to 1.
#' @param s         The number of samples to be generated.
#' 
#' @return 
#' An \eqn{s \times m} matrix of samples, each row stands for one sample from the PMD with success probability matrix \code{pmat}.
#' 
#' @examples 
#' pp <- matrix(c(.1, .1, .1, .7, .1, .3, .3, .3, .5, .2, .1, .2), nrow = 3, byrow = TRUE)
#'  
#' rpmd(pmat = pp, s = 5)
#' 
#' @export
rpmd = function(pmat, s=1){
  chck = pmat.check(pmat)
  if(chck!=1){stop(chck)}
  mm = ncol(pmat)
  rnd = matrix(NA,nrow = s,ncol = mm)
  for(i in 1:s){
    rnd[i,] = t(rpmd_arma(pmat))
  }
  return(rnd)
}



