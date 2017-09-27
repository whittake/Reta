#======================================================================
# auxiliary R functions for real and simulated data
#======================================================================

#####------------------------------------------------
#' @title  xamplesDataJ
#' @name xamplesData
#'
#'
## @include wine.R
#'
#' @description
#'  Auxiliary examples  of well known small data sets.
#'
#' @examples
#'  
#'  require(reta) 
#'  xdata <- reinisfn()   # Type the function name, to get more detail
#'  head(xdata)   
#'  xdata <- mathfn()
#'  head(xdata)
#'
#' @export
#'
xamplesData <- function(){NULL}
##------------------------------------------------


##------------------------------------------------
#' @rdname  xamplesData
#' @description reinisfn
#' Makes the  Reinis data on CHD  available.\n\n
#' Converts a contingency table to a  binary data matrix.
#'
reinisfn = function(){
    require(gRbase) ;  data(reinis) ; dim(reinis)  # table
    # colnames(reinis)=c("smo","mwk","pwk","sys","pro","fam")
    dimnames(reinis)
    xdata = (arr2df(sum(reinis),reinis)=="y")+0; rownames(xdata)=NULL
    dimnames(xdata)
    dim(xdata) # 1829    6
    xdata
}
##------------------------------------------------

##------------------------------------------------
#' @rdname  xamplesData
#' @description mathfn: Makes the  mathmarks  data  available.
#' out = mathfn(); dimnames(out)  ; class(dimnames(out))
#'
mathfn = function(){
    require(gRbase) ;  data(mathmark); dim(mathmark) ;
    head(mathmark)  # 88 5
    mathmark
}
##------------------------------------------------

##------------------------------------------------
#' @rdname  xamplesData
#' @description carcassfn: Makes the  carcass  data  available.
carcassfn =  function(){
    require(gRim) ;   data(carcass) ; str(carcass)
    head(carcass)
    carcass
}
##------------------------------------------------

##------------------------------------------------
#' @rdname  xamplesData
#' @description testCatDat: Makes a trivial binary  data matrix  available.
testCatDat = function(){
#' testCatDat()
    xdf = expand.grid(c(0,1),c(0,1),c(0,1),c(0,1))
    xdf = rbind(xdf,c(0,1,0,1))
    xdf = rbind(xdf,c(0,1,0,0))
    xdf = rbind(xdf,c(0,0,1,1))
    xdf = rbind(xdf,c(0,0,0,1))
    xdf = rbind(xdf,c(0,0,0,1))
    xdf = rbind(xdf,c(0,1,0,0))
    colnames(xdf)= c('a','b','c','d')
    return(xdf)
}#end testCatDat
##------------------------------------------------


#####------------------------------------------------
#' @title  xamplesSimJ
#' @name xamplesSim
#'
#' @description
#'  Auxiliary examples  for generation of input : graph and data.
#'
#' @examples
#' 
#'  require(reta) ;  require(Rgraphviz) 
#'  plot( as(jsimGraph(p=5,edgeprob=.3),"graphNEL"))
#'  xdata <- jsimAR(n=400,p=6,alpha=.4)
#'  head(xdata)   
#'  out = exactVar(process="treeAv", p=4,alpha=-.4)
#'  names(out)
#'  jsimBinDG(10)
#'
#' @export
#'
xamplesSim <- function(){NULL}
##------------------------------------------------

##-----------------------------------------------
#' @rdname xamplesSim
#' @description jsimGraph: Select a random graph.
#' Bernoulli random edge model with option for fixing the number of edges.
#' Edges are chosen independently with the given edge probability.
#' Sampled without replacement if fixed is FALSE,
#' sampled with replacement if fixed is TRUE.
#' Returns the adjacency matrix.
#'#  plot( as(jsimGraph(p=5,edgeprob=.3),"graphNEL"))
#'#  g1 <- graph::randomEGraph(1:5, .2)
#'
jsimGraph = function(p,edgeprob=.3,fixed=FALSE){
    "
"
    edges = expand.grid(1:p,1:p)        # has row and col names
    edges = edges[edges[,1]>edges[,2],] # all undirected edges
    rownames(edges) = NULL
    nedges =  p*(p-1)/2
    if (fixed){
        ind = sample(1:nedges,size=ceiling(edgeprob*nedges),replace=FALSE)
    }else{
        ind = which(runif(nedges) <= edgeprob)
    }
    return( edges2adjmat(p,edges[ind,]) )
}
##----------------------------------------------------

##------------------------------------------------
#' @rdname xamplesSim
#' @description  jsimBinDG: Simulation of dependent binary data matrix
#'  from a Bayes net.
# 'Input:   n sample size.  
#' Output: rectangular matrix nx6 with independent rows.
#' #  jsimBinDG(20)
#'
jsimBinDG=function( n = 4000){
  x1 = as.numeric(runif(n)<.5)
  x2 = as.numeric(runif(n)<.4)  # (x1==1)*   always give a 0 
  x3 = as.numeric(runif(n)<.2) #(x2==1)*as.numeric(runif(n)<.2)+  
       # (x2==0)*as.numeric(runif(n)<.8)
  x4 = as.numeric(runif(n)<.5)
  x5 = (x4==1)*as.numeric(runif(n)<.5)+
       (x3==1)*as.numeric(runif(n)<.7)  #  now take three levels
  # x6 = (x5==1)*as.numeric(runif(n)<.3)
  x6 = (x1+x2+x3)>2 ; 
  x = cbind(x1,x2,x3,x4,x5,x6)
  return( x )  
}#end    jsimBinDG
#-----------------------------------------------

##------------------------------------------------
#' @rdname  xamplesSim
#' @description  jsimAR:   realisations of an AR(1) process.
#' Input: n integer, number of independent realistations; p
#' integer,length of each realisation; alpha real, 1-dimensional vector.
#' Output: list of data and adjacency matrix.
#'  The true adjacency matrix has all edges on the first  off-diagonal,
#' ie tri-diagonal.
#' Implements this by using a Cholesky decomposition of the
#' inverse variance.
#' ##  tri-diagonal inv var matrix, alpha is correlation
#' # jsimAR(n=200,p=12,alpha=.1)
jsimAR = function(n,p,alpha=0.1,adjmat=FALSE){
    mat = abs(outer(1:p,1:p,"-"))
    invvar = 1*(mat==0) + as.numeric(alpha)*(mat==1)
    x = solve(chol(invvar)) %*% matrix(rnorm(p*n), nrow=p)
    if (adjmat){
        return( list(xdata=t(x),adjmat.true= (mat==1)*(!(alpha==0)) ))
    }else{  return(t(x)) }
}#end jsimAR
#----------------------------------------------------

##------------------------------------------------
#' @rdname  xamplesSim
#' @description jsimARdev:   realisations of an AR(1) process.
#' Input: n integer, number of independent realistations; p
#' integer,length of each realisation; alpha real, 1-dimensional vector.
#' Output: list of data and adjacency matrix.
#'  The true adjacency matrix has all edges on the first  off-diagonal,
#' ie tri-diagonal.
#' A development version.
#' The  inverse variance matrix is tri-diagonal.
#' 
jsimARdev = function(n,p,alpha=1){
    onepp = 1+0*eye(p)
    up = solve(1*!lower.tri(onepp))
    invvar = alpha*eye(p) + up + t(up)
    x = mvrnorm(n, mu=rep(0,p), Sigma=solve(invvar))
    return(x)
}#end jsimARdev
#----------------------------------------------------

#----------------------------------------------------
  # old require(MASS)
  # n=5; p=8; alpha = .2
  # eye = function(k){diag(rep(1,k))} ;   onepp = 1+0*eye(p)
  # up = solve(1*!lower.tri(onepp)) ;   invvar = alpha*eye(p) + up + t(up)
  # invvar = cov2cor(invvar)
  # x = mvrnorm(n, mu=rep(0,p), Sigma=solve(invvar))
#----------------------------------------------------

##------------------------------------------------
#' @rdname  xamplesSim
#' @description  jsimAR2: Markov chain with a 2 state memory,  an AR(2) process.
#' Input: n integer, number of independent realistations; p
#' integer,length of each realistation; alpha real, 2-dimensional vector.
#' Output: list of data and adjacency matrix.
#'  The true adjacency matrix has all edges on the first two off-diagonals;
#'  but nowhere else.
#'
jsimAR2 = function(n,p,alpha=c(.1,.1),adjmat=FALSE){
    mat = abs(outer(1:p,1:p,"-"))
    invvar = 1*(mat==0) + as.numeric(alpha[1])*(mat==1)+
        as.numeric(alpha[2])*(mat==2)
    if (min(eigen(invvar,symmetric=TRUE,only.values=TRUE)$values)<10E-4){
        return(-1)}
    x = solve(chol(invvar)) %*% matrix(rnorm(p*n), nrow=p)
    if (adjmat){
        return( list(xdata=t(x),adjmat.true=(mat==1)+(mat==2) ))
    }else{  return(t(x)) }
}#end jsimAR2
#----------------------------------------------------

##------------------------------------------------
#' @rdname  xamplesSim
#' 
#' @description jsimMA: realisations of an MA(1) process.
#' Input: n sample size, p length of realisation.
#' The true adjacency matrix is complete.
#' ##  n = 20 ;  p = 4
jsimMA = function(n,p){
    x = matrix(rnorm(n*p),p,n)
    ones = 1+0*diag(p)
    x = solve(1*!lower.tri(ones),x)
    return(t(x))
}#end jsimMA
#----------------------------------------------------




#####------------------------------------------------
#' @rdname xamplesSim
#'
#' @description jsimGaussianGM: Simulates a Gaussian sample from a given independence graph.
#'
#' Multivariate Gaussian observations are simulated with the given
#' independence graph.  The inverse variance parameters for the edges in
#' the graph are uniformly distributed, and truncated away from 0. A
#' diagonal parameter is added, and a check is made that this matrix is
#' positive definite. The data are returned.
#'
#' adjmat = jsimGraph(p=5,edgeprob=.3)
#'  jsimGaussianGM(10,adjmat)
#'  
jsimGaussianGM = function(nss,adjmat){
    edges = adjmat2edges(adjmat)   # TODO: use /igraph?
    ## Generate precision matrix from matrix of edge pairs
    p = dim(adjmat)[1]
    invvar.true = matrix(0,nrow=p,ncol=p)  # reset
    null =  apply(edges,1,function(e){
        r = runif(1,min=.2,max=1)* sign(rbinom(1,1,prob=.5)-.5)
        invvar.true[e[1],e[2]] <<- r  # need the global assign??
        invvar.true[e[2],e[1]] <<- r
    }) # modifies
    # Augment and find inverse ie varmat 
    diag(invvar.true) = rowSums(abs(invvar.true)) + 1e-2
    invvar.true = cov2cor(invvar.true)
    if(min(eigen(invvar.true)$values)<1e-9){
        stop("Precision matrix not positive definite")
    }
    # if(! all.equal(invvar.true,t(invvar.true)) ){
    #   stop("precision matrix not symmetric")    }
    # dimnames(invvar.true) = list(n.labels,n.labels) paste("X",1:p,sep='')
    # Simulate data
    var.true = solve(invvar.true)
    return(xdata = MASS::mvrnorm(nss, mu=rep(0,p), Sigma=cov2cor(var.true)))
}#end jsimGaussianGM
#------------------------------------------------



####-------------------------------------------------
#' @rdname xamplesSim
#' @description jsimCount : binary data matrix  from a log-linear model.
#' 
#' @param size,adjmat,strength=1
#' 
#' @examples
#'  jsimCount(p,adjmat)
#'  adjmat <- jsimGraph(p=5,edgeprob=.3)
#'  jsimCount(size=10,adjmat)
#'  jsimCount(size=10,adjmat,strength=2)
#'  
#'
#'  jsimCount:  Poisson counts giving log-linear data matrix.
#' Generates realisations from a log-linear model with 
#' dependencies  determined from the graph and an overall strength index.
#' Input:
#' size, approximate sample size;
#' adjmat: adjacency matrix specifying graphical model;
#' strength: scalar, strength=0 independence, =1 mild, =2+ strong dependency.
#' Output: realisations from a log-linear model.
#' 
#' Details:
#' A is the  2^p x p design matrix that lists all cells of the p-way table
#' and carries the main effects (defaulted to 0).
#' The linear predictor is built from the two way interactions
#' composed of the strength, columns of the design matrix specified by the
#' graph and an additional random strength uniformly distributed away from 0.
#' A poisson count is generated from the linear predictor.
#' The realisations are the rows of A replicated by the counts.
#'# if want to store twi design  put B = matrix(0,nrow=k,ncol=nedge)
#'# and  assigns  B[,r] = A[,i]*A[,j] in loop
#'
jsimCount = function(size,adjmat,strength=1){ 
## initialise
    if(is.matrix(adjmat)){
        p = dim(adjmat)[1]
    }else{
        p=adjmat # if argument is scalar, it must be the dimension
        adjmat = randomEdge(p)
    }
    nedge = sum(adjmat)/2
    k = 2^p 
    A = matrix(0,nrow=k,ncol=p)
    # iterate for A 
    r = k
    for (i in 1:p){  #  node i
        r = r/2      #  block size
        A[,i] = as.numeric(gl(2,r,length=k))-1
    }
    lamA = runif(p)*(rbinom(p,size=1,prob=.5)*2-1) # main effect paras
    lam = A %*% lamA
    ## initialise and iterate to build linpred lam,
    ## do not construct the design matrix for the twi
    lamB =  runif(nedge,.2,1)*(rbinom(nedge,size=1,prob=.5)*2-1)* strength
    r = 0 
    for (i in 2:p){
        for (j in 1:(i-1)){                     # cat(c(i,j),"\n")
          if  (adjmat[i,j]){   # only edges in lower tri of adjmat
              r = r+1          # edge index
              lam = lam + (A[,i]*A[,j])*lamB[r] 
          }
        }# end j
    }# end i  lam
    mu = exp(lam)*size/sum(exp(lam))
    count = rpois(k,lambda=mu)
    # replicate x-values by count 
    out =  sapply(1:k,function(row){ A[rep(row, each=count[row]),]
                                    },simplify=FALSE)

    out <- do.call("rbind", out)
    out <- out[sample.int(nrow(out)),]
    xdf <- as.data.frame(out)
    return(xdf)
}
#------------------------------------------------




#####------------------------------------------------
#' @rdname xamplesSim
#' @description  jsimTreeAv :   realisations of a  tree averaging
#'  stuctured  process.
#' n numeric, depth of tree,
#' alpha  real coefficient.
#' 
#' jsimTreeAv(4,.3)
#'
jsimTreeAv = function(n=4,alpha=.2){
    p = 15
    A = matrix(0,nrow=p,ncol=p)
    j = 1
    half = (p-1)/2
    for (i in 1:half){
        ind=c(i,i+j,i+j+1)
        A[i,ind] = c(-1,alpha,alpha)
        j = j+1 ;#    cat(ind,"\n")
    }
    diag(A[(half+1):p,(half+1):p]) = 1
    return(t(solve(A) %*% matrix(rnorm(p*n), nrow=p)))    
}#end jsimTreeAv
#---------------------------------------------------

####------------------------------------------------
#' @rdname xamplesSim
#' @description yfsetup: 
#'
#' Simulations from a Bernoulli edge graph.
#' 
#' @param
#' paras : a list of several items.
#'
#' ranseed : set
#'
yfsetup = function(paras,  ranseed=NULL){
  #  input: paras is a dataframe dim 1 by nparas
    
  #{Random seed}
  if(!is.null(ranseed)){set.seed(ranseed)}
  #  extract parameters
  p = paras$p              # number of nodes
  npRatio = paras$npRatio  # ratio sample to dimensions 
  eDens = paras$eDens      # prop edges  ie not sparse
  n = npRatio*p
      
  #{Select edges for $Gtrue$}
  alle = expand.grid(1:p,1:p)     # has row and col names
  alle = alle[alle[,1]>alle[,2],] # all undirected edges
  nall.e = nrow(alle)
  e.sel = sample(1:nall.e,size=ceiling(eDens*nall.e),replace=FALSE)
  e.true = alle[e.sel,]

  #{Generate precision matrix}
  invcov.true = matrix(0,nrow=p,ncol=p)
  nullvalue = apply(e.true,1,function(e){
     r = runif(1,min=-1,max=1);
     invcov.true[e[1],e[2]]<<-r;
     invcov.true[e[2],e[1]]<<-r}
    ) # modifies invcov.true
  if(! all.equal(invcov.true,t(invcov.true)) ){
     stop("precision matrix not symmetric")}
  diag(invcov.true) = rowSums(abs(invcov.true))+1e-4
  invcov.true = cov2cor(invcov.true)
  if(min(eigen(invcov.true)$values)<1e-9){
     stop("Precision matrix not positive definite")}
  # A = invcov.true ; R=chol(A) ; t(R) %*% R - A ; R %*% t(R) ; solve(R)
  # chol2inv(R)
  
  #{Generate data}
  X.sim = solve(chol(invcov.true)) %*% matrix(rnorm(p*n), nrow=p)
  return(  list(invcov.true=invcov.true, X.sim=t(X.sim)))
}# end jnewsetup
#----------------------------------------------------



#####------------------------------------------------
#' @title exactVar  specifies a theoretical variance
#'
#' @description exactVar
#' Specifies the  variance from a theoretical model including:
#' a Markov chain, 
#' intra-class correlation, with and without random loadings
#' 4-cycle with equal weights, thoujg
#' all coefficient vectors (loadings) are set to alpha with a random signs
#' There is  an additional perturbation if random=TRUE.
#' 
#' @param
#' process : switching variable.
#' 
#' p : integer, number of variables.
#' 
#' alpha: numeric, scalar parameter of the inverse variance.
#'
#' random : random perturbations and sign changes.
#' @return
#' A list with   the exact variance and the true adjacency matrix.
#' @examples
#'  out = exactVar(process="treeAv", p=4,alpha=-.4)
#'  # out = exactVar(process="mc1",p=12,alpha=-.4) FAILS ??
#'  out$cor
#'  out$Amat+0
#' 
exactVar = function(process,p=4,alpha=.4,random=FALSE){
    # if (!(p>0&is.integer(p))){return("**  invalid p **")}
    # p=4;alpha=.4
    mat = abs(outer(1:p,1:p,"-"))
    load = alpha* rep(1,p)
    if (random){
        load = alpha * sign(runif(p,-1,1))
        load = load + runif(p,-1,1)/10
    }
    invvar = 0*mat
    switch(process,
           cluster = {  # inverse var of process
               invvar[2,1] = load[1]; 
               invvar[3,2] = load[2]; 
               invvar[4,2] = load[3]; 
               invvar = 1*(mat==0) + invvar + t(invvar)
               Vxx = solve(invvar)
           }, # end cluster
           chain = {  #  inverse var of MC1 process
               invvar = load*(mat==1)
               invvar[upper.tri(invvar)] = 0
               invvar = 1*(mat==0) + invvar + t(invvar)
               Vxx = solve(invvar)
           }, # end chain
           decomp = { 
               invvar[2,1] = load[1]; 
               invvar[3,2] = load[2]; 
               invvar[4,2] = load[3]; 
               invvar[4,1] = load[4]; 
                invvar = 1*(mat==0) + invvar + t(invvar)
               Vxx = solve(invvar)
           }, # end decompl
           pcycle = {  # inverse var of p-cycle process
               invvar = load*((mat==1)+(mat==p-1))
               invvar[upper.tri(invvar)] = 0
               invvar = 1*(mat==0) + invvar + t(invvar)
               Vxx = solve(invvar)
           }, # end pcycle
           bayesNetA = { B = diag(rep(1,4))
               B[1,] = c(1,0,0,0)
               B[2,] = c(-load[1],1,0,0)
               B[4,] = c(-load[2],0,0,1)
               B[3,] = c(0,-load[3],1,-load[4])
               invB=solve(B); Vxx= invB %*% t(invB)
              #  BBT =  B %*% t(B) ; Vxx=solve(BBT) wrong 
           },  # end immorA
           bayesNetB = {  B = diag(rep(1,4))
               B[2,] = c(0,1,0,0)
               B[4,] = c(0,0,0,1)
               B[1,] = c(1,-load[1],0,-load[2])
               B[3,] = c(0,-load[3],1,-load[4])
               invB=solve(B); Vxx= invB %*% t(invB)
           },  # end immorB
           bayesNetC = { B = diag(rep(1,4))
               B[1,] = c(1,0,0,0)
               B[2,] = c(-load[1],-load[2],1,-load[3])
               B[3,] = c(0,0,1,0)
               B[4,] = c(0,0,0,1) 
               invB=solve(B); Vxx= invB %*% t(invB)
           },  # end immorC
           flag = { B = diag(rep(1,4))
               B[1,] = c(1,0,0,0)
               B[2,] = c(-load[1],1,0,0)
               B[4,] = c(-load[2],0,0,1)
               B[3,] = c(0,-load[3],1,-load[4]) # not right
               invB=solve(B); Vxx= invB %*% t(invB)
              #  BBT =  B %*% t(B) ; Vxx=solve(BBT) wrong 
           },  # end flag
          
           latent = { # intra-class corr latent factor random loadings
               Vxx =  load %*% t(load)  + diag(ones(p)) 
           },  # end latent
           causal = { # so called causal model
               B = diag(rep(1,4))
               B[3,2]= load[1]; B[4,1]= load[2]; B[4,3]= load[3]; B[3,4]= load[4]
               Vxx = solve(t(B)%*%B)
           }, # end causal 
           cauCom = { # so called causal model
               B = diag(rep(1,5))
               B[3,2]= load[1]; B[4,1]= load[2];
               B[3,5]= load[3]; B[4,5]= load[4];
               Vxx = solve(t(B)%*%B)[1:4,1:4]
               # solve(Vxx)  
           }, # end  cauCom
           cauSel = { # so called causal model
               B = diag(rep(1,5))
               B[3,2]= load[1]; B[4,1]= load[2];
               B[5,3]= load[3]; B[5,4]= load[4];
               V = solve(t(B)%*%B)
               Vxx = V[1:4,1:4]-as.matrix(V[1:4,5])%*%V[5,1:4]/V[5,5]
               # as.matrix(V[1:4,5]) ; solve(Vxx)  really Vxx|y
               #  solve(Vxx)   this is a chain
           }, # end cauSel 
           treeAv = { # tree-averaging: depth p
               pvbles = 2^p -1  # p changes meaning 
               A = matrix(0,nrow=pvbles,ncol=pvbles)  # AX=Z
               j = 1
               half = (pvbles-1)/2
               for (i in 1:half){
                   ind=c(i,i+j,i+j+1)
                   value = c(-1,alpha,alpha)
                   if (random){
                       value=value+c(0,runif(1,-1,1),runif(1,-1,1))/10
                   }
                   A[i,ind] = value
                   j = j+1 ;#    cat(ind,"\n")
               }
               diag(A[(half+1):pvbles,(half+1):pvbles]) = 1
               inv = solve(A)
               Vxx = inv %*% t(inv)
           }# end treeAv
           )# end switch 
    Vxx = cov2cor(Vxx) # solve(Vxx) solve(Vxx[c(1,2,4),c(1,2,4)])
    if (all(eigen(Vxx)$values>0)){
        Atrue = (abs(solve(Vxx))>.00000001)
        diag(Atrue)=FALSE
        return(list(corr=Vxx,Amat=Atrue)) 
    }
    return('not pos def\n') 
}# end exactVar
#------------------------------------------------

#======================================================================
