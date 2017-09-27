###############################################################################
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.

# Author: Joe Whittaker 
###############################################################################
#======================================================================
# auxiliary R functions for triples
#======================================================================

#----------------------------------------------------------------------
## this comment is included
# this is not
#----------------------------------------------------------------------

#####------------------------------------------------
#' @title triples: the package
#'
#' 
#' @description
#' 
#' Triples is an R package that implements graphical model search
#' based on evaluating entropy interactions of the empirical or
#' theoretical entropy function, and in particular, 3-way
#' interactions.
#'
#' The package is object orientated and written using R5 reference
#' classes. See John Chambers ..
#'
#' Examples are given in the documentation for each reference class.
#' A good place to start is \code{\link{wine}}.
#' 
#' @usage    library(reta)
#' @title   GM search using the triples package
#' @description
#' The function \code{mytriples} 
#' @param
#' xdata : nxp data matrix.
#' @param
#' alpha : numeric, common significance level for all MI tests.
#' @return
#' Adjacency matrix determined by the search procedure of \code{triples}.
#' @examples
#'  library("reta")
#'  source("retaPkg/R/xtrafns.R")
#'  set.seed(1234)
#'  xdata <- jsimAR(n=100,p=12,alpha=.4)
#'  myadjmat <- myTriples(xdata)
#'  myadjmat + 0 
#'
#'  set.seed(1237)   ## yangs example
#'  X.data <- matrix(rnorm(3*40), ncol=3)
#'  myTriples(X.data)
#'  #   alpha micut: 0.05 70.93829
#'
#'  xdata1 <- jsimAR(n=100,p=12,alpha=.4)
#'  xdata2 <- jsimAR(n=100,p=12,alpha=.4)
#'  myCombine(xdata1,xdata2)+0
#'
myTriples = function(xdata, alpha=0.05){ 
    obj <- Updater$new(xdata)   
    obj$basicTriples(alpha)
    cat("alpha micut:",alpha ,obj$micut,"\n")       
    return(obj$adjmat)
}
#----------------------------------------------------

#####------------------------------------------------
#' @title   myCombine
#' @examples
#' require(reta)
#' set.seed(1234)
#'  n1=200 ; n2=100  ; p=6
#'   out1 <- jsimAR(n1,p,alpha=.4,adjmat=TRUE)
#'   xdata1 <- out1$xdata
#'   xdata2 <- jsimAR(n2,p,alpha=.4)
#' myCombine(xdata1,xdata2) # default type is source
#' myCombine(xdata1,xdata2,pi=.4,type="source")
#' myCombine(xdata1,xdata2,pi=.4,type="shrink")
#' myCombine(xdata1,xdata2,pi=2,type="penalty")
#' 
myCombine = function(xdata1,xdata2,   #top copy
                     pi= nrow(xdata2)/(nrow(xdata1)+nrow(xdata2)),
                     type = c("source","shrink","penalty")
                     ){
    "
myCombine(xdata1,xdata2) # default type is source
  2sources  (1-pi) d1 + pi d2 ; micut set empirically
  shrinkage (1-pi) d1 + pi d2 ; micut set as single source
  penalise         d1 + pi d2 ; micut set as single source
"
    type= match.arg(type)
    obj  <- Combine$new(xdata1,xdata2,type)
    obj$comb(pi)
    return(obj$objc$adjmat)
}
#----------------------------------------------------

#----------------------------------------------------
#' ##1. Writing additional methods for Updater cf 2.Rnw
#' 
#'  ##2.  match triples in two lists # TODO
#'  obj$listDelta3(4)            # 3rd-order fwd.diff conditioned on col4
#'  # get the same triple
#' triple = as.integer( c(7,8,9))
#' match(triple,obj$listDelta3(3))
#' 
#' ?match
#'  obj$listDelta3(3)[43,]       # 68.68mbits
#'  obj$listDelta3(4)[34,]       # reduced to 46.02mbits by conditioning on vble10
#----------------------------------------------------

#####------------------------------------------------
#' @title jpcalg: calls PCalg for  comparisons.
#' 
#' @param
#' x : data matrix.
#' 
#' gamma : optional scaling of the threshold.
#' 
#' @return
#' Adjacency matrix of graph.
#' 
#' @examples
#' #todo
#' require(pcalg)
#'   suffStat = list(C = cor(X), n = nrow(X)) #  sufficient stats
#'   indepTest=gaussCItest
#'   skel = skeleton(suffStat, indepTest, p=ncol(X), alpha=0.05)
#'   # nodes(skel AT graph) = as(gfobj$Amat,"graphNEL") AT nodes
#'   plot(skel)  
## get skeleton
#'  require(Rgraphviz) ; #  require(gRim); #  ls(package:Rgraphviz)
#'  suffStat = list(C = cor(xdata), n = nrow(xdata)) #  sufficient stats
#'       method = "stable")   #"original","stable.fast"
#'   # skelAmat = as( skel AT graph, "matrix" ) ; sum(skelAmat)/2 # 21 edges
#'   choose(11,2) ; choose(11,3) # 55 165 
jpcalg = function(x,gamma=0.5){
    out.skel =  pcalg::skeleton(   # from pcalg
        suffStat = list(C = var(x), n = nrow(x)),
        indepTest=gaussCItest, p=ncol(x),
        alpha=0.05*exp(3*(.5-gamma)))
    return(as(out.skel@graph,"matrix"))
}#end jpcalg
#------------------------------------------------

#'####------------------------------------------------
#' @title   Normal quantiles
#' @description
#' Ranks a one dimensional vector, and returns quantiles of the normal distribution.
#' @param
#' x : numeric vector.
#' @examples
#'  hist( qnrank(runif(1000)))
#'  \dontrun{xqn = apply(wine,2,qnrank)}
#' 
qnrank=function(x){ # usage 
    n = length(x) 
    qn = qnorm(seq(1:n)/(n+1))
    return(qn[ rank(x ,ties.method ="random")])
}
#'----------------------------------------------------

#'-----------------------------------------------
#' @title xlogx : deals with 0 in  x*log(x)
#'
#' @description 
#' aim: makes x*log(x) cope with x=0
#'
#' @param
#'  x  scalar or numeric, ignores negative entries
#'
#' @value  scalar, log to base 2 ie bits, x 1024 to get millibits
#'
#' @examples
#'  x=c(3,.5) ; xlogx(x) ; xlogx(0); xlogx(-.2);
#'  x = c(.3,.2) ; x = c(.3,.2,0)       ; xlogx(x)  
#'  xlogx(c(2,0)) 
#'
xlogx = function(x){
    nat2mbit <- 2^(10)/log(2)  # depends if a method or a function
    if(!is.numeric(x)){stop}
    x = x[x>0]
    if(is.null(x)){h=0}else{h = nat2mbit*sum(x*log(x))}
    return(h)
}#end xlog
#----------------------------------------------



#####-----------------------------------------------
#' @title  jcpr: cpr for a 2x2 table (input as array)
#' 
#' @param
#' qq : 2x2 table.
#' @examples
#' 
jcpr = function(qq){
    if (!is.array(qq)){ stop('**  not array **') }
    if (!length(dim(qq))==2){ stop('**  not 2 dim **') }
    if (qq[1,2]==0||qq[2,1]==0){ stop('** off diag zeros **') }
    return(qq[1,1]*qq[2,2]/(qq[1,2]*qq[2,1]))
} #end jcpr
#------------------------------------------------

#####-----------------------------------------------
#' @title  jallcpr
#' 
#' @param
#' pp : 2x2x2 table?array.
#' 
jallcpr = function(pp){
    if (!is.array(pp)){ stop('**  not array **') }
    if (!length(dim(pp))==3){ stop('**  not 3 dim **') }
    outsto = matrix(0,nrow=3,ncol=3)
    qq = pp[,,1] ; outsto[1,1] = jcpr(qq) 
    qq = pp[,,2] ; outsto[1,2] = jcpr(qq) 
    qq = apply(pp,c(1,2),sum) ; outsto[1,3] =jcpr(qq)
    qq = pp[,1,] ; outsto[2,1] = jcpr(qq) 
    qq = pp[,2,] ; outsto[2,2] = jcpr(qq) 
    qq = apply(pp,c(1,3),sum) ; outsto[2,3] =jcpr(qq)
    qq = pp[1,,] ; outsto[3,1] = jcpr(qq) 
    qq = pp[2,,] ; outsto[3,2] = jcpr(qq) 
    qq = apply(pp,c(2,3),sum) ; outsto[3,3] =jcpr(qq)
    # evaluate if marg lies in conditional range
    outsto = apply(outsto,1,function(row){
        indic =(row[3] < min(row[1],row[2]))||
               (row[3] > max(row[1],row[2]))
        return(c(row,1*indic))})
    # evaluate simpson :  marg<1 and min>1 or marg > 1 and max <1
    outsto = t(outsto)
    outsto = apply(outsto,1,function(row){
        indic = ((row[3]<1) & (min(row[1],row[2])>1 )) ||
                ((row[3]>1) & (max(row[1],row[2])<1 ))
        return(c(row,1*indic))})
    return(t(outsto))
}# end jallcpr
#------------------------------------------------


#####------------------------------------------------
#' @title  jmcc
#' @param
#' p12 : 2x2 table of true positives, etc
#'
#' @return
#' So called Matthews correlation coefficient.
#' 
#' @examples
#'  #todo, for now a description 
#'  # en.wikipedia.org/wiki/Matthews_correlation_coefficient
#'  #  0.5*sum(abs(Amat1-Amat2)) ;     FP+FN
#'  #         1   0 pred
#'  #  true 1 tp fn
#'  #       0 fp tn
#'  #  +.01 fudge factor
jmcc <- function(p12){
  p12 = as.array(p12,dim=c(2,2))
  TP = p12[1,1]
  TN = p12[2,2]
  FP = p12[2,1]
  FN = p12[1,2]   
  MCC = (TP*TN - FP*FN)/sqrt((TP+FP+.01)*(TP+FN+.01)*(TN+FP+.01)*(TN+FN+.01))
  return(MCC)
}#end jmcc(p=12)
#------------------------------------------------


#####-----------------------------------------------
#' @title jallcomb All combinations.
#'
#' @description
#' Builds integer matrix of all combinations.
#' Auxiliary function used internally to generate and hold all
#' combinations.
#' It is an alternative to looping over the combinations and
#' generating each on the fly.
#' There is a C version. 
#' 
#' @param
#' labels : integer vector, the set of possible entries,
#' @param
#' r : integer, the size of each combination.
#' 
#' @return
#'   Integer matrix,  all combinations, \code{(nCr x r)}.
#' 
#' @examples
#'   \dontrun{ .Call('jallcombCR',c(8,4,5,6,7),3)}
#'   t(combinat::combn(c(8,4,5,6,7),3))  # R version
#' 
jallcomb = function(labels,r){NULL}
#----------------------------------------------------


#####------------------------------------------------
#' @title  lt2full
#' 
#' @param
#' ltdata :
#' 
lt2full = function(ltdata){
 # lt of cov matrix including diagonal ie p(p+1)/2 elements
 # held as a vector
   p = (sqrt(1+8*length(ltdata))-1)/2
   aa <- diag(p)
   aa[upper.tri(aa, diag=TRUE)] <- ltdata
   aa <- aa + t(aa) - diag(diag(aa)) 
}
#------------------------------------------------

#####------------------------------------------------
#' @title  concatCov
#' 
#' @param
#'  xdata1: n1xp data matrix/frame or list n1, S1
#'  xdata2: n2xp data matrix/frame or list n2, S2
#'
#' 
#' @value
#'  (n1+n2)xp data matrix/frame  when both are data matrix/frame
#'  otherwise list: n1+n2, (n1S1+n2S2)/(n1+n2)
#' @examples
#'   xdata1 <<- jsimAR(n=100,p=4,alpha=.4)
#'   xdata1 <<- xdata1*0
#'   xdata2 <<- jsimAR(n=200,p=4,alpha=.1)
#' xdata2 = list(n=5,S=diag(diag(var(xdata1))))
#' concatCov(xdata1,xdata2)

concatCov = function(xdata1,xdata2){
    #'
    if ( (is.matrix(xdata1) || is.data.frame(xdata1)) & 
         (is.matrix(xdata2) || is.data.frame(xdata2)) ){
        xBoth = rbind(xdata1,xdata2)
        return(xBoth)
    }
    if (is.list(xdata1) & length(xdata1)==2){
        n1 = xdata1$n ; S1 = xdata1$S
    }else{
        n1 = nrow(xdata1); S1 = var(xdata1)
    }
    if(is.list(xdata2) & length(xdata2)==2){
        n2 = xdata2$n ; S2 = xdata2$S
    }else{
        n2 = nrow(xdata2); S2 = var(xdata2)
    }
    nc = n1 + n2
    Sc = (n1*S1+n2*S2)/(n1 + n2)
    xBoth = list(n=nc , S=Sc)
    return(xBoth)
}
#------------------------------------------------

####------------------------------------------------
#' @title  ones: a vector of ones.
#' 
#' A vector of ones of given length.
#' @param
#' k : length required.
#' @return
#' Integer Vector of ones.
#' @examples
#' ones(5)
#' 
ones = function(k){
  rep(1,k)
}#end ones
#------------------------------------------------

#####------------------------------------------------
#' @title adjmat2edges edges2adjmat: graph utilities
#' @examples
#'  p=4
#'  adjmat =  matrix(FALSE,nrow=p,ncol=p)
#'  adjmat[2,1]=TRUE ; adjmat[1,2]=TRUE
#'  adjmat[3,2]=TRUE ;  adjmat[2,3]=TRUE
#'  adjmat[4,3]=TRUE ;  adjmat[3,4]=TRUE
#'  adjmat
#'  adjmat2edges(adjmat)
#'  a = adjmat2edges(adjmat)
#'  edges2adjmat(4,a)
#'
#' convert  adj matrix to matrix of edge pairs see/igraph?
#' convert  matrix of edge pairs to adj matrix  see/igraph?

adjmat2edges = function(adjmat){
    adjmat.lt = adjmat[lower.tri(adjmat)]
    p = dim(adjmat)[1]
    edges = expand.grid(1:p,1:p)          # possible edges
    edges = edges[edges[,1]>edges[,2],]   # choose lower triangle
    rownames(edges) = NULL
    return(edges[adjmat.lt,])
}#end adjmat2edges

edges2adjmat = function(p,edges){
    adjmat = matrix(FALSE,nrow=p,ncol=p)
    null =  apply(edges,1,function(e){  #cat(e,"\n")
        adjmat[e[1],e[2]] <<- TRUE  # need the global assign
        adjmat[e[2],e[1]] <<- TRUE
    })
    return(adjmat)
}#end edges2adjmat
#----------------------------------------------------

#####-----------------------------------------------
#' @title Random edge adjacency matrix
#'
#' @description
#' Auxiliary function to handle the mat="random" option in
#' the  setAdjMat method of the AdjMat class.
#'
#' @param
#' p : the dimension of the matrix.
#' @return
#' A symmetric Boolean matrix.
#' 
#' @examples
#'  require(reta)
#'  xdata <- jsimAR(n=400,p=6,alpha=.4)
#'  obj = AdjMat$new(xdata)
#'  obj$randomEdge
#'  obj$setAdjMat(mat="random")
#'  obj$adjmat
#'  obj$randomEdge(5)
#'  randomEdge <- function(p) { NULL }
#'
randomEdge = function(p){
    s = matrix(sample(c(0,1),size=p*p,replace=TRUE),nrow=p)
    s[lower.tri(s,diag=FALSE)] = 0
    s = s + t(s)
    diag(s) = 0
    return(matrix(as.logical(s),p))
} #end randomEdge
#-----------------------------------------------


#####------------------------------------------------
#' @title wrongEdges
#'
wrongEdges = function(amat, adjmat.true){""
    p = dim(adjmat.true)[1] ; totedges = p*(p-1)/2
    nedges = sum(adjmat.true)/2
    we = amat - adjmat.true
    FP = sum(we>0)/2
    FN = sum(we<0)/2  
    TP = sum(amat*adjmat.true)/2
    TN = totedges - (FP+FN+TP)
    out = c(FP,FN,TP,TN)
    names(out) = c("FP","FN","TP","TN")
    return(out)
} #end wrongEdges
#'-----------------------------------------------------

####------------------------------------------------
#'  oldjwronge: counts number of edges differ comparisons.
#' 
#' @param
#' Amat1 : adjacency matrix
#'
#' Amat2 : adjacency matrix.
#' 
#' mcc : if True returns the table.
#' 
#' @return
#' Adjacency matrix of graph.
#' 
#' @examples
#' #todo
oldjwronge = function(Amat1,Amat2,mcc=FALSE){
#'  check ?nargs mcc = TRUE
#' jwronge(Amat1,Amat2,mcc=TRUE) jwronge(Amat1,Amat2)
    lt = lower.tri(Amat1)
    if (mcc){
    p12 = table(as.data.frame(cbind(Amat1[lt],Amat2[lt])))
    return(p12) ; #return(jmcc(p12)) ;
  }else{
    return( sum(abs(Amat1[lt]-Amat2[lt])) ) ;
  }
}#end oldjwronge
#------------------------------------------------

#####-----------------------------------------------
#' @title getKwayArray
    
#' aim: example of a 2^k table
#' input: k number of variables
#' output: 2^k table


#'  pp= getKwayArray(4) ; dimnames(pp)
#'  as.data.frame(pp)
getKwayArray=function( k=4){
  V = 1:k
  dimV = rep(2,k)
  p = 1:prod(dimV)  
  p = p/sum(as.numeric(p))

  # labels
  y = matrix(c(0,1),nrow=2,ncol=k)
  pp = table(data.frame(y))
  pp = 0*pp+array(p,dim=dimV)  # names(dimnames(pp))

  return(pp)
}#end  getKwayArray
#-----------------------------------------------


####-------------------------------------------------
#' @title   arr2df
#' @examples
#'    pp = array(1:8,dim=c(2,2,2,2)) 
#'    pp = array(1:6,dim=c(2,3)) 
#'    arr2df(n=sum(1:16),pp)
#'    # arr2df(96, as.array(obj$xtab)) ; str(obj$xtab)
#'    length(dim(pp))  
#'    x =   arr2df(sum(1:8),pp) ; class(x)  
#'    arr2df(8,pp)  # will give result but not right as need to choose good n
#'    as.data.frame(pp) 
#'    as.data.frame(as.table(pp))
#'
arr2df = function(n,pp){"
aim: Convert array of frequencies to data matrix
     with vector of frequencies
     and then to replicate rows of the matrix
     to match each frequency. Approximate due to rounding.
Note:     arr2df: need to deal with >3 levels; and with rounding 
input: n, upper limit to the number of rows in the df
input: pp  a table, mway array of proportions.
global modified: xdf,  n by p data frame
"
            if (!is.array(pp)){ return('input is not an array/table') }
            pp = pp/sum(as.numeric(pp))
            out = as.data.frame(as.table(pp))  #  archive  out = arr2datmat(pp)
            x = out[,1: length(dim(pp))]   # as.numeric(x)  str(x) dimnames(x)
            if (x[1,1]=="A"){
                x = apply(x,c(1,2),function(item){
                    if(item=="A"){item=0}else{if(item=="B"){item=1}else{item=2}}
                })}
            freq = out$Freq
            freq = as.integer(n*freq)
            # if (!(sum(freq)>1)){ return('input is not a count') }
            xx = sapply(1:nrow(x),
                        function(i){x[rep(i, each=freq[i]),]},simplify=FALSE)
            #  modify global
            xdf <- as.data.frame(do.call("rbind", xx))
            return(xdf)
        }
# end arr2df
#-------------------------------------------------

#####------------------------------------------------
#' @title   fakeCov: fakes data with the exact covariance
#'
#' Generates a (p+1) x p data(!) matrix whose empirical variance
#' is identical to the theoretical variance.
#' Uses the Cholesky decomposition of the theoretical variance,
#' and augments to reverse engineer the mean adjustment made by computing the
#' empirical variance.
#' 
#' @param
#' vxx : theoretical variance matrix.
#' 
#' @examples
#'
#'    X = jsimAR(n=200,p=4,alpha=.4)  # data matrix for testing
#'    vxx = cor(X) 
#'    xdata = fakeCov(vxx)
#'    dim(xdata)
#'    all.equal(cor(xdata),vxx)
#' 
fakeCov = function(vxx){
  p = nrow(vxx)
  lambda = sqrt(p+1)+1    #  lambda^2-(p+lambda)^2/(p+1)
  B = solve(chol(solve(vxx)))   # round(100.*solve(B%*%t(B)))  # good
  # augment matrix to get rid of mean in stats::cor
  return( t(cbind(B,lambda*apply(B,1,mean))) )
}# end fakeCov  
#------------------------------------------------

#------------------------------------------------

#======================================================================
