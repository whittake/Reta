###############################################################################
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.

# Author: Joe Whittaker 
###############################################################################

#####-----------------------------------------------
#####
#' @title  CorMat reference class for continuous data
#'
#' @description 
#' 
#' A method in this RC read the data matrix and compute the
#' correlation matrix or read a list of sufficient statistics, and
#' extracts the correlation matrix, and the dimension of the data.  A
#' method is given that calculates an entropy from the correlation
#' \code{getHcor}.  The method is an application of the analytic
#' formula for the entropy of the multivariate Normal distribution.
#'
#' The fields of this RC are n = "numeric", p = "integer", xcor =
#' "matrix".
#'
#' @seealso \code{\link{InitData}}.
#'
#'
#' @examples
#' 
#' ## Input rectangular data matrix
#'  source("retaPkg/R/aCorMat.R") 
#'  library(reta)
#'  set.seed(1234)
#'  xdata <- jsimAR(n=400,p=4,alpha=.4)
#'  # xdata <-  matrix(rnorm(20),nrow=5,ncol=4)
#'  obj <- CorMat$new()
#'  str(obj) ; obj$xcor2mimat
#'
#'  ## make data input explicit
#'  obj$readContDat(xdata)
#'  str(obj)
#'  obj$p           # dimension
#'  obj$xcor     
#'  obj$xcor2mimat()
#'  -obj$nat2mbit * log(1-obj$xcor^2 + diag(rep(1, obj$p)) )/2
#'  # obj$mimat 1st row
#'  #   0.0000 129.7716  29.334   5.403   1.348   2.0906   0.3964
#'  tuple = c(1,4) 
#'  obj$getHcor(tuple)  #  -0.2978
#'  
#' ## OR input sufficient statistics
#'  #  n = nrow(xdata)
#'  xin = list(n=obj$n, S= var(xdata))
#'  xin
#'  obj$readVarMat(xin)
#'  obj$xcor2mimat()
#'
#'  obj <- CorMat$new()
#'  obj$readVarMat(list(n=14,S=diag(1:4)))
#'  str(obj)
#'  obj$p           # dimension
#'  obj$xcor
#'  obj$xcor2mimat()
#'  obj$nat2mbit    # conversion factor to millibits 1477 
#'
#' @export
#'
#' 
####-----------------------------------------------
CorMat <- setRefClass("CorMat",
                      
#-----------------------------------------------
    fields = list(
        n = "numeric",
        p = "integer",
        xcor = "matrix",
        nat2mbit = "numeric"
    ),
#-----------------------------------------------

#-----------------------------------------------
    methods = list(

#-----------------------------------------------
        readContDat = function(xin){
            "
Reads data as a matrix and finds the correlation matrix.
"
            if (!(is.matrix(xin) || is.data.frame(xin))) {return("error3\n")}
            #   if (min(svd(xin)$d)<1.0e-10){stop('Data singular')}
            n <<- as.integer(nrow(xin))
            p <<- as.integer(ncol(xin))
            xcor <<- cor(xin)
            nat2mbit <<- 2^(10)/log(2)
        } #end readRawData
, 
#-----------------------------------------------

#-----------------------------------------------
         readVarMat = function(xin){
            "
Reads the sample size and variance matrix as the list(n=n,S=S)
and finds the correlation matrix.
"
            if (!(is.list(xin)&&(length(xin)==2))) {return("error2\n")}
            n <<- as.integer(xin$n)
            p <<- as.integer(dim(xin$S)[1])
            xcor <<- cov2cor(xin$S)
            #  dimnms <<- dimnames(xin$S) fails here
            # cat("Read sufficient statistic\n")
            nat2mbit <<- 2^(10)/log(2)
       } #end readVarMat
,
#-----------------------------------------------

#-----------------------------------------------
getHcor =  function(tuple){
            "
An entropy is associated with a tuple (ordered subset) of nodes in
the graph.
The method evaluates the entropy for a given tuple and stores both
the tuple and the entropy.
It  implements  the correlation matrix formula appropriate to the
Gaussian distribution.
It makes a call to RccpEigen to evaluate the appropiate determinant.
Invocation requires the argument
tuple, an integer vector of the nodes forming the tuple.

A related method {getHfly} first checks whether the entropy
has been previously calculated, and organises storage.
"
#'  nat2mbit=obj$nat2mbit ; xcor=obj$xcor ; tuple = c(1,2,3)
#'  nat2mbit*log(det(xcor[tuple,tuple] ))/2  # Rversion
            if (max(tuple)>p){return("tuple invalid\n")}
            h = nat2mbit*.Call('jlogdetCR',xcor[tuple,tuple] )/2  # .self$
            return(h)
        }# end getHcor
,  
#'-----------------------------------------------


#'-----------------------------------------------
xcor2mimat =  function(){
    "
Uses the analytic formulae for the Gaussian entropy to
compute the matrix of pairwise MIs from the correlation matrix.
"
#' mimat=obj$mimat
            # cat(p,"here5\n")
            p = dim(xcor)[1]  # work around because cannot access global p
            # cat(n,"here6a\n")
            h1vec <-  rep(0,p)   # .self$p
            # cat(n,"here6b\n")
            mimat <- -nat2mbit*log(1-xcor^2+diag(rep(1,p)))/2
            # cat(n,"here7\n")
            h12vec <- -mimat[lower.tri(mimat)]
            return(list(h1vec=h1vec,h12vec=h12vec,d12vec=h12vec, mimat=mimat))
        }# end xcor2mimat
,
#'-----------------------------------------------


#'-----------------------------------------------------
develop = function(){""
}# end develop
#'-----------------------------------------------------


) )# end CorMat

#####-----------------------------------------------
