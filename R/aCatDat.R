###############################################################################
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.

# Author: Joe Whittaker 
###############################################################################

#####-----------------------------------------------
#' @title  CatDat reference class for categorical data entry.
#'
#' @description
#'
#' A method in this RC reads the data as a data matrix or as a data
#' frame.  A method is given that calculates an entropy from the data
#' frame \code{getHcat}.  The method is an application of the analytic
#' formula for the entropy
#' 
#' The methods are applied to compute the matrix of mutual
#' informations (MIs) by tabulations on the data frame, and secondly,
#' to compute an arbitrary entropy from the data frame \code{getHcat}.
#'
#' 
#' The fields of this RC include n = "numeric", p = "integer", xdf =
#' "data.frame", xtab = "array", xlevels = "integer".
#'
#' @seealso \code{\link{InitData}}.
#' 
#' #-----------------------------------------------
#' @examples
#' ## source("retaPkg/R/aCatDat.R")
#' ##          
#'  library(reta)
#'  set.seed(1234)
#'  
#'  ##  Enter data from rectangular nbyp dataframe
#'  xdata <- data.frame(jsimBinDG(n=4000))
#'  obj <- CatDat$new()
#'  str(obj)  # essentially shows fields held by AdjMat and HDlis
#'  obj$readCatDat(xdata)
#'  obj$n
#'  obj$p           # dimension
#'  head(obj$xdf)        
#'  obj$xlevels
#'  obj$getHcat(c(2,4))   # 1759
#'  obj$getHcat(2) 
#'  out = obj$xdf2mimat() ;   out 
#'  tuple = c(2,4) 
#'  table(obj$xdf[, tuple])  
#'  tuple=as.integer(3)
#'  table(obj$xdf[, tuple])  # no heading : code needs workaround
#'  
#'  
#' @export
#'

####-----------------------------------------------
CatDat <- setRefClass("CatDat",
    fields = list(
        n = "numeric",
        p = "integer",
        xdf = "data.frame",
        xtab = "table",   # would like this to be  array?
        xlevels = "numeric",
        nat2mbit = "numeric"
    ),
#-----------------------------------------------
    methods = list(
#-----------------------------------------------


#-----------------------------------------------
        readCatDat = function(xin){
            "
Reads the data from either a data frame or a matrix,
but not from an array or a table.
It determine the dimensions of the data and computes the
number of levels (attained values) of each variable.

The functions modifies the global fields: 
  p, integer, dimension of variables,
  xdf,  data frame of cases by variables,
  xtab,  array giving counts of variables
  xlevels,  vector of number of levels
  nat2mbit,   numeric constant, conversion rate of natural entropies to millibits.
"
            if (!(is.matrix(xin) || is.data.frame(xin))) {return("error1\n")}
            n <<- as.integer(nrow(xin))
            p <<- as.integer(ncol(xin))
            xdf <<- as.data.frame(xin)
            xlevels <<- apply(xin,2,function(col){
                lev = unique(col)
                length(lev)
            })
            nat2mbit <<- 2^(10)/log(2)
        }
,#end readRawData
#-----------------------------------------------




#-----------------------------------------------
        jallPairs = function(k){""
            if (!(k>0)){"k not positive"; return() }
            a = cbind(gl(k,k),gl(k,1))
            return(a[a[,1]<a[,2],])
        }
,#end jallpairs
#-----------------------------------------------
 
# this is where xlogx used to be 

#-----------------------------------------------
        getHcat =  function(tuple){"
 bottom up computation 
does not check if h(tuple) has already been computed
requires: xdf  nxp data frame,
protects against table(1dim df) where loses name
input: tuple, integer vector of dimensions
output: h, numeric entropy in millibits
calls xlogx
"
#' p=obj$p ; xdf =  obj$objb$xdf ; getHcat=obj$getHcat
#' tuple=c(3,4)  # tuple=5   tuple=c(1,2)  getHcat(tuple)

            tuple = as.integer(tuple)
            x = xdf[,tuple]   # cannot condense these  .self objb$
            # x = xdf[,tuple]   # cannot condense these  .self$
            x = as.data.frame(x)    # two lines             str(x)
            pp = table(x)
            pp = as.numeric(pp)
            pp = pp/sum(pp) 
            h =  -xlogx(pp)
            # hlis[[ltuple]][[hash]] <<- list(ss=tuple,h=h)
                                        # hlis storage gTriples
            return(h)
        }#end getHcat
,
#-----------------------------------------------
 
#-----------------------------------------------
        xdf2mimat =  function(){"
  mimat  numeric matrix,   mutual information between variable pairs
  measured in millibits.
  kmax integer, maximum order of the algorithm; the order determines
the size of conditioning sets, and the magnitude of the list storage requirements.
calculates mimat for AdjMat RC
would like to do this in lapply as in Cormat by global def 
"
#'
#' p=obj$p ; xdf =  obj$xdf ;  getHcat=obj$getHcat ; jallPairs=obj$jallPairs
#' obj$xdf2mimat ; obj$mimat
#' source("retaPkg/R/aCatDat.R")
#'  obj <- CatDat$new() ;  obj$readCatDat(xdata)  
#'            
#' one dimensional entropies
#'  p  =  obj$p
#'  xdf =  obj$xdf
#'  getHcat = obj$getHcat  
#'  jallPairs = obj$jallPairs
#   = obj$
             p = dim(xdf)[2]
            # cat("when1\n")
            # xlevels <<- as.integer(rep(0,p))  # todo
            entList = lapply(1:p,function(node){ # gives list output
                # h = .self$getHcat(item) 
                h = getHcat(node)    # node = as.integer(4)
            })
            # reorganise as vectors    entList[[2]][2]
            h1vec = unlist(entList)
            
            #' two dimensions 
            pairs = jallPairs(p)
            npairs = nrow(pairs)
            milist = lapply(1:npairs,function(item){ # item=2  output list
                pair = pairs[item,]  # 2 dim  # pair = c(1,2)#print(pair)
                h = getHcat(pair)   # old NaN due to xlogx at 0
                d = h  - h1vec[pair[1]] - h1vec[pair[2]]
                list(h,d)
            })
            h12vec = unlist(milist)[2*(1:npairs)-1]    #
            d12vec = unlist(milist)[2*(1:npairs)]    #

            mimat = 0*diag(p)   # reorganise as matrix
            mimat[lower.tri(mimat, diag=FALSE)] = -d12vec
            mimat = mimat + t(mimat)
            dimnames(mimat) = c(nodenms,nodenms)
  
            return(list(h1vec=h1vec,h12vec=h12vec,d12vec=d12vec,mimat=mimat))
        }# end  xdf2mimat
,
#-----------------------------------------------


#'-----------------------------------------------------
develop = function(){"" # no comment stops it entering web

#-----------------------------------------------
        readArray = function(n,pp){"
NOT IN USE
 Reads the sample size and the mway array of counts.
 Evaluates the dimension of the array.
"            
#'        array(1:8,dim=c(2,2,2))
#'        array(1:6,dim=c(2,3))
#'        as.table(array(1:6,dim=c(2,3)))
#'        ?array  ?table
            n <<- as.integer(n)
            xtab <<- as.array(pp)
            # xlevels <<- dim(xtab)  # class()
            p <<-  length(xlevels)
        }
#end readArray
#-----------------------------------------------
 

#-----------------------------------------------
        getHmway =  function(tuple){"  TODO
  requires: xtab,  k-way array,
  input: tuple, integer vector of dimensions
  output: h, numeric entropy in millibits
  calls xlogx
"
#' pp = array(1:8,dim=c(2,2,2)) ;  pp=pp/sum(as.numeric(pp))
#' getH(c(1,2,3))
#' tuple=as.integer(c(1,3)) ;  getH(tuple) ; class(tuple)
#' getH(2)
#' getH(integer(0)) # error 
#'
#'   x=c(.1,.2,.3,.4)
#'   pp=array(x,dim=c(2,2))
#'   xlogx(pp)   # -1891
#'
            tuple = as.integer(tuple)
            h = -xlogx(apply(xtab,tuple,sum))  # top down pway array
    #  bottom up from nxp xdf 
    # hlis[[ltuple]][[hash]] <<- list(ss=tuple,h=h)  # hlis storage 
            return(h)
        }
#end getHmway
#-----------------------------------------------
}# end develop
#'-----------------------------------------------------


   
#####-----------------------------------------------
) )#end methods CatDat
#####-----------------------------------------------
#####-----------------------------------------------
#' @title  All pairs of integers.
#'
#' @description
#' An auxiliary function used internally
#' to produce  all  pairs in the integers \code{1:k}.
#' 
#' @param  k  number of integers.
#'
#' @return
#' Integer matrix  \code{(kC2 x 2)}.
#' 
#' @examples
#'  xdata <-  jsimAR(n=100,p=6,alpha=.4)
#'  obj <- HDlis$new(xdata)
#'  obj$jallPairs(4)    
#'  system.time(obj$jallPairs(1000)) 
#'  system.time(combinat::combn(1000,2))  #  somewhat slower.
jallPairs <- function(k){ NULL }

#####-----------------------------------------------
