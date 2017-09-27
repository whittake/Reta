#####-----------------------------------------------
#####
#' @title  InitData reference class for  data entry 
#'
#' @include aCatDat.R   aCorMat.R 
#'  
#' @description
#'
#' The RC reads the input data, determines the type of data
#' (continuous, categorical), chooses the appropriate function to
#' compute the entropy, and uses this to initialise entropy storage and
#' the matrix of pairwise MIs.  
#'
#' @details
#'
#' The constructor method (\code{initialize}) of this class determines
#' the data type of the input and reads the data accordingly.  The
#' format of the data is a data matrix or a list of sufficient
#' statistics.  The current data types handled invoke either a CorMat
#' object or a CatDat object. This choice of type also determines the
#' method \code{getH} to compute the entropy of any subset of
#' variables from either the correlation matrix or from the data frame
#' of observations.
#'
#' The method creates the pairwise MI matrix \code{mimat} and one and
#' two dimensional entropies.
#'
#' Among the reference classes of the triples package \code{InitData}
#' is the only one with an \code{initialize} method: consequently this
#' is called whenever an object from a superclass is invoked.
#'
#' @field  n   numeric, number of independent units;
#' @field   p   integer, number of variables;
#' @field   nat2mbit   numeric, conversion ratio from natural logs to millibits;
#' @field   h1vec   numeric, container for entropies;
#' @field   mimat   matrix, pairwise MI measures in millibits;
#' @field   getH   function, method to extract entropy from sufficient statistic.
#'
#' @seealso \code{\link{CorMat}}, \code{\link{CatDat}}.
#'
#' 
#' @examples
#' 
#' ## Input continuous data matrix 
#'  require(reta)
#'  help(package="reta")
#'  obj <- InitData$new()
#'  str(obj)
#'  set.seed(1234)
#'  xdata <- jsimAR(n=400,p=4,alpha=.4)  
#'  length(unique(xdata[,1]))
#'  obj <- InitData$new(xdata)
#'  obj$dimnms
#'
#' ## Inspect elements of the object
#'  obj$nat2mbit             # conversion nat2mbit
#'  obj$getH
#'  obj$getH(c(1,3))
#'  obj$xcor
#'  names(InitData$fields())
#'  obj$mimat
#'  obj$h12vec  ; obj$h1vec ; obj$d12vec
#'  CorMat$methods("getHcor")
#'
#' ## Input categorical data
#'  xdata <- data.frame(jsimBinDG(n=400))
#'  tail(xdata)
#'  obj <- InitData$new(xdata)
#'  str(obj)
#'  obj$getH
#'
#' ## Inspect elements of the object
#'  obj$n
#'  obj$p           # dimension
#'  obj$nat2mbit
#'  obj$getH(c(1,3))
#'  xlogx    # obj$xlogx  fails 
#'  head(obj$xdf)  
#'  obj$xlevels
#'  obj$h1vec
#'  obj$h12vec
#'  obj$mimat
#'  obj$xdf2mimat()
#'
#' ## Check labels
#'  require(reta)
#'  xdata <- reinisfn()  ; head(xdata)  ; class(xdata)  
#'  obj <- InitData$new(xdata)
#'  str(obj)
#'  obj$xlevels
#'  obj$nodenms
#'  obj$mimat
#' 
#' @export   # do i need this 
#'  
#-----------------------------------------------
InitData <- setRefClass("InitData",
    contains = list("CorMat","CatDat"),
    fields = list(
        nodenms = "list",
        h1vec = "numeric",
        h12vec = "numeric",
        d12vec = "numeric",
        mimat = "matrix",
        getH = "function"
    ),
#-----------------------------------------------
    methods = list(
#-----------------------------------------------
        
#-----------------------------------------------
       initialize =  function(xin=0){
           "
Determine the type of the data.
Read the dimnames.

Cases:
type in (2,3)  categ df data frame/matrix
type     > 3   cont df data frame/matrix
type    == 0   cont suff stat  list(n, S)

"
#' obj$readCatDat
           type=0  # default 
           # check unique values of first column
           if (is.matrix(xin)||is.data.frame(xin)){ # data frame/matrix
               type=length(unique(xin[,1]))
               nms = dimnames(xin)[2]
               if( !is.null(nms)){ nodenms<<- nms}
           }
           # cat(type,"type\n")
           if (type>1 && type<4){   # categ if have  2 or 3 levels
               readCatDat(xin)
               # cat("Read categorical n") 
               getH = .self$getHcat       # must have .self
                    # .self$callSuper(getHcat)  fails 
               out = xdf2mimat()
            }else{   # type =0,>3   continouos
                if (type>3 ){     
                    readContDat(xin)  #  data matrix/frame
                    # cat("Read continuous data frame/matrix\n") 
                }else{  #  type==0  taken to be cont, suff stat=cormat
                    readVarMat(xin)
                    # nodenms evaluated in readVarMat
                    # if(!is.null(xcor))nodenms{ <<- dimnames(xcor)} fails
                    # cat("Read continuous data suff stat\n") 
                 }
                getH = .self$getHcor
                out = xcor2mimat()
            }
           getH <<- getH
           h1vec <<- out$h1vec
           h12vec <<- out$h12vec
           d12vec <<- out$d12vec
           mimat <<- out$mimat
        }
, # end initialize

#-----------------------------------------------
 
#'-----------------------------------------------------
develop = function(){""
}# end develop
#'-----------------------------------------------------

   
#####-----------------------------------------------
) )#end methods InitData
#####-----------------------------------------------
