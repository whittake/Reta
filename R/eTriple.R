#####-----------------------------------------------
#####
#' @title Triple reference class
#'
#' @include   dHDlis.R
#' 
#' @description

#' The class contains the methods that find the entropy and the 3-way
#' interaction associated with a triple using one additional entropy
#' evaluation and a conditional MI evaluation.  The difference,
#' \code{delta}, is stored in \code{dlis}.
#' 
#' @details
#'
#' The method \code{getHfly} computes the entropy, while 
#'  \code{getDelta} computes the interaction. 
#'
#' The private variable \code{ltuple} is determined by the length of
#' the tuple or of tuple-pair and is used to differentiate lists.
#'
#' @examples
#'
#'  require(reta)
#' help(package='triples')
#'   # require(testthat)
#'   #  source("retaPkg/inst/tests/testabcfgh.R")
#'  set.seed(1234)
#'  xdata <- jsimAR(n=400,p=8,alpha=.4)
#'  #  xdata <- data.frame(jsimBinDG(n=400))
#'  obj <- Triple$new(xdata)
#'
#'  obj$initHDlis() ; obj$getallH12()
#'  obj$dlis[[2]] ; obj$mimat 
#'  obj$hlis[[2]] ; unlist(obj$hlis[[1]])
#'
#'  triple = c(1,2,3)
#'  obj$getH(triple)
#'  obj$getHfly(triple)
#'  obj$hlis[[length(triple)]]
#'  obj$dlis[[2]]
#'  obj$listHD(opt='H',kord=2) ;obj$listHD(kord=3);obj$listHD(opt='H',kord=3)
#'
#'
#'  obj$getDelta(triple) #  RC-Triple: warning 
#'  obj$listHD(kord=3)
#'  ## increments list at each new call:
#'  obj$getHfly(c(1,3,4)) ;  length(obj$hlis[[3]])
#'  obj$getHfly(c(1,2,4)) ;  length(obj$hlis[[3]])
#'  obj$getHfly(c(1,3,4,2)) ; obj$hlis[[4]]
#'  obj$getHfly(sort(c(1,3,4,2))) ; obj$hlis[[4]]  # default 
#'
#'  obj$setAdjMat("complete") 
#'  obj$getDelta(c(1,2,3)) 
#'  obj$dlis[[3]]             # includes both conditional d23|1 and d123
#'  length(obj$hlis[[3]])
#'
#'  tuple  = c(2,4,5)
#'  obj$getD(tuple)  # all entropy interactions
#'  obj$mimat[tuple,tuple]  # 2-way interactions agree with mimat
#'  obj$getDelta(tuple)    # so does 3-way interaction
#'
#'##  kord=4 is  not done with this version yet
#'  obj$getHfly(c(1,2,3,4)) ; obj$hlis[[4]]  # ok 
#'  obj$getDelta(c(1,2,3,4)) #  RC-Triple: warning 
#'  obj$hlis[[4]]  # ok 
#'  obj$hlis[[3]]
#'  obj$dlis[[3]]  # contains 1234 34|12
#'  obj$dlis[[4]]  # empty
#'
#' ## initialisation  requires adjmat, update options, and MI statistics
#'  set.seed(1234)
#'  xdata <- jsimAR(n=400,p=3,alpha=.4)
#'  obj <- Triple$new(xdata) ; str(obj)   # start again 
#'  obj$mimat 
#'  obj$micut ; obj$setMIcut() ; obj$micut 
#'  obj$adjmat+0 ; obj$setAdjMat() ; obj$adjmat+0 
#'   obj$score12way() ; obj$score[2]
#'  obj$dlis  ;  obj$initHDlis() ; obj$dlis
#'  obj$dlis[[2]] ;  obj$getallH12() ; obj$dlis[[2]]
#'
##----------------------------------------------
Triple <- setRefClass("Triple",
    contains = list("HDlis"),
    methods = list(
        
##----------------------------------------------
        getHfly = function(tuple){" 
Checks if the entropy h for the tuple is available, and if not
calls the function getH() to evaluate and store.
"
            tuple = as.integer(tuple)
            ltuple = length(tuple)
            if (ltuple<1){ return(0) } # entropy of empty set to zero 
            hash = digest(tuple)
            if (!is.null(hlis[[ltuple]][[hash]])){
                # cat("\n fly ",tuple," look up h ")
                return(hlis[[ltuple]][[hash]]$h) # $h fails  [[2]] works
            }else{
                # cat("\n geth ",tuple," compute h now  ")
                h = .self$getH(tuple)
                hlis[[ltuple]][[hash]] <<- list(ss=tuple,h=h)
                return(h)
            }
        }# end getHfly
,  
##----------------------------------------------



##----------------------------------------------
        getDelta = function(triple){"
The 3rd-order interaction is evaluated for the triple
(tuple of length 3) from its entropy and the entropies of its subsets.
Specifically it is calculated as the difference between
the marginal MI of a pair in the triple, and the conditional MI
of the pair given the 3rd element of the triple.
It is stored with the triple, with any conditioning set, in dlis.
Requires arguments: triple;
Requires globals:  dlis;
Calls methods: getHfly,;
Checks are made along the way to ensure lower order entropies and 
interactions are available.
"
#' triple=c(1,2,3) ; dlis=obj$dlis ; getHfly=obj$getHfly
            ltuple = length(triple)
            if (!(ltuple==3)){
                return(cat("RC-Triple: argument not of length 3"))}
            triple = as.integer(triple)
            hasht = digest(triple)
            # check if delta for required subsets  already there
            if(!is.null(dlis[[ltuple]][[hasht]])){
                #cat("RC-Triple: delta available\n")
                return(dlis[[ltuple]][[hasht]]$d)}
            #  take first pair from triple, triple is sorted 
            pair =  as.integer(triple[1:2])   # is.integer(pair)
            hashp = digest(pair)
            mip = -dlis[[2]][[hashp]]$d     #  mi(i\indep{}j)
            # get mi(i\indep{}j|k)= -h[ijk]+ h[ik]+ h[jk]- h[k]
            k =   as.integer(triple[3])  # k=Ak for the cond set 
            mip.k = -getHfly(triple) + 
                     getHfly(c(pair[1],k)) +
                     getHfly(c(pair[2],k)) -
                     getHfly(k)
## disable writing cond test to dlis 
##            hashp.k = digest(list(pair,k)) # list storage, different to triple
##            dlis[[ltuple]][[hashp.k]] <<- list(ss=list(pair,k),d=-mip.k)  #<<-
            # evaluate delta3
            delta3 = mip - mip.k                        
            dlis[[ltuple]][[hasht]] <<-  list(ss=triple,d=delta3) #<<- 
            return(delta3)
        }#end getDelta3
,
##----------------------------------------------



##----------------------------------------------
allH = function(tuple){"
Computes the entropies of all subsets of the tuple.
Input: integer tuple;
Output: numeric vector of entropies in standard order
(eg . 1 2 12 3 .. 123) of the sorted tuple.
Details:
Calls getHfly to evaluate the entropy.
Uses base::intToBits to generate the indicator vectors.
"
## 2nd go : generate row wise using intToBits
## tuple = c(3,4,5,6)
## tuple =   numeric(0)
## tuple = c(4,7,8)
##   allH(tuple)
##  obj$getHfly(tuple)
##  ind = c(TRUE,TRUE)  ; tuple = c(4,7) ; tuple[ind]  ! magic
##  intToBits(0)
    k = length(tuple)
    if (k<1){return("nothing to do")}
    sapply(0:(2^k-1), function(i){       # i=3 i=0 i=7
        ind = as.logical(intToBits(i))[1:k]  # first k elements of 32bit vector
        getHfly(tuple[ind]) 
    })
}# end allH
,
##-------------------------------------------------



##-------------------------------------------------
getD = function(tuple,all=FALSE){"
Mobius transform of the entropy on the powerset of the tuples
gives the entropy interactions upto k, the length of the tuple.
Calls allH
Input:  integer tuple;
Output: numeric vector of entropy interactions in standard order
(eg . 1 2 12 3 .. 123) of the sorted tuple.
"
## generate col wise using floor, same as remainder(n,k)
## tuple =  c(1,2,3)
## tuple =  c(2,3,4)
##   obj$allH(tuple)
##    obj$getD(tuple,all=TRUE)
## allH = obj$allH
    k = length(tuple)
    if (k<1){return("nothing to do")}
    D = allH(tuple)
    all2 =0:(2^k-1)        # straightforward sequence
    for (i in 1:k){
        xind =  (all2 - 2^i*floor(all2/2^i)) >= 2^(i-1)
        D[xind] =  D[xind] - D[!xind] 
    }
    if (!all){
        allbutone = 1:(length(D)-1)
        return(D[- allbutone])
    }else{
        return(D)
    }
}# end getD
,
###----------------------------------------------------


###----------------------------------------------------
develop = function(){""
}# end develop
###----------------------------------------------------


))

######----------------------------------------------
