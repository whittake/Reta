#####-----------------------------------------------
#####
#' @title AdjMat reference class to hold  adjacency matrix manipulations.
#'
#' @include bInitData.R
#'
#' @description
#'
#' This class holds the graph as a Boolean adjacency matrix;
#' including options to set and update the matrix.
#' 
#' @details
#'
#' It determines the
#' MI cut point for thresholding the observed MIs which initialises
#' most graphical search procedures.  The method \code{removeEdge} is
#' held in this reference class.  The adjacency matrix can be reset at
#' any stage, but a knock on consequence is that additional entropy
#' evaluations corresponding to the new graph may be required.
#' 
#' Each edge of the graph is listed together with its MI by the method
#' \code{initEtovisit()}; the next weakest edge is returned by a call to
#' \code{nextEdge()} and expunged from \code{etovisit}.  The potential
#' conditioners of the edge are the common neighbours of the edge
#' nodes; these are evaluated by a call to \code{getCset()} and the set
#' is returned in rank order of strength, strongest first.
#' 
#' The method \code{testEdgeLoop()} provides a dummy loop over the graph
#' indexed by edge and then its potential conditioning set. The edge
#' and the elements of the set are printed.  The method
#' \code{score12way()} initialises the (triangular) score of the current
#' graph.  
#' 
#' @field  adjmat, logical symmetric matrix, of dimension \code{p}.
#' @field  initadjmat, character, designating initialising \code{adjmat}.
#' @field  etovisit, matrix of edges in current graph and MIs.
#' @field  score[2], numeric, the (2way interaction) score of the graph.
#' @field  bic, numeric, the BIC multiplier to penalise the score.
#' @field  micut, numeric, cut point for thresholding the MIs.
#'
#'
#' @examples
#'
#'  library(reta)
##  source("retaPkg/R/cAdjMat.R")
#'  obj <- AdjMat$new()
#'  str(obj)
#'  # start again with data
#'  set.seed(1234)
#'  xdata <- jsimAR(n=400,p=6,alpha=.4)
#'  obj <- AdjMat$new(xdata)
#'  str(obj) ; obj$nodenms
#'  obj$mimat ;  obj$h12vec   # available from initialization 
#'  
#' ## The field adjmat holds the adjacency matrix:
#'  obj$setAdjMat("complete") ;   obj$adjmat +0
#'  obj$setAdjMat("random")  ;  obj$adjmat +0
#'  AdjMat$methods()
#'
#' ## Choose a random adjacency matrix, and incorporate in object
#'  amat <- jsimGraph(p=5,edgeprob=.5)  ; amat+0
#'  obj$setAdjMat(amat)  ; obj$adjmat +0
#'
#' ## Set the MI cutpoint
#'  obj$nat2mbit
#'  obj$n
#'  obj$setMIcut()
#'  obj$micut
#'  obj$bic
#'
#' ## Initialise the adjacency matrix
#'  obj$setAdjMat("thresh")
#'  obj$adjmat +0
#'  obj$setAdjMat()      # thresholded by default
#'  obj$adjmat+0          
#'
#' ## translate the adjacency matrix into a graphNEL object to use a
#' ## function from the graph package:
#'  graph::degree(as(obj$adjmat,"graphNEL")) 
#'  Rgraphviz::plot(as(obj$adjmat,"graphNEL"))  
#'
#' ##   Initialise the score
#'  obj$score12way() ;   obj$score[2,]
#'
#' ## Initialising the adjacency matrix (continued)
#'  obj$adjmat ; obj$setAdjMat() ; obj$adjmat +0
#'  obj$setAdjMat(mat="complete") ; obj$adjmat +0
#'  obj$setAdjMat(diag(6)) ; obj$adjmat +0
#'
#' ## Set the MI cutpoint (continued)
#'  obj$micut                #  the value of the field
#'  obj$setMIcut(alpha=.01)  #  change the significance level
#'  obj$micut                #  
#'  obj$field('micut',100)   #  set micut directly
#'  obj$micut
#'
#' 
#' ## Sample categorical data with a given independence graph
#'  set.seed(.56365)
#'  p=5
#'  amat <- jsimGraph(p,edgeprob=.5)  ; amat+0  # Here take a random graph
#'  xdata <- jsimCount(size=1000,amat)
#'  nrow(xdata) ;  head(xdata)        # look at sample of data
#'  obj <- AdjMat$new(xdata)
#'  obj$setMIcut()    
#'  obj$setAdjMat() ; obj$adjmat +0  # graph from thresholded MIs
#'  obj$mimat
#'
#' ## Score calculations
#'  obj$score12way()  ; obj$score
#'     p1hat = apply(xdata,2,mean)
#'     -sapply(p1hat,xlogx) -sapply(1-p1hat,xlogx)
#'     obj$h1vec 
#'     c(sum(obj$h1vec), obj$score[1])
#'
#'   sum(obj$adjmat * obj$mimat)/2 #  marginal MIs if in graph
#'   c(obj$nat2mbit*(1/2*log(obj$n)/obj$n),obj$bic)  # BIC penalty per  term
#'   sum(obj$adjmat)/2        #'   number of 2way interactions
#'   c(-sum(obj$adjmat*obj$mimat)/2, obj$score[2,1])
#'   c( obj$bic*sum(obj$adjmat)/2, obj$score[2,2])
#'
#'
#'# Notes on search ---------------
#' ##   Initialise the edges to visit to reduce the complexity of the graph
#'  obj$setMIcut() ; obj$setAdjMat() 
#'  obj$initEtovisit()
#'  obj$etovisit   # there are two dummy rows at the end 
#'  # Select the weakest edge and eliminate from the list 
#'  obj$nextEdge() ; obj$etovisit ; obj$nextEdge() ; obj$etovisit
#'
#' ## Find the set of potential conditioning nodes for an edge.
#'   obj$getCset(c(3,4))
#'
#' ##  Print the edge-conditioning sets (the triple) in the loop
#'  obj$initEtovisit()
#'  obj$testEdgeLoop() # the first column indexes the edge, the next 3
#'                    # list the triple, and 0 should follow
#'----------

####-----------------------------------------------
AdjMat <- setRefClass("AdjMat",
    
    contains = list("InitData"),  # requires mimat but no hlis  etc
                      
#-----------------------------------------------
    fields = list(
        adjmat = "matrix",
        # initadjmat = "character",
        etovisit = "matrix",
        score = "matrix",  # intorderx2 of sum, bic
        bic = "numeric",
        micut = "numeric"
    ),
#-----------------------------------------------
        
#-----------------------------------------------
    methods = list(
        
#-----------------------------------------------
      
#-----------------------------------------------
        setMIcut = function(alpha=.05, df=1, gamma=0.5){
            "
Sets the common cut point for all CMI tests.
The default threshold is the 0.95 quantile of
the chi-squared 1df distribution, divided by sample size and 
expressed in millibits,
The cut point can be scaled, or set to a given value.

Global argument required: n the number of rows of the data matrix,
(number of independent observations).
The optional arguments are:  
alpha, significance level;
df, degree of freedom;
gamma, optional scaling of the threshold.

Globals modified: micut, bic, numeric value in millibits.
"
             bic <<-  nat2mbit * (1/2 * log(n)/n)
             micut <<- qchisq(1-alpha,df=1)*nat2mbit/(2*n)
             micut <<- micut*exp(3*(gamma-.5))
         }
,#end setMIcut
#-----------------------------------------------

#-----------------------------------------------
        mbits2chisq = function(n, ent){" 
 input:  n integer sample size
         ent numeric,  entropy information divergence deviance in millibits
 output: numeric, chisq equivalent
"
#' mbits2chisq(n=40, ent=300)    # corresponds to obs chisq = 16.25
            return(2*n*ent/nat2mbit)
        }
,# end mbits2chisq
#-----------------------------------------------
      
#-----------------------------------------------
        setAdjMat = function(mat="thresh"){"
Initialises the boolean adjacency matrix to various matrices
determined by its options.
Options are specified by the mat argument, currently taking values
\n
 thresh: by thresholding (the first pass graph of the PC algorithm);
\n
 random: random edges;
\n
 indep: the completely independent graph;
\n
 complete: the complete graph;
\n
 A: a given pxp matrix.
"
            if(is.null(dim(mat))){
                switch(mat,
                       thresh={adjmat <<- (mimat > micut)},
                       random={adjmat <<- randomEdge(p)},
                       indep={adjmat <<- matrix(rep(FALSE,p*p),nrow=p)},
                       complete={adjmat <<- matrix(rep(TRUE,p*p),nrow=p)}
                       )
                diag(adjmat) <<- FALSE
            }else{
                if(is.matrix(mat)){
                    mode(mat) = "logical"
                    adjmat <<- mat 
                } 
            }
            dimnames(adjmat) <<- dimnames(mimat) 
        }#end setAdjMat
,
#-----------------------------------------------
        

#-----------------------------------------------
        initEtovisit = function(){"
Candidate edges.
Contains the edge and its strength (MI) for each edge in the current graph.
Globals called: adjmat, mimat.
Globals modified: etovisit.
"
            out = cbind(jallPairs(p), (mimat*adjmat)[lower.tri(mimat)] ) 
            out = out[out[,3]>0,]       # eliminates pairs with zero MI
            out = rbind(out, rep(0,3))  # two dummy rows to stop crashing 
            etovisit <<- rbind(out, rep(0,3))
        }#end initEtovisit
, 
#-----------------------------------------------

#-----------------------------------------------
        nextEdge = function(){"
The method returns the weakest edge in estrength, and eliminates it from
etovisit
Ties are handled by the R function {which.max}.
Globals called: etovisit.
Globals modified: etovisit.
Returns: wedge, the vector of nodes, miedge, associated strength.
"
            ix = which.min(etovisit[etovisit[,1]>0,3])
            # index of weakest edge in etovisit  excluding last 2 dummy rows
            wedge = etovisit[ix,1:2]
            miedge =  etovisit[ix,3]
            etovisit <<- etovisit[-ix,]  # eliminates edge
            return( list(wedge=wedge, miedge=miedge) )
        }
, #end nextEdge
#-----------------------------------------------


#-----------------------------------------------
        getCset = function(edge){"
Takes an edge and finds the cset as the  common neighbours of the edge.
If there are at least 2 such nodes then the nodes are reported
in the order of  the relative strengths of each node.
The code here requires the symmetry of mimat;
no checks is made for an empty set.
Globals called: adjmat, mimat, etovisit.
Returns the cset in ranks order of strongest first. 
"
            cset = which(adjmat[edge[1],]&adjmat[edge[2],])
            if (length(cset)<2){return(cset)}
            s = unlist(lapply(cset,
                     function(item){mimat[edge[1],item]+mimat[edge[2],item]}))
            return(cset[order(-s)])
        }#end getCset
, 
#-----------------------------------------------

#-----------------------------------------------
        removeEdge = function(pair){"
One edge in the graph is removed by updating 
{adjmat}.
"
            adjmat[pair[1],pair[2]]  <<- FALSE
            adjmat[pair[2],pair[1]]  <<- FALSE
        }#end  removeEdge
, 
#-----------------------------------------------

#-----------------------------------------------
        addEdge = function(pair){"
One edge in the graph is add by updating 
{adjmat}.
"
            adjmat[pair[1],pair[2]]  <<- TRUE
            adjmat[pair[2],pair[1]]  <<- TRUE
        }#end  addEdge
, 
#-----------------------------------------------

#-----------------------------------------------
        score12way = function(){
            "
Score set to sum of pairwise MI for all edges in adjmat.
Perhaps we should really take MIs from dlis.
BIC penalty: Used bic recalling that we work in avg entropy, to
account for the divisor of n.
We  have not adjusted for triples here - but  do so in score3way
Globals called: mimat, adjmat.
Globals modified: score.
"
#'  score=obj$score            
#'  h1vec=obj$h1vec  ; adjmat=obj$adjmat  ; mimat=obj$mimat ; bic=obj$bic
            score <<- matrix(0,nrow=3,ncol=3) # <<- 
            score[1,1] <<- sum(h1vec)
            score[1,2] <<- p
            score[1,3] <<- bic
            score[2,1] <<- -sum(adjmat * mimat)/2  # neg dij
            score[2,2] <<-  sum(adjmat)/2
            score[2,3] <<- bic
            dimnames(score) <<- list(c("1way","2way","3way"),
                                     c("entrop","df","bic"))
        }#end score12way
,
#-----------------------------------------------

#-----------------------------------------------
        testEdgeLoop = function(){"
Test loop with no access to entropy or mention of updating, to check
methods initEtovisit and nextEdge work well.
Globals called: mimat, etovisit.
Globals modified: none.
"
#' obj$testEdgeLoop()
            initEtovisit()  # 
            nedges = nrow(etovisit)-2  # two dummy 0 rows to maintain matrix
            for (e in 1:nedges){ # loop
                out =  nextEdge()  # check etovisit not empty - still a matrix
                pair = out$wedge
                mip =  out$miedge
                cset = getCset(pair) 
                for (c in cset){ # loop on elements of cset
                    mip.cset = mimat[pair[1],pair[2]]-mip  # should be 0
                    cat(c(e," ",pair," ", c," ",  mip.cset,"\n"))
                }
            }
        }#end testEdgeLoop
,
#-----------------------------------------------

##-----------------------------------------------------
develop = function(){""
}# end develop
##-----------------------------------------------------

    ))# end AdjMat
#####-----------------------------------------------

