###############################################################################
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.

# Author: Joe Whittaker 
###############################################################################

####-----------------------------------------------
#' @title Updater reference class
#'
#' @include    eTriple.R
#' 
#' @description
#'
#' The \code{Updater} reference class holds the \code{thin} method
#' that executes the search through the edges and updates on the basis
#' of approximate MI tests. It holds related methods for
#' initialising, displaying and setting options.
#'
#' @details
#'
#' The ordering of  the edges and of the conditioning
#' sets is deterministic, based on measures of strength calculated
#' from the MI matrix. These are handled in the class
#' \code{\link{adjmat}}.  Exception handling for instance, nodes
#' that have insufficient neighbours, or conditioning sets that are
#' too large singleton nodes, is also described.
#'
#'
#' The method \code{thin} executes a single pass through the
#' list of edges and their associated conditioning sets.
#'
#' This reference class holds the \code{score3way} for computing the
#' 3way interaction score of a given graph The score itself is
#' composed of two terms. The first term  \code{score[2,1]} is the
#' sum of the negative marginal MI of each each in the graph calculated from
#' \code{adjmat} and \code{mimat}.
#' The second \code{score[3,1]} is the adjustment for
#' any complete triangle of edges computed from the sum of 3-way
#' entropy interaction terms held in \code{dlis[[3]]}.  The score is
#' computed afresh from the current \code{adjmat} and is not updated
#' during the search \code{thin}.
#' 
#'
#' @seealso \code{\link{CorMat}}, \code{\link{CatDat}},
#'   \code{\link{InitData}}, \code{\link{HDlis}}, \code{\link{Triple}},
#'
#' @examples
#' 
#'  getwd()
#'  require(reta)  ;  require(Rgraphviz)
#'  require(testthat) 
#'  source("retaPkg/inst/tests/testabcfgh.R")
##  source("retaPkg/R/fUpdater.R")
#'
#'
#'  require(reta) ;  require(Rgraphviz)
#'  obj <- Updater$new()
#'  str(obj)
#'  set.seed(1234)
#'  xdata <- jsimAR(n=400,p=7,alpha=.4)
#'  obj <- Updater$new(xdata)
#'  obj$n ; obj$p ; obj$mimat 
#'
#' ## search for a graph
#'  obj$setMIcut()  ;  obj$micut ; obj$bic
#'  obj$setAdjMat() ;  obj$adjmat+0 # thresholded 
#'  obj$score12way() ;  obj$score   # 1129
#'  obj$initHDlis()
#'  obj$getallH12() 
#'  obj$thin()      # could use obj$skeleton  
#'  obj$adjmat+0    # pretty much like that of the AR1 process
#'  obj$score12way()  ;  obj$score  # 1010  reduced but with fewer terms
#'
#' ## alternative initialisations
#'  obj$setAdjMat("thresh") ; obj$adjmat+0    # starting adjmat 
#'  obj$thin() ;   obj$adjmat+0            # finishing adjmat
#'  obj$setAdjMat("complete") ; obj$adjmat+0  # starting adjmat 
#'  obj$thin() ;   obj$adjmat+0            #  gets there
#'  obj$setAdjMat("random") ; obj$adjmat+0    # starting adjmat 
#'  obj$thin() ;   obj$adjmat+0            # fails to arrive!!
#'  # because wrong 0's are always maintained in thinning
#'
#' ## moralization and further search
#'  obj$setAdjMat("random") ; obj$adjmat+0    # starting adjmat
#'  obj$thin() ; obj$moralize() ; obj$thin() ; obj$adjmat+0
#' 
#' ## basicTriples combines initialisation, thinning, and moralization
#'  obj$basicTriples() ;   obj$adjmat+0
#' 
#' ## scoring, take another example 
#'  set.seed(123456)
#'  amat <- jsimGraph(p=5,edgeprob=.5)  ; amat+0   # a butterfly  graph
#'  xdata <- jsimGaussianGM(400,amat)
#'  obj <- Updater$new(xdata)
#'  obj$basicTriples() ;   obj$adjmat+0
#'  plot(as(obj$adjmat,"graphNEL"))
#'  obj$mimat
#'  obj$score12way() ; obj$score
#'  # 3way interaction and score 
#'  obj$listHD()   
#'  obj$listHD("D",kord=3,sorted="byvalue")
#'  obj$score3way() ; obj$score 
#'  # the 3way score has to be adjusted down by the two triangles in graph
#'  c(-2109 - 2017, 2*obj$bic)  #  c(- d234 - d145, 2*bic)
#'
#'
#' ## categorical data, from a graphical model with a random graph 
#'  set.seed(25784)
#'  amat <- jsimGraph(p=5,edgeprob=.6)  ; amat+0   # the sampled graph
#'  xdata <- jsimCount(size=400,amat)
#'  obj <- Updater$new(xdata)
#'  obj$basicTriples()
#'  obj$adjmat-amat               
#'  wrongEdges(obj$adjmat,amat)   # not great  # 1 FP 4 FN
#'
#'  obj$micut
#'  obj$mimat    # rather small
#'  obj$listHD("D",kord=2,sorted="byvalue")  
#'  obj$listHD("D",kord=3,sorted="bytuple")
#'  obj$listHD("D",kord=3,sorted="byvalue")
#'
#' ## example  that shows thinning after moralizing can give graph identical to
#'#  skeleton
#'  set.seed(1234)
#'  xdata <- data.frame(jsimBinDG(n=4000))
#'  obj <- Updater$new(xdata)
#'  obj$skeleton()  
#'  Rgraphviz::plot(as(obj$adjmat,"graphNEL"))
#'  obj$listHD()
#'  obj$moralize() 
#'  obj$listHD()
#'  Rgraphviz::plot(as(obj$adjmat,"graphNEL"))
#'  obj$thin() ;    
#'  Rgraphviz::plot(as(obj$adjmat,"graphNEL"))
#'  obj$basicTriples()  
#'  obj$listHD()
#'  Rgraphviz::plot(as(obj$adjmat,"graphNEL"))
#'  obj$setAdjMat("complete") ; obj$adjmat+0 ; obj$micut
#'  obj$thin() 
#'  Rgraphviz::plot(as(obj$adjmat,"graphNEL"))
#'  obj$score12way()  
#'  obj$score3way()  ;  obj$score  # no triangles
#'
#' #-------------------------------------------------------------
#'  ## examples with negative delta TODO
#'  out=obj$listHD(kord=3,cond=FALSE) ;  delta3=out[,4]
#'  plot(ecdf( delta3 )) ; grid(12); abline(v=0)
#'                                     #  separation in MI small large
#'  plot(ecdf( obj$listHD(opt="D")[,4] )) ; grid(12); abline(v=0)
#'                                     # a range of values, some negative
#'
#'   obj$listTailsD3()
#'   obj$listSynergies(deltacut=20)
#'
#' 
#' ##  match triple in two lists - see xamples - fix sort bytuple
#' ##  obj$matchTriples()  # get from archive it wanted 
#' #-------------------------------------------------------------
#'
#' @export
#'
##----------------------------------------------
Updater <- setRefClass("Updater",
     contains = list( "Triple"),
                       
     methods = list(
##----------------------------------------------

##----------------------------------------------
         basicTriples = function(alpha=0.05){
             "
"
             setup(alpha) 
             thin()
             moralize()
             thin()
             score3way()
         }#end basicTriples
,
##----------------------------------------------

##----------------------------------------------
         moralize = function(){
             "
"
             nodeLoop(option='moralize')
         }#end basicTriples
,
##----------------------------------------------

##----------------------------------------------
         score3way = function(){
             "
"
             nodeLoop(option='score')
         }#end basicTriples
,
##----------------------------------------------

##----------------------------------------------
         skeleton = function(alpha=0.05){
             "
follows the independence
edge testing routine based on approximate conditional mutual information 
to thin the given graph.
If the underlying truth is a Bayes network, this generates its skeleton.
"
             setup(alpha) 
             thin()
         }#end setup
,
##----------------------------------------------
      
##----------------------------------------------
         setup = function(alpha=0.05){
             "
uses the default settings to setup  the triples code and find the twi score
incorporates the initial pass of the PC algorithm
"
             setMIcut(alpha) 
             setAdjMat() 
             score12way() 
             initHDlis()
             getallH12() 
         }#end setup
,
##----------------------------------------------

##----------------------------------------------
         thin = function(){
             "
The loop visits each edge of the current graph finding conditioning sets, and building an approximate conditional MI test statistic from the sequence of potential covariates.     The edge when the test statistic falls below \\code{micut}.

Globals modified:
\\code{adjmat}, 
\\code{hlis,dlis}, teh data bases for the entropies, entropy interactions.
"
#'  obj$setAdjMat(complete) ; # obj$adjmat+0
#'  obj$initEtovisit
             initEtovisit()  # 
             if (nrow(etovisit)==2){return('Stop no edges to inspect')}
             nedges = nrow(etovisit)-2 # two dummy 0 rows to maintain matrix
             # cat(nedges)
             for (e in 1:nedges){   # loop e=1
                 out =  nextEdge() # check etovisit not empty - still matrix
                 pair = out$wedge  # is.integer(pair)
                 mip =  out$miedge
                 cset = getCset(pair)
                 mip.cset = mip
                 if (mip.cset<micut) {
                     removeEdge(pair)
                     # print('edgeMgone')
                     next(e)
                 }# if mi
                 for (c in cset){ # loop on elements of cset  c=3
                     triple = sort(c(pair, c))
                     delta = getDelta(triple)
                     if(delta>0){
                         mip.cset = mip.cset - delta # decrement 
  ## to store each test made, need to adapt getDelta code, hacking 
  ## hashp.k = digest(list(pair,k)) # list storage, different to triple
  ## dlis[[ltuple]][[hashp.k]] <<- list(ss=list(pair,k),d=-mip.k)  #<<-
                         #  cat(c(pair,"xxx",c,888,mip.cset,delta,"\n"))
                         #  mip.cset can go negative ie poor approx?
                         if(mip.cset<micut){
                             removeEdge(pair)
                             next(e)
                         }# if mi
                     }# if delta
                 }# for c 
             }# for e
         } # #end thin
,
##----------------------------------------------

##----------------------------------------------
    nodeLoop = function(option='score'){
        "

The  call \\code{nodeLoop()} passes through each node of the graph finding pairs and then triples=(node,pair)=(i23) say.
A check is made to see if the  triple is new to this search.
If so, the 3way entropy interaction di23 is evaluated or retrieved. 
If the score option is set the 3way score is updated by
incrementing  \\code{score[3]} by di23 and penalising by \\code{bic}.
If the moralize option is set and if  di23 is 
 negative, the edge between 2-3 is set to added to the adjacency matrix.

Requires: adjmat, HDlis via call to getDelta, micut
Globals modfied: adjmat.

"
#'  getDelta = obj$getDelta 
#'  adjmat = obj$adjmat ; p=obj$p ; score =obj$score
#'  plot(as(amat,"graphNEL"))
#'  plot(as(obj$adjmat,"graphNEL"))
        tvisited = list()
        if (option=='score'){
            score[3,1:3] <<- rep(0,3)
            score[3,3] <<- bic
        }
        for (inode in 1:p){  # inode = 4
            pairsmat = as.matrix(integer(0))
            tt = adjmat[inode,]*(1:p)
            neighbours = tt[tt>0]    # of inode # print(neighbours) 
            lnei = length(neighbours)  
            if ((lnei<2)){ next(inode) }
            pairsmat = .Call('jallcombCR',neighbours,2)
            for (ipair in 1:nrow(pairsmat)){ # ipair=1
                pair = pairsmat[ipair,]
                # print(c(pair[1],pair[2]))
                triple = as.integer(sort(c(pair, inode)))
                hasht = digest(triple)
                if(!is.null(tvisited[[hasht]])){
                    next(ipair)
                }
                tvisited[[hasht]] = TRUE  # triple
                delta3 =  getDelta(triple)
                if ((option=='score')&(sum(adjmat[triple,triple])==6)){  
                     score[3,1] <<- score[3,1] + delta3   #<<-
                     score[3,2] <<- score[3,2] +1 
                }
                if (option=='moralize'){
                    if (delta3 <  -micut){ addEdge(pair)}
                }
            }# ipair
        }# inode
    }#end  nodeLoop   
,     
##----------------------------------------------




##-----------------------------------------------------
develop = function(){""
#' ## Processing the output from the updating loop
#'  set.seed(1234)
#'  xdata <- jsimAR(n=400,p=5,alpha=.4)
#'  obj <- Updater$new(xdata)
#'  obj$basicTriples()   # the minimal statments to run a search
#'  
#'  plot(as(obj$adjmat,"graphNEL"))
#'  head( obj$listHD(kord=2) )          # edge strengths, all 2way interactions
#'  plot(ecdf( -obj$listHD(kord=2)$D ),
#'             main='ecdf') ; grid(12); abline(v=0) ; abline(h=0)
#'
#'  obj$listHD(opt="H")                 # triples with evaluated entropies
#'  obj$listHD("D",kord=3)              # all positive here 
#'  obj$getD(c(1,2,3))                  # interaction of specific tuples
#'  obj$getD(c(1,2,3,4))                # interaction of specific tuples
    #'


#'  ## Now separate the marginal and conditional entries
#'   set.seed(1234)
#'   xdata <- jsimAR(n=400,p=4,alpha=.4)
#'   obj <- Updater$new("cont",xdata) 
#'   obj$setMIcut() ; obj$micut
#'   obj$setAdjMat() ; obj$adjmat
#'   obj$initHDlis() ;  obj$getallH12()
#'   obj$thin()  # occasional error
#'   obj$listHD("D",kord=3)
#'
#'   obj$listHD("D",kord=3,cond=FALSE)   #   TODO? now appears ok
#'   obj$listHD("D",kord=3,cond=TRUE) 
#'   #  3-way marginal interaction
#'    #  2-way interaction x1-x2 conditional on x3
#'                     # ie x1\indep{}x2|x3 is +166.41
#' ## comparison of  marginal and conditional deltas
## UPTO here : not using 4wya at the moment
#'  head(  obj$listHD(opt="D",kord=3,sorted="bytuple")) 
#'  head(  obj$listHD(opt="D",kord=4,sorted="bytuple"))

    ##----------------------------------------------
     score3wayOLD = function(){
         "
The  call \\code{score3way()} passes through each node of the graph finding pairs and then triples=(node,pair)=(i23) say.
A check is made to see if the  triple is new to the search.
If so, the 3way entropy interaction di23 is checked for negativity, and 
if negative, the edge between 2-3 is set to TRUE.

Requires: adjmat, HDlis via call to getDelta, micut
Globals modfied: adjmat.

"
#'  obj$  
#'  adjmat = obj$adjmat ; p=obj$p
#'  plot(as(amat,"graphNEL"))
#'  plot(as(obj$adjmat,"graphNEL"))
         tvisited = list()
         for (inode in 1:p){  # inode = 6
            pairsmat = as.matrix(integer(0))
            tt = adjmat[inode,]*(1:p)
            neighbours = tt[tt>0]    # of inode # print(neighbours) 
            lnei = length(neighbours)  
            if ((lnei<2)){ next(inode) }
            pairsmat = .Call('jallcombCR',neighbours,2)
            for (ipair in 1:nrow(pairsmat)){ # ipair=1
                pair = pairsmat[ipair,]
                # print(c(pair[1],pair[2]))
                triple = as.integer(sort(c(pair, inode)))
                hasht = digest(triple)
                if(!is.null(tvisited[[hasht]])){
                    next(ipair)
                }
                delta3 =  getDelta(triple)
                tvisited[[hasht]] = TRUE  # triple
                if (delta3<  -micut){   # if (delta3<0){ replaced
                    adjmat[pair[1],pair[2]] <<- TRUE
                    adjmat[pair[2],pair[1]] <<- TRUE
                }# if 
            }# ipair
        }# inode
    }#end   junk  
     
##----------------------------------------------

    
}# end develop
##-----------------------------------------------------


) )# end Updater

#####--------------------------------------------------------
