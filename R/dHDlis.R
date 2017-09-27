###############################################################################
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.

# Author: Joe Whittaker 
###############################################################################

#####-----------------------------------------------

#' @title HDlis container class.
#'
#' @include cAdjMat.R
#'
#' @description
#' 
#' The intermediate computations required for graphical model search
#' are held in the fields of this class.
#' 
#' @details
#'
#' The values are held in a vector list, itself composed of several lists,
#' indexed by a tuple (ordered subset of nodes) or by a tuple-pair.
#' The lists are arranged so that tuples of the same order (size)
#' are in the same list.
#' The tuple is hashed using \code{digest}.
#' 
#' The content at each tuple of the \code{hlis} container of a given
#' order is a list of the tuple and its entropy.
#' 
#' The content at each tuple of the \code{dlis[[2]]} container is a
#' list of the tuple and its 2-way entropy interaction.
#' The content at each tuple or tuple-pair of the \code{dlis[[3]]}
#' container is a list of the tuple and 3-way entropy interaction;
#' or a list of the  tuple-pair  and 2-way entropy interaction conditional on
#' the other node.
#'
#' An  element of a typical  list is  formed from
#' \code{lis[[ltuple]][[digest(tuple)]] <<- list(ss=tuple,h=h)},
#' where \code{ltuple} is the
#' length of the \code{tuple}. 
#' 
#' The methods that initialise the list, reset the list or elements of it,
#' and access the contents of the list, are included here.
#' 
#' @field  hlis  list of lists of tuples, entropies,
#' @field  dlis  list of lists of tuple, tuple-pairs, and their interactions,
#'               and conditional interactions.
#' 
##----------------------------------------------
#' @examples
#'
#'  library(reta)
#'  #  source("retaPkg/R/dHDlis.R")
#'
#'  set.seed(1234)
#'  xdata <- data.frame(jsimBinDG(n=400))
#'  obj <- HDlis$new(xdata)
#'  # xdata <- jsimAR(n=400,p=4,alpha=.4)
#'  # obj <- HDlis$new(data)
#'  str(obj)   
#'  obj$mimat ;  obj$h12vec   # available from AdjMat 
#'  obj$initHDlis()
#'  str(obj)   # kmax hlis dlis exist but uninteresting 
#'
#' ## Initialise MI matrix and storage 
#'  obj$getallH12()  
#'  str(obj)   # kmax hlis dlis exist but uninteresting 
#'
#'  obj$xlevels     # just for catdat
#'  obj$hlis[[1]]  # non-zero and positive
#'  obj$hlis[[2]]  # positive 
#'  obj$dlis[[2]]  # negative
#'  obj$listHD(kord=2)  # all negative
#'  obj$listHD(opt="H",kord=2)   #  extracted entropies of tuples of size 2
#'  obj$listHD(opt="D",kord=2)   #  extracted 2-way interactions (-MIs)
#'  obj$listHD(opt="D",kord=2,sorted="byvalue")  
#'
#' ## component parts of the list  hlis:
#'  vector("list",3)
#'  list(ss=c(4,6),h=3.4)
#'  digest(as.integer(c(4,6)))    # unique identifier of the tuple
#'
#' ## voiding parts of hlis:
#'  HDlis$methods()
#'  obj$hlis[[2]] ; obj$resetHDlis(2)
#'  obj$hlis[[2]] 
#'  obj$field("hlis",vector("list")) ; # kills ALL stored infor
#'  obj$hlis
#'
#' ## relating HDlis entries of order 2 with  mimat
#'   set.seed(12349)
#'   xdata <- jsimAR(n=400,p=4,alpha=.4)
#'   obj <- Updater$new(xdata) 
#'   obj$mimat   # no labels as  xdata nnot a data frame
#'   obj$initHDlis() ;  obj$getallH12()
#'  obj$listHD(kord=2) ; obj$nodenms
#'   mi <- 0*diag(obj$p)
#'   mi[lower.tri(mi, diag=FALSE)] <- -obj$listHD(opt='H',kord=2)[,3]
#'   mi <- mi + t(mi) ; mi  
#'
#' ## Categorical data check that node names reproduced in lists 
#'  library(reta)
#'   xdata <- testCatDat()
#'   obj <- Updater$new(xdata) 
#'   obj$mimat # is labelled 
#'   obj$initHDlis() ;  obj$getallH12()
#'   obj$listHD(kord=2)  # is labelled 
#'   obj$listHD(kord=2,sorted='byvalue')
#'   obj$listHD(kord=2,sorted='bytuple')
#'   obj$listHD(kord=2,opt='H')
#'   obj$listHD(kord=2,opt='H',sorted='byvalue')
#'   obj$listHD(kord=2,opt='H',sorted='bytuple')
#'   obj$listHD()
#'   obj$setMIcut() ; obj$setAdjMat()      
#'   obj$adjmat   # is labelled 
#'
#'
#'
#####----------------------------------------------
HDlis <- setRefClass("HDlis",
    contains = list("AdjMat"),  
    fields = list(
        hlis = "list",
        dlis = "list",
        kmax = "integer"
    ),
##----------------------------------------------
    methods = list(
        
#####----------------------------------------------
        initHDlis = function(){"
Sets up the data base to store the entropies (hlis),
and their interactions (dlis).
The lists are subdivided by order 1,..,kmax.

NOT here  The 2nd-order entropies  (the marginal MIs) are computed from the MI matrix.

Globals required:   kmax, p, mimat.
Globals modified:   hlis,  dlis.
"
            kmax <<- as.integer(6)
            hlis <<- vector("list") 
            dlis <<- vector("list") 
            null = lapply(1:kmax,function(i){
                              hlis[[i]] <<- vector("list")
                              dlis[[i]] <<- vector("list")
                          })
        }#end initHDlis
,
##----------------------------------------------

##----------------------------------------------
        getallH12 =  function(){"
requires jallPairs mimat
different if type=cont or type=cat
if type=cont want to capitalise on xcor
one dimensional entropies - only need to compute if type=cat
"
#'  mimat=obj$mimat ; h1vec=obj$h1vec;  h12vec=obj$h12vec ; d12vec=obj$d12vec  
#'  p=obj$p ;  hlis=obj$hlis ;  dlis=obj$dlis ; jallPairs=obj$jallPairs
            out = lapply(1:p,function(item){ # computation with list output
                tuple = as.integer(item)
                hash = digest(tuple)
                hlis[[1]][[hash]] <<- list(ss=tuple,h=h1vec[item])
            })
            ## two dimensional entropies 
            pairs = jallPairs(p)
            out = lapply(1:nrow(pairs),function(item){
                                        # item=2 comp to output list
                tuple = as.integer(pairs[item,])  # 2 dim
                hash = digest(tuple)
                dlis[[2]][[hash]] <<- list(ss=tuple,d=d12vec[item])  
                hlis[[2]][[hash]] <<- list(ss=tuple,h=h12vec[item])
            })
        }#end getallH12
,
##----------------------------------------------

##----------------------------------------------
       resetHDlis = function(kord=3:kmax){"
The method removes sub hid lists  of given order(s).

Arguments: kord integer, or range of integers;
Globals fields modified: hlis,  dlis.
"
            lapply(kord,function(i){ hlis[[i]] <<- vector("list")}) 
            lapply(kord,function(i){ dlis[[i]] <<- vector("list")}) 
       }#end resetHDlis
,
##----------------------------------------------


#####----------------------------------------------
        listHD = function(opt="D",kord=3,sorted="none"){
            " 
## disabled the option ,cond=FALSE
The method returns the entropies or the entropy interactions
of a given order in matrix format.
The entropies are indexed by the tuple of nodes in the first 2 or 3 columns.
The interactions are either marginal or conditional.
## If marginal, flagged by 0 in the final column, the tuple is given first.
## If conditional, flagged by 1 in the final column, the interaction
## relates to the pair in the first 2 columns conditional on the node in the
## 3rd column.

The list may be sorted by value or tuple or not at all.

Arguments:
  opt, with values 'H'  or 'D' indicating entropy or interaction;
  kord, integer the length of the tuple plus conditioning set;
  sorted, character with values 'none', 'byvalue', 'bytuple'.
##  cond, logical indicating if just the marginal interaction is returned.

Returns:
dataframe with columns corresponding to the nodes of the tuple the value;
##  and a flag for a marginal or conditional interaction if cond=TRUE,
##  while just the marginal interaction if cond=FALSE.

Globals modified: nodenms .
"
#' dlis=obj$dlis ; kord=2 ; opt="D" # sorted="byvalue"
#' nodenms= obj$nodenms
#'  p = obj$p ; out
#'  obj$listHD(opt="D",kord=2)   #  extracted 2-way interactions (-MIs)
            ## set up out as matrix 
            if (opt=="H"){
                out = hlis[[kord]]
            }else{
                out = dlis[[kord]]
            }
            if (length(out)==0) {return("List is empty.")}
            names(out) = NULL
            # rownms: must run this before next apply
            if(length(nodenms)==0){ nodenms <<- list(1:p)}
            rownms = lapply(out,function(item){
              noquote(paste(nodenms[[1]][item[[1]]],sep='',collapse='-')) 
            })
            out = lapply(out,function(item){
                c(item[[1]],item[[2]])
            })
            out = t(matrix(unlist(out),nrow=kord+1))
            rownames(out) = rownms
            ;
            ## modify out to satisfy options
            switch(sorted,
                   none = {},
                   byvalue = {
                           out = out[order(out[,ncol(out)]),] 
                   },
                   bytuple = {
                       out[  order(rownames(out)), ]
                   }
                  )
            out = data.frame(out) #  head(out)
            switch(opt,  # adjust header 
                   H = {
                           colnms = c(paste("x",1:kord,sep=''),"H")
                   },
                   D = {
                           colnms = c(paste("x",1:kord,sep=''),"D")
                   }
                   )# switch
            colnames(out) = colnms          # head(out)
            return(out)           
        }#end listHD
,  
##----------------------------------------------


#####----------------------------------------------
plotEcdf  = function(order="2way"){""
    # todo: refresh plot 
      if(order=="2way") {
          d2way <- listHD(kord=2,opt='D')     # extract  2-way interactions
          mbits <- -d2way[,3]                     #change sign for MIs
          main = 'Pairwise MIs: ecdf'
      }else{
          d3way <- listHD(kord=3,opt='D')     # extract 3-way interactions
          if (is.data.frame(d3way)){
              mbits <- d3way[,4]
              main = '3way interactions: ecdf'
          }else{
              return("List is empty.")
          }
      }
    Fn <- ecdf(  mbits)                     # empirical distribution
    # dev.cur()  # plot.new()  ?device # to reinitialise, do but not here
      plot.stepfun(Fn,pch='.',cex=4,col="blue",
           main=main,xlab='millibits',ylab='ecdf')
      abline(h=0,v=0,col="gray60")
      abline(h=1,v=0,col="gray60")
      grid(17)
}#end plotEcdf
,  
##----------------------------------------------


#####----------------------------------------------
        listSynergies = function(deltacut=0,kord=3){"
The method identifies the synergies, the negative 3rd-order differences.  To extract marginal differences kord=3; for conditioning on sets of size one, kord is incremented by one, and so on.  The collider is determined by the smallest correlation (in absolute value).  It should be replaced by conditional (partial) correlations, when kord exceeds 0.  (and eventually by MIs)

Arguments:
 deltacut, numeric cutpoint for deciding the magnitude of the interaction;
 kord, integer, specifying the order of the difference.

cor????

Returns:
 data frame of triples, collider, synergies, and correlation coefficients.

No global fields modified.
deltacut=20 ; kord=3
"
            negdelta = listHD(opt="D",kord=3,cond=FALSE)
            ind = (negdelta[,4] < deltacut)
            sind = sum(ind)
            if (sind==0){ return("There are no synergies.")  }
            negdelta = negdelta[ind,]
            syn = t( apply(negdelta[,1:3],1,function(row){  # just the triple
                # row=c(4,6,7)
                lt = lower.tri(xcor[row,row])
                corrs = xcor[row,row][lt]
                j = which.min(abs(corrs))
                collider = row[c(3:1)[j]]  # back to vble numbering 
                c(row,collider,corrs) # mimat[row,row][lt])
            } ))
            syn = data.frame(syn)
            syn = cbind(syn[,1:4],negdelta[,4],syn[,5:7])
            nms = c("t1","t2","t3","coll","delta","cor12","cor13","cor23")
            colnames(syn) = nms
            return(syn)   
        }#end listSynergies
, 
##----------------------------------------------

#####----------------------------------------------
listTailsD3 = function(topx=2,kord=3){
            "
The method  identifies the topx largest and smallest 3rd-order interactions.

Arguments: topx, integer number.

Returns: data frame of triples, and differences.
"
            delta = obj$listHD(opt="D",kord, sorted="byvalue",cond=FALSE)
            if (is.null(delta)) {return(NULL)}
            if (nrow(delta)<4) {return(NULL)}
            ind = 1:nrow(delta)
            delta = data.frame(delta[c(head(ind,topx),tail(ind,topx)),])
            nms = c("t1","t2","t3","delta")
            colnames(delta) = nms
            return(delta)
        }#end listTailsD3
,
##----------------------------------------------

##----------------------------------------------------
develop = function(){""
}# end develop
##----------------------------------------------------

) )#end HDlis

######----------------------------------------------
