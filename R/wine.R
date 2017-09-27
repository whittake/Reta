###############################################################################
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.

# Author: Joe Whittaker 
###############################################################################
####-----------------------------------------------
#' White wine from portugal.
#'
#' A dataset containing 11 measured explanatory variables for wine
#' quality on  almost 5,000 cases.
#'
#' @format A data frame with 4898 rows and 11 variables:
#' \describe{
#'   \item{fixed.acidity}{}
#'   \item{volatile.acidity}{}
#'   \item{citric.acid}{}
#'   \item{residual.sugar}{}
#'   \item{chlorides}{}
#'   \item{free.sulfur.dioxide}{}
#'   \item{total.sulfur.dioxide}{}
#'   \item{density}{}
#'   \item{pH}{}
#'   \item{sulphates}{}
#'   \item{alcohol}{}
#' }
#'
#' @source \url{http://www.vinhoverde.pt/en/}
#' 
#' @examples
#' 
#'  library(reta)
#'  ls(package:triples)
#'  str(wine)    #    'data.frame':	4898 obs. of  11 variables:
#'  obj <- DisplayG$new(wine)  
#'  obj$initUpdater() 
#'  obj$showOptions() 
#'  obj$updateLoop() ; obj$edgechange  # 13
#'  obj$field("kappa",as.integer(4))   
#'  obj$updateLoop() ; obj$edgechange  # 27
#'  obj$field("kappa",as.integer(5))
#'  obj$updateLoop() ; obj$edgechange  # made no difference to edgechange
#'  obj$showOptions()                  # left with G having 22 edges
#'  obj$adjmat+0
#'
#'  obj$adjmat2gNEL()
#'  plot(obj$gr)
#'  obj$listSynergies(deltacut=-10)        # only one collider in this example
#'  obj$renderSynergies(-30)
#'  out = obj$listHid(opt="d3",kord=3) ; head(out)
#'  plot(ecdf(out[,4])) ; grid(12); abline(v=0); abline(v=-10)  # deltas
#'
#'    obj$listTailsD3(3)
#' DisplayG$methods()
#' 
"wine"
