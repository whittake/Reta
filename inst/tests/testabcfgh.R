#-----------------------------------
#'  detach("package:triples", unload=TRUE)
#'  require(testthat)
#'   require(triples) #   hadley  testthat http://r-pkgs.had.co.nz/tests.html
#'   source("/home/whittake/research/13fegm/15triples/triplesPkg/inst/tests/testabcfgh.R")
#'   test_dir("triplesPkg/inst/tests/")
#'   expect_that(5 * 2, equals(11))
#' str(obj)
#-----------------------------------

test_that("CorMat RC", {
  cat("\n from aCorMat.R ")
  set.seed(1234)
  xdata <- jsimAR(n=100,p=6,alpha=.4)
  obj <- CorMat$new()  # str(obj)
  obj$readContDat(xdata) 
  expect_that( obj$p==6,  is_true())
  expect_that( dim(obj$xcor)[1]==6,  is_true())
  out = obj$xcor2mimat()
  expect_that( obj$nat2mbit>1476,  is_true())
  expect_that( dim(out$mimat)[1]==6,  is_true())
  tuple = c(1,6) 
  expect_that( obj$getHcor(tuple) < -0.5,  is_true())
  # suff stat
  obj <- CorMat$new()
  obj$readVarMat(list(n=14,S=diag(1:4)))
  obj$xcor2mimat()
  expect_that( obj$n==14,  is_true())
  expect_that( obj$p==4,  is_true())
})

test_that("CatDat RC", {
  cat("\n from aCatDat.R ")
  set.seed(1234)
  xdata <- data.frame(jsimBinDG(n=4000))
})

test_that("InitData", {
  cat("\n from bInitData.R ")
  set.seed(1234)
  # cont
  xdata <- jsimAR(n=400,p=6,alpha=.4)
  obj <- InitData$new(xdata)
  # names(InitData$fields())
  expect_that( obj$getH(c(1,3))>-25,  is_true())
  expect_that( dim(obj$xcor)[1]==6,  is_true())
  expect_that( max(obj$h12vec)> -1,  is_true())
  # cat
  set.seed(1234)
  xdata <- data.frame(jsimBinDG(n=4000))
  obj <- InitData$new(xdata)
  expect_that( obj$nat2mbit>1476,  is_true())
  expect_that( obj$p==6,  is_true())

  expect_that( round(obj$getH(c(1,3)))==1748,  is_true())
  # expect_that( obj$xlogx(4)==8192,  is_true())
  expect_that( dim(head(obj$xdf))[1]==6 ,  is_true())     
  # obj$xlevels
  expect_that( round(obj$h1vec[1])==1024,  is_true())
  expect_that( round(obj$h12vec[1])==2014,  is_true())     
})

test_that("AdjMat", {
  cat("\n from cAdjMat.R ")
  set.seed(1234)
  xdata <- jsimAR(n=400,p=6,alpha=.4)
  obj <- AdjMat$new(xdata)
  expect_that(  round(obj$h12vec[1])== -125,  is_true()) 
  obj$setMIcut()
  expect_that( round(obj$micut)==7,  is_true()) 
  obj$setAdjMat()     
  expect_that( sum(obj$adjmat[c(1,2),c(5,6)])==0,   is_true())      
  obj$score12way()
  expect_that( round(obj$score[1,3])==11,  is_true())
  obj$initEtovisit()
  expect_that(  sum(dim( obj$etovisit )-c(11,3)==0)==2,   is_true()) 
  expect_that( sum(obj$nextEdge()$wedge-c(1,3))==0,   is_true())
  expect_that( sum(obj$getCset(c(3,4))-c(2,5))==0,   is_true())
})

test_that("HDlis", {
  cat("\n from fHDlis.R ")
  set.seed(1234)
  xdata <- jsimAR(n=400,p=4,alpha=.4)
  obj <- HDlis$new(xdata)
  obj$h12vec   
  obj$initHDlis()
  obj$getallH12()  
  expect_that( length(obj$hlis[[2]])==6,  is_true())
  expect_that( length(obj$dlis[[2]])==6,  is_true())
  expect_that( round(obj$hlis[[2]][[6]]$h) == -146,  is_true())
  expect_that( round(obj$dlis[[2]][[6]]$d) == -146,  is_true())
  #
  set.seed(1234)
  xdata <- data.frame(jsimBinDG(n=400))
  obj <- HDlis$new(xdata)
  obj$initHDlis()
  obj$getallH12()  
  #  obj$xlevels     # TODO - just for catdat
  expect_that( round( obj$hlis[[1]][[5]]$h ) ==  1153,  is_true())
  expect_that( round(obj$hlis[[2]][[5]]$h) ==  1257,  is_true())
  expect_that( round( obj$dlis[[2]][[5]]$d ) == -49,  is_true())
  # hash problem
  pair =  c(4,5)
  expect_that( !(digest(pair)==digest(as.integer(pair))),  is_true())
})



test_that("Triple", {
   cat("\n first from gTriple.R ")
  set.seed(1234)
  xdata <- jsimAR(n=400,p=4,alpha=.4)
  obj <- Triple$new(xdata)
  obj$initHDlis() ; obj$getallH12()
  triple = c(1,2,3)
  expect_that( round( obj$getH(triple) )==-281,  is_true())
  expect_that( obj$getHfly(triple)-obj$getH(triple)==0,  is_true())
  expect_that( length(obj$hlis[[length(triple)]])==1,  is_true())
  expect_that( round( obj$getDelta(triple) )==24,  is_true())
  expect_that( length(obj$dlis[[3]])==1,  is_true())
  expect_that( round( obj$hlis[[3]][[1]]$h  )==-281,  is_true())
  out = class( obj$listHD(opt="H",kord=2) )
  expect_that(match(out,"data.frame")==1,is_true())

  set.seed(1234)
  xdata <- data.frame(jsimBinDG(n=400))
  obj <- Triple$new(xdata)
  obj$initHDlis() ; obj$getallH12()
  triple = c(1,2,3)
  expect_that( round(obj$getH(triple)) ==2777 ,  is_true())
  expect_that( obj$getHfly(triple)-obj$getH(triple) ==0,  is_true())
  expect_that( length(obj$hlis[[length(triple)]]) ==1,  is_true())
  expect_that( round(obj$getDelta(triple)) == -3,  is_true())
  expect_that( length(obj$dlis[[3]]) ==1,  is_true())
  expect_that( round(obj$dlis[[3]][[1]]$d) == -3,  is_true())
  out = class(obj$listHD(opt="H",kord=2))
  expect_that(match(out,"data.frame")==1,is_true())
  # inspect hlis ilis
  pair =  as.integer(c(1,2))        # vital
  hashp = digest(pair)  #  names(obj$hlis[[2]])
  expect_that( round( obj$hlis[[2]][[hashp]]$h) == 1999,  is_true())
  # hashp.A = digest(list(pair, as.integer(3)))  # str(obj)
  # expect_that(round( obj$dlis[[3]][[hashp.A]]$d )== -192,  is_true())
 })

test_that("Updater", {
 cat("\n from hUpdater.R ")
  set.seed(1234)
  xdata <- jsimAR(n=400,p=6,alpha=.4)
  obj <-  Updater$new(xdata) 
  obj$setMIcut()  
  obj$setAdjMat() 
  obj$score12way()  
  expect_that( round( obj$score[1,3] ) == 11,  is_true())
  obj$initHDlis()
  obj$getallH12() 
})

test_that("Updater", {
 cat("\n  All done.  ")
})
