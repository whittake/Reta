// ############################################################################
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
// 
// Author: Joe Whittaker 
// ############################################################################
#ifndef _FE3TEST_H
#define _FE3TEST_H

#include <Rcpp.h>
using namespace Rcpp;
#include <RcppEigen.h>
using namespace RcppEigen;

RcppExport SEXP jallcombCR(SEXP, SEXP) ;
// IntegerMatrix jallcombCC(IntegerVector, int, int) ;
RcppExport SEXP jfindimaxCR(SEXP, SEXP) ;
RcppExport SEXP jupdateAmatD3CR(SEXP, SEXP, SEXP, SEXP, SEXP) ;//for R version
RcppExport SEXP jupdateAmatD3C(SEXP, SEXP, SEXP, SEXP, SEXP) ;//for gFinder
RcppExport SEXP jh2deltaCR(SEXP, SEXP, SEXP, SEXP) ;
// double jh2deltaCC(IntegerMatrix, NumericVector, IntegerVector, double) ;
RcppExport SEXP jwhichSubsetCR(SEXP, SEXP) ;
 int jwhichSubsetCC(IntegerMatrix, IntegerVector) ;
RcppExport SEXP jlogdetCR(SEXP) ;
// double jlogdetc(NumericMatrix) ;

// tochuck
// RcppExport SEXP getDelta3c(SEXP, SEXP, SEXP, SEXP) ;
// RcppExport SEXP jedgeordC(SEXP, SEXP, SEXP) ;  // just for backwd compat


#endif
