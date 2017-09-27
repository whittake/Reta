

//==================================================================
//  #include "/home/whittake/prog/cpp/fe3code.h"   
//include "/home/whittake/research/13fegm/14gFinder/toyRCpkg/src/fe3code.h"
#include "fe3code.h"
//==================================================================

//==================================================================
//// jfindimaxCR
/* 
jfindimaxCR: finds node with largest degree.
See man pages for jfindimax.
*/ 

// [[Rcpp::export]]
SEXP jfindimaxCR(    SEXP degree_ptr, SEXP Imax_ptr    ) {
  IntegerVector deg(degree_ptr)  ;      // degree 
  IntegerVector Imax(Imax_ptr) ;        // nodes already visited 
  int p ; p = deg.size() ;                     
  // std::cout << p << std::endl ;      //  testing

  int imax ;
  if (sum(deg)<1){ return wrap(-2) ; } 
  int test=0 ;                                // indicates imax not new 
  for( int index=0; index<p ; index++){          
       imax = which_max(deg)  ; deg(imax) = 0 ;
       imax += 1 ;                            // +1 Cpp2R node
       test = min(abs(Imax-imax)) ;
       if (test>0){ return wrap(imax) ; }     // continue until new imax
  }
  return wrap(-1);                            // no new max found 
}
  
//==================================================================
// IntegerMatrix jallcombCC(IntegerVector, int) ;

//==================================================================
//// jallcombCR
/*
See manual page for jallcomb.
*/
// [[Rcpp::export]]
SEXP jallcombCR( SEXP labels_ptr, SEXP r_ptr) {
  const IntegerVector  labels(labels_ptr) ;
  int   llabels = labels.size() ; 
  int   r = as<int>(r_ptr) ;
  int   ncr = ::Rf_choose(llabels,r)  ;
  IntegerMatrix      allcomb(ncr,r) ;  // output 
  IntegerVector      comb(r) ;         // working 
  std::vector<bool>  v(llabels) ; 
  std::fill( v.begin() + llabels - r, v.end(), true) ;
  int ii, jj , kk=0 ;
  do{ jj = 0 ;                          // index of number of trues
      for (ii = 0; ii < llabels ; ++ii){  // std::cout << v[ii] ;
        if (v[ii]){ comb[jj] = labels[ii] ; ++jj ; }
      }
      for (ii = 0; ii < r ; ++ii){      // std::cout<< " " << comb[ii] << " ";
        allcomb(kk,ii) = comb[ii] ;     // [] on LHS does not work here
      }
      kk += 1 ;
  } while (std::next_permutation(v.begin(), v.end())) ;
  return wrap(allcomb) ;
}

//==================================================================
//// jupdateAmatD3CR 
/*
See jupdateAmatD3.
Updates adjacency matrix as a function, so that outputs the matrix.
*/
    
// [[Rcpp::export]]
SEXP jupdateAmatD3CR(
    SEXP Amat_ptr, SEXP MImat_ptr, SEXP micut_ptr,
    SEXP subset_ptr, SEXP delta3_ptr
    ) {
  LogicalMatrix   Amat(Amat_ptr) ;
  NumericMatrix   MImat(MImat_ptr) ;
  double micut = as<double>(micut_ptr) ;
  IntegerVector  subset(subset_ptr) ;
  double delta3 = as<double>(delta3_ptr) ;

  int npair=2,  lsubset=3 ; 
  int ii, jj ; 
  double  minew ;
  IntegerVector pair(npair) ;
  
  std::vector<bool> v(lsubset) ; 
  bool edge ;  
 
  minew = micut + delta3 ;
  std::fill(v.begin() + lsubset - npair, v.end(), true) ;
  do{
      jj = 0 ;
      for (ii = 0; ii < lsubset ; ++ii){
        if (v[ii]){ pair[jj] = subset[ii]-1 ; ++jj ;}
      }
      //---  update edge from  pair ;
      edge = (MImat(pair[0],pair[1]) > minew ) ;  
      // edge=4  testing 
      Amat(pair[0],pair[1]) = edge ;
      Amat(pair[1],pair[0]) = edge ;  
  } while (std::next_permutation(v.begin(), v.end())) ;
  return    wrap(Amat) ;  // wrap(Amat)
}//end jupdateAmatD3CR

//==================================================================
//// jupdateAmatD3C
/*
Mutated version.
See the manual page for ?jupdateAmatD3.
*/
    
// [[Rcpp::export]]
SEXP jupdateAmatD3C(
    SEXP Amat_ptr, SEXP MImat_ptr, SEXP micut_ptr,
    SEXP subset_ptr, SEXP delta3_ptr
    ) {
  LogicalMatrix   Amat(Amat_ptr) ;
  NumericMatrix   MImat(MImat_ptr) ;
  double micut = as<double>(micut_ptr) ;
  IntegerVector  subset(subset_ptr) ;
  double delta3 = as<double>(delta3_ptr) ;

  int npair=2,  lsubset=3 ; 
  int ii, jj ; 
  double  minew ;
  IntegerVector pair(npair) ;
  
  std::vector<bool> v(lsubset) ; 
  bool edge ;  
 
  minew = micut + delta3 ;
  std::fill(v.begin() + lsubset - npair, v.end(), true) ;
  do{
      jj = 0 ;
      for (ii = 0; ii < lsubset ; ++ii){
        if (v[ii]){ pair[jj] = subset[ii]-1 ; ++jj ;}
      }
      //---  update edge from  pair ;
      edge = (MImat(pair[0],pair[1]) > minew ) ;  
      // edge=4  testing 
      Amat(pair[0],pair[1]) = edge ;
      Amat(pair[1],pair[0]) = edge ;  
  } while (std::next_permutation(v.begin(), v.end())) ;
  return    wrap(-1) ;  // wrap(Amat)
}//end jupdateAmatD3C

//==================================================================
//// jh2deltaCR
/*
See man for jh2delta in fe3RtoCfns.R.

Note: Heading for non-exportable C-C function would be 
 double jh2deltaCC(IntegerMatrix Subsets2, NumericVector H2, 
                  IntegerVector trip, double h3){ 
*/ 
// [[Rcpp::export]]
SEXP jh2deltaCR( SEXP Subsets2_ptr, SEXP H2_ptr,
		SEXP trip_ptr, SEXP h3_ptr
		) {
  IntegerMatrix Subsets2(Subsets2_ptr) ;
  NumericVector H2(H2_ptr) ;
  IntegerVector trip(trip_ptr) ;
  double h3 = as<double> (h3_ptr) ;
  int npair = 2 ; int ntrip = 3 ; 
  IntegerVector pair(npair) ;
  std::vector<bool> v(ntrip) ;
  double  h2, sumh2, delta3 ; 
  int ii, jj ;
  sumh2 = 0 ;
  std::fill(v.begin() + ntrip - npair, v.end(), true) ;
  do{
      jj = 0 ;  // index of local pair
      for (ii = 0; ii < ntrip ; ++ii){
          if (v[ii]){ pair[jj] = trip[ii] ; ++jj ;}
      }     
      sumh2 += H2[jwhichSubsetCC(Subsets2,pair)-1] ; // C indexing
  } while (std::next_permutation(v.begin(), v.end())) ;
  delta3 = h3 - sumh2 ;
  return wrap(delta3)   ;  // delta3
}

//==================================================================
//// jwhichSubsetCC
// This version is for C-C usage.
int jwhichSubsetCC(IntegerMatrix Subsets, IntegerVector subset){
    int ncol = Subsets.ncol() ;
    if (ncol != subset.size()) {return -2 ;}
    // Iterate on matrix rows
    for(int i=0; i < Subsets.nrow(); i++) {
        if (sum(Subsets.row(i)==subset) == ncol){ return i+1 ;}
        //+1 for R indexes, in C++ context, return i 
    }
    return -1 ;
}
//==================================================================
//// jwhichSubsetCR
/*
This version is called by R, see the help for jwhichSubset
in fe3RtoCfns.R.

@param Inputs: Subsets, subset.   eg Subsets3, trip
@return Integer: Rindex of the subset in Subsets if subset present.
Error messages: 
   -2 if Subsets and subset different sizes
   -1 if subset not present.
@examples # jwhichSubsetCR(Subsets3, trip)
*/
// [[Rcpp::export]]
SEXP jwhichSubsetCR( SEXP Subsets_ptr, SEXP subset_ptr){
  IntegerMatrix Subsets(Subsets_ptr) ;
  IntegerVector subset(subset_ptr) ;
    int ncol = Subsets.ncol() ;
    if (ncol != subset.size()) {return wrap(-2) ;}
    // loop on matrix rows
    for(int i=0; i < Subsets.nrow(); i++) {
        if (sum(Subsets.row(i)==subset) == ncol){ return wrap(i+1) ;}
        //+1 for R indexes, in C++ context, return i 
    }
    return wrap(-1) ;
}

//==================================================================
//// jlogdetCR
/*
See man for jlogdet in fe3RtoCfns.R.
Can internalise code with call double jlogdetCR(NumericMatrix A_ptr).
*/
// [[Rcpp::export]]
SEXP jlogdetCR(SEXP A_ptr) {
  using Eigen::MatrixXd ;
  using Eigen::VectorXd ;
  using Eigen::Map ;
  typedef Map<MatrixXd> MapMatd ;   // matrix2array jiggery
  MatrixXd  A( as<MapMatd>(A_ptr) );
  // std::cout << A << std::endl ;      // for testing
  VectorXd  Dvec ;
  Dvec = A.ldlt().vectorD() ;           // chol LDL^T decomp
  return wrap( Dvec.array().log().sum() ) ;
}

//==================================================================

