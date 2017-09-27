// ############################################################################
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
// 
// Author: Sylvain Gubian, PMP SA
// ############################################################################
#include <Rcpp.h>
#include "retatuple.h"

RCPP_MODULE(RetaTupleWrapper){
    using namespace Rcpp ;

    class_<RetaTuple>( "RetaTuple" )

    .default_constructor()

    // read and write property
    .property( "H", &RetaTuple::getH, &RetaTuple::setH )
    ;
}

