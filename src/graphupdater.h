// ############################################################################
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
// 
// Author: Sylvain Gubian, PMP SA
// ############################################################################
#include "Rcpp.h"
#include "retacontainer.h"
#include "node.h"
#include "graph.h"

using namespace Rcpp;

class GraphUpdater
{
    private :
        double _m ;
        NumericMatrix* _cormatPtr ;
        Graph _graph ;
        RetaContainer _container ;
    public :
        GraphUpdater();
        void init() ;
        void updateGraph(unsigned int order) ;
        void setCorMat(NumericMatrix* mPtr) {_cormatPtr = mPtr;}
        Node* getNonVisitedMaxDegreeNode();
};

