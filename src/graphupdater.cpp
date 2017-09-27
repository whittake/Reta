// ############################################################################
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
// 
// Author: Sylvain Gubian, PMP SA
// ############################################################################
#include "graphupdater.h"
#include <math.h>

GraphUpdater::GraphUpdater()
{
    _cormatPtr = 0;
}

void GraphUpdater::init()
{
    if (!_cormatPtr)
    {
        return;
    }
    for(int i=0; i < _cormatPtr->nrow(); ++i)
    {
        Node* currentNode = new Node(i);
        _graph.addNode(currentNode);
        for(int j=0; j < i; ++j)
        {
            int arr[] = {i, j};
            vector<int> vec(arr, arr + sizeof(arr) / sizeof(int));
            RetaTuple* t = new RetaTuple(vec);
            // Not sure about -log
            t->setH(-log(1.0 - (*_cormatPtr)(i,j) * (*_cormatPtr)(i,j)) / 2.0);
            _container.add(t);
            Node* node = new Node(j);
            _graph.addNode(node);
            currentNode->connectTo(node);
        }
    }
}

void GraphUpdater::updateGraph(unsigned int order)
{
    /*
    NodeMapIt it;
    NodeMap* nodes = _graph.getNodes();
    Node* n = getNonVisitedMaxDegreeNode();
    do {
        vector<int> v;
        n->getSubset(v);
        RetaTuple* t = _container.getWeakestTuple(v);
        RetaTupleSet* ts = _container.getCombinations(t,t->getOrder()-1);
        _container.computeForwardDiff(ts);
        n->setVisited(true);
        n = getNonVisitedMaxDegreeNode();
    } while (!n);
    */
}

Node* GraphUpdater::getNonVisitedMaxDegreeNode()
{
    Node* n;
    n = _graph.getNonVisitedMaxDegreeNode();
    return n;
}

