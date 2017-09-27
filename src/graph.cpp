// ############################################################################
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
// 
// Author: Sylvain Gubian, PMP SA
// ############################################################################

#include "graph.h"

Node* Graph::getNodePtr(unsigned int index)
{
    NodeMapIt it;
    it = _nodes.find(index);
    if (it != _nodes.end())
    {
        return it->second;
    }
    else {
        return NULL;
    }
}

void Graph::addNode(Node* n)
{
    pair<unsigned int, Node*> p;
    p = make_pair(n->getIndex(), n);
    _nodes.insert(p);
}

Graph::~Graph()
{
    NodeMapIt it;
    for(it = _nodes.begin(); it != _nodes.end(); it++)
    {
        if (it->second)
        {
            delete(it->second);
        }
    }
}

// Weakest tuple definition:
// If we have tuple {1,2,3,j}, the weakest tuple is the one
// that has h(i,j) + h(k,j) lowest value
// Say the weakest is {I,K,J}, the test will be:
// h(i,k) + delta(I,K,J) > cutoff for all combinations of i,k
// of {I,K,J}
// deltas will be calculated with:
//

Node* Graph::getNonVisitedMaxDegreeNode()
{
    NodeMapIt it;
    unsigned int maxDegreeNode = 0;
    Node* selected = 0;
    for(it = _nodes.begin(); it != _nodes.end(); it++)
    {
        Node* n = it->second;
        if (!n->visited())
        {
            if (n->degree() > maxDegreeNode)
            {
                maxDegreeNode = n->degree();
                selected = n;
            }
        }
    }
    return selected;
}
