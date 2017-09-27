// ############################################################################
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
// 
// Author: Sylvain Gubian, PMP SA
// ############################################################################
#ifndef GRAPH_H
#define GRAPH_H

#include <map>
#include "node.h"
using namespace std;

class Graph {
    private:
        map<unsigned int, Node*> _nodes;
    public:
        Graph() {};
        virtual ~Graph();
        Node* getNodePtr(unsigned int);
        void addNode(Node*) ;
        NodeMap* getNodes() {return &(_nodes);}
        Node* getNonVisitedMaxDegreeNode();
};

#endif
