// ############################################################################
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
// 
// Author: Sylvain Gubian, PMP SA
// ############################################################################
#ifndef NODE_H
#define NODE_H

#include <map>
#include <vector>
using namespace std;

class Node
{
    private:
        unsigned int _index;
        bool _visited;
        map<unsigned int, Node*> _connections;
    public:
        Node(unsigned int index) { _index = index; }
        void connectTo(Node*);
        unsigned int getIndex() { return _index; }
        void disconnect(Node*);
        bool visited() {return _visited;}
        void setVisited(bool v) { _visited = v;}
        unsigned int degree() {return _connections.size();}
        void getSubset(vector<int>&);
};

typedef map<unsigned int, Node*>::iterator NodeMapIt;
typedef map<unsigned int, Node*> NodeMap;

#endif
