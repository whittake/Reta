// ############################################################################
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
// 
// Author: Sylvain Gubian, PMP SA
// ############################################################################
#include "node.h"
#include <utility>
void Node::connectTo(Node* n)
{
    pair<unsigned int, Node*> p;
    NodeMapIt it = _connections.find(n->getIndex());
    if (it == _connections.end())
    {
        p = make_pair(n->getIndex(), n);
        _connections.insert(p);
        n->connectTo(this);
    }
}

void Node::disconnect(Node* n)
{
    NodeMapIt it = _connections.find(n->getIndex());
    if (it != _connections.end())
    {
        _connections.erase(it);
        Node* found = it->second;
        found->disconnect(n);
    }
}


void Node::getSubset(vector<int>& v)
{
    // From a given node 5 with neighbours {2,3,4}, it created a tuple
    // {2,3,4,5}
    // Return a vector of indexes: node + connections indexes
    NodeMapIt it;
    for (it = _connections.begin(); it != _connections.end(); ++it)
    {
        v.push_back(it->first);
    }
    v.push_back(_index);
}
