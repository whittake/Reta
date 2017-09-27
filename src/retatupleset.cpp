// ############################################################################
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
// 
// Author: Sylvain Gubian, PMP SA
// ############################################################################
#include "retatuple.h"
#include "retatupleset.h"
#include <assert.h>

RetaTupleSet::RetaTupleSet()
{
    _order = 0; // No order defined yet
}

RetaTupleSet::~RetaTupleSet()
{
    RetaTupleSetTypeIter it;
    for(it=_tuples.begin(); it!=_tuples.end(); it++)
    {
        if (*it != NULL)
        {
            delete (*it);
        }
    }
}

void RetaTupleSet::add(const RetaTuple* rtp)
{
    // Here, would be better to check that if there is no element
    // yet, set the order of the subset
    // If there is already an element, check that the inserted one
    // is about the same order.
    if (_tuples.size() == 0)
    {
        _order = rtp->getOrder();
    }
    if (rtp->getOrder() != _order)
    {
        throw RetaTupleSetException("Cannot add tuple. Order is incorrect");
    }
    _tuples.insert((RetaTuple*)rtp);
}

bool RetaTupleSet::exists(const RetaTuple& rt)
{
    RetaTuple* rtp = (RetaTuple*)&rt;
    RetaTupleSetTypeIter it = _tuples.find(rtp);
    return(it != _tuples.end());
}

double RetaTupleSet::computeCMI()
{
    // This need to be implemented. I do not remember
    // if CMI has to be computed on tuples within the same order
    // If not, this method has to be inside RetaContainer
    return 0.0;
}

void RetaTupleSet::clear()
{
    _tuples.clear();
    _order = 0;
}

