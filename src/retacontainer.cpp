// ############################################################################
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
// 
// Author: Sylvain Gubian, PMP SA
// ############################################################################
#include <algorithm>
#include <assert.h>
#include "retacontainer.h"
#include "retatupleset.h"

RetaContainer::RetaContainer()
{
}

RetaContainer::~RetaContainer()
{
    RetaTupleSetMapIter it;
    for(it=_sets.begin(); it!=_sets.end(); it++)
    {
        delete it->second;
    }
}

void RetaContainer::add(const RetaTuple* rt)
{
    assert(rt!=NULL);
    unsigned int order = rt->getOrder();
    RetaTupleSetMapIter it = _sets.find(order);
    if (it == _sets.end())
    {
        RetaTupleSet* rsp = new RetaTupleSet();
        rsp->add(rt);
        RetaTupleSetPair p = std::make_pair(order, rsp);
        _sets.insert(p);
    }
    else
    {
        RetaTupleSet* rsp = it->second;
        rsp->add(rt);
    }
}

float RetaContainer::getH(const RetaTuple& rt)
{
    unsigned int order = rt.getOrder();
    RetaTupleSetMapIter it = _sets.find(order);
    if (it == _sets.end())
    {
        return -1.0; // Guessing that a H negative value is impossible, so, it means error
    }
    else
    {
        RetaTupleSet* rsp = it->second;
        RetaTuple* rtp = rsp->find(rt);
        if (rtp!=NULL)
        {
            return rtp->getH();
        }
        else
        {
            return 0;
        }

    }

}

unsigned int RetaContainer::getMaxOrder()
{
    return _sets.rbegin()->first;
}

void RetaContainer::getGlobalCombinations(
        RetaTupleSet* toFill,
        unsigned int order,
        unsigned int size)
{
    assert(toFill != NULL);
    RetaTupleSetMapIter it = _sets.find(order);
    vector<int> buffer;
    if (it != _sets.end())
    {
        RetaTupleSet* rs = it->second;
        RetaTupleSetTypeIter it2;
        for(it2=rs->start(); it2 != rs->end(); it2++)
        {
            RetaTuple* current = *it2;
            const vector<int>& data = current->getData();
            std::vector<bool> v(order);
            std::fill(v.begin() + size, v.end(), true);
            do {
                for(vector<bool>::size_type i = 0; i != v.size(); i++)
                {
                    if (!v[i])
                    {
                        buffer.push_back(data[i]);
                    }
                }
                RetaTuple* tp = new RetaTuple(buffer);
                toFill->add(tp);
                buffer.clear();
            } while(next_permutation(v.begin(), v.end()));
        }
    }
}

/*
void RetaContainer::getCombinations(

        RetaTupleSet* toFill,
        RetaTuple* input,
        unsigned int size)
{
    // TO DO
}
*/

bool RetaContainer::exists(const RetaTuple& rt)
{
    unsigned int order = rt.getOrder();
    RetaTupleSetMapIter it = _sets.find(order);
    if (it == _sets.end())
    {
        return false;
    }
    else
    {
        RetaTupleSet* rsp = it->second;
        return (rsp->exists(rt));
    }
}

RetaTuple* RetaContainer::find(const RetaTuple& rt) const
{
    RetaTupleSetMapIter it = _sets.find(rt.getOrder());
    if (it == _sets.end())
    {
        throw RetaContainerException("tuple order not found in container");
    }
    RetaTupleSet* found = it->second;
    if (!found->exists(rt))
    {
        //throw RetaContainerException("tuple not found in container");
        return NULL;
    }
    return (it->second)->find(rt);
}
