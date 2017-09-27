// ############################################################################
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
// 
// Author: Sylvain Gubian, PMP SA
// ############################################################################
#include <algorithm>
#include "retatuple.h"

RetaTuple::RetaTuple()
{
    _h = 0.0;
    _computed = false;
}

RetaTuple::RetaTuple(const vector<int>& v)
{
    _data = v;
    sort(_data.begin(), _data.end());
    _h = 0.0;
    _computed = false;
}

RetaTuple::RetaTuple(const RetaTuple& rt)
{
    _data = rt._data;
    _h = rt._h;
    _computed = rt._computed;
}

RetaTuple::~RetaTuple()
{
}

void RetaTuple::setData(const vector<int>& v)
{
    _data = v;
    sort(_data.begin(), _data.end());
}

const vector<int>& RetaTuple::getData() const
{
    return _data;
}

float RetaTuple::getH()
{
    return _h;
}

void RetaTuple::setH(float v)
{
    _h = v;
}

float RetaTuple::getDelta()
{
    if (getOrder() == 2)
    {
        return _h;
    }
    else {
        // This is the forward difference that will
        // be computed at the RetaContainer level
        return _delta;
    }
}

const unsigned int RetaTuple::getOrder() const
{
    return _data.size();
}


void RetaTuple::randomize(unsigned int order)
{
    vector<int> v(order);
    for(unsigned int i=0; i < order; ++i)
    {
        v[i] = rand() % 1000000 + 1; // 1 to 10000
    }
    setData(v);
}



