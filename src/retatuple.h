// ############################################################################
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
// 
// Author: Sylvain Gubian, PMP SA
// ############################################################################
#ifndef RETATUPLE_H
#define RETATUPLE_H

#include <vector>
#include <iostream>
#pragma once
using namespace std;

#define MAX_TUPLE_DATA_DISPLAY  10

class RetaTuple
{
    private:
        vector<int> _data;
        float _h;
        float _delta;
        bool _computed;

    public:
        RetaTuple();
        RetaTuple(const vector<int>&);
        RetaTuple(const RetaTuple&);
        virtual ~RetaTuple();
        friend ostream& operator<<(ostream& os, const RetaTuple& rt);
        void setData(const vector<int>&);
        const vector<int>& getData() const;
        float getH();
        void setH(float);
        float getDelta();
        void setDelta();
        const unsigned int getOrder() const;
        bool operator==(const RetaTuple& rt) const
        {
            // Two tuples will be equal if the sorted
            // _data vectors will be equal
            // I changed the strategy, they are sorted at insert
            // so no need to sort in this operator definition
            return _data == rt._data;
        };
        RetaTuple& operator=(const RetaTuple& rt)
        {
            _data = rt._data;
            _h = rt._h;
            _computed = rt._computed;
            return *this;
        }
        // For testing purpose only. Randomizing a tuple after being stored will make it
        // impossible to find.
        void randomize(unsigned int order);
};


inline ostream& operator<<(ostream&os, const RetaTuple& rt)
{
    os << "RetaTuple-O" << rt._data.size() << " (" ;
    unsigned int i;
    for(i=0; i < MAX_TUPLE_DATA_DISPLAY && i < rt._data.size(); ++i)
    {   if (i > 0)
        {
            os << ", ";
        }
        os << rt._data[i];
    }
    if (i==MAX_TUPLE_DATA_DISPLAY)
    {
        os << ", ...";
    }
    os << ")";
    return os;
}

#endif

