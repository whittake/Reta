// ############################################################################
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
// 
// Author: Sylvain Gubian, PMP SA
// ############################################################################
#ifndef RETASUBSET_H
#define RETASUBSET_H

#include <utility>
#include <functional>
#include <unordered_set>
#include <exception>
#include <string>
#include "retatuple.h"
using namespace std;
#define MAX_SUBSET_DATA_DISPLAY  10

// Overloading hasher of unordered_map with hashing the
// RetaTuple object
// This hash_combine code is using reciprocal of the golden ratio
// http://stackoverflow.com/questions/4948780/magic-number-in-boosthash-combine
typedef vector<int>::const_iterator IVIter;

struct int_hash_combine {
    size_t operator()(const RetaTuple *rtp ) const
    {
        // RetaTuple handles vector of integers
        hash<int> hasher;
        size_t seed = 0;
        vector<int> v = rtp->getData();
        for(IVIter it = v.begin(); it != v.end(); ++it)
        {
            seed ^= hasher(*it) + 0x9e3779b9 + (seed<<6) + (seed>>2);
        }
        return seed;
    }
};


struct int_vec_equal {
    bool operator()(const RetaTuple *rtp1, const RetaTuple* rtp2) const
    {
        return (rtp1->getData() == rtp2->getData());
    }
};



// Just alias the type for better readability
typedef unordered_set< RetaTuple*, int_hash_combine, int_vec_equal> RetaTupleSetType;
typedef RetaTupleSetType::const_iterator RetaTupleSetTypeIter;


class RetaTupleSetException: public exception
{
    private:
        string _message;
    public:
        RetaTupleSetException() {
            _message = "RetaTupleSetException occured";
        }
        RetaTupleSetException(const string& msg)
        {
            _message = msg;
        }
        virtual const char* what() const throw()
        {
            return _message.c_str();
        }
        ~RetaTupleSetException() throw() {};
};


class RetaTupleSet
{
    public:
        RetaTupleSet();
        virtual ~RetaTupleSet();
        friend ostream& operator<<(ostream& os, const RetaTupleSet& rs);
        void add(const RetaTuple*);
        bool exists(const RetaTuple&) ;
        void setH(RetaTuple*, float);
        void clear();
        size_t size()
        {
            return _tuples.size();
        }
        double computeCMI();
        RetaTuple* find(const RetaTuple& rt)
        {
            RetaTupleSetTypeIter it = _tuples.find((RetaTuple*)&rt);
            if (it == _tuples.end())
            {
                //throw RetaTupleSetException("tuple not found in subset");
                return NULL;
            }
            else {
                return *it;
            }
        };
        RetaTupleSetTypeIter start() {
            return (_tuples.begin());
        }
        RetaTupleSetTypeIter end() {
            return (_tuples.end());
        }

    private:
        RetaTupleSetType _tuples;
        unsigned int _order;
} ;


inline ostream& operator<<(ostream&os, const RetaTupleSet& rt)
{
    RetaTupleSetTypeIter it = rt._tuples.begin();
    os << "RetaTupleSet-O" << rt._order << "[" << rt._tuples.size() << "]:" << endl;
    unsigned int i = 0;
    for(i=0; i < MAX_SUBSET_DATA_DISPLAY && it != rt._tuples.end(); it++, ++i)
    {
        RetaTuple* tp = *it;
        os << *tp << endl;
    }
    if (i==MAX_SUBSET_DATA_DISPLAY)
    {
        os << "..." << endl;
    }
    return os;
}


#endif
