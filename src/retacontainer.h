// ############################################################################
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
// 
// Author: Sylvain Gubian, PMP SA
// ############################################################################
#ifndef RETACONTAINER_H
#define RETACONTAINER_H

#include <map>
#include <utility>
#include "retatupleset.h"
using namespace std;

typedef std::map<unsigned int, RetaTupleSet*> RetaTupleSetMap;
typedef std::map<unsigned int, RetaTupleSet*>::const_iterator  RetaTupleSetMapIter;
typedef std::pair<unsigned int, RetaTupleSet*> RetaTupleSetPair ;

class RetaContainerException: public exception
{
    private:
        string _message;
    public:
        RetaContainerException() {
            _message = "RetaContainerException occured";
        }
        RetaContainerException(const string& msg)
        {
            _message = msg;
        }
        virtual const char* what() const throw()
        {
            return _message.c_str();
        }
        ~RetaContainerException() throw() {};
};

class RetaContainer
{
    private:
        // A hashmap of RetaTupleSets stored based on the order (key)
        // This is the equivalent of the RetaLis that we have implemented in R
        RetaTupleSetMap _sets;
        RetaTupleSet _buffer;

    public:
        RetaContainer();
        virtual ~RetaContainer();
        void add(const RetaTuple*);
        float getH(const RetaTuple&);
        unsigned int getMaxOrder();
        const RetaTupleSet& getSubset(unsigned int order);
        void getGlobalCombinations(RetaTupleSet*, unsigned int order, unsigned int size);
        //void getCombinations(RetaTupleSet*, RetaTuple*, unsigned int size);
        bool exists(const RetaTuple&);
        RetaTuple* find(const RetaTuple&) const;
};

#endif
