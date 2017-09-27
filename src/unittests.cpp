// ############################################################################
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
// 
// Author: Sylvain Gubian, PMP SA
// ############################################################################
#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "retatuple.h"
#include "retatupleset.h"
#include "retacontainer.h"
#include "graph.h"
#include "node.h"
#include "graphupdater.h"
#include <vector>
#include <iostream>
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */

using namespace std;

class RetaContainerFixture
{
    protected:
        RetaContainer* c;
        RetaTupleSet* s;
    public:
        RetaContainerFixture()
        {
            srand (time(NULL));
            c = new RetaContainer();
            s = new RetaTupleSet();
        }
        ~RetaContainerFixture()
        {
            delete s;
            delete c;
        }
};

class RetaGraphFixture
{
    protected:
        Node* node1;
        Node* node2;
        Node* node3;
        Graph* graphPtr;
    public:
        RetaGraphFixture()
        {
            node1 = new Node(1);
            node2 = new Node(2);
            node3 = new Node(3);
            graphPtr = new Graph();
        }
        ~RetaGraphFixture()
        {
            delete graphPtr;
        }
};

class RetaGraphUpdaterFixture
{
    protected:
        GraphUpdater* guptr;
    public:
        RetaGraphUpdaterFixture()
        {
            guptr = new GraphUpdater();
        }
        ~RetaGraphUpdaterFixture()
        {
            delete guptr;
        }
};

TEST_CASE_METHOD(RetaGraphUpdaterFixture, "GraphUpdater is a working component", "[node]" )
{
    SECTION("We can initialize a Graphupdater using a NumericMatrix")
    {
        NumericMatrix* mptr = new NumericMatrix(3,3);
        //NumericMatrix m(3,3);

        //guptr->setCorMat(mptr);
    }
}

TEST_CASE_METHOD(RetaGraphFixture, "Node is a working component", "[node]" )
{
    SECTION("Check that we can add nodes into the Graph and retrieve them")
    {
        graphPtr->addNode(node1);
        graphPtr->addNode(node2);
        Node* node1test = graphPtr->getNodePtr(1);
        Node* node3test = graphPtr->getNodePtr(3);

        REQUIRE(node1test != NULL);
        REQUIRE(node3test == NULL);
        REQUIRE(node1test->getIndex() == 1);
    }
    SECTION("Check that we can connect and disconnect nodes")
    {
        node1->connectTo(node2);
        node1->connectTo(node3);
        REQUIRE(node1->degree() == 2);
        node1->disconnect(node3);
        REQUIRE(node1->degree() == 1);
    }
}


TEST_CASE_METHOD(RetaContainerFixture, "RetaTuple is a working component", "[retatuple]" )
{
    SECTION("Initialization of a random  RetaTuple...")
    {
        RetaTuple* t = new RetaTuple();
        unsigned int order = 4;
        t->randomize(order);
        cout << *t << endl;
        REQUIRE(t->getData().size() == order);
    }
}


TEST_CASE_METHOD(RetaContainerFixture, "RetaTupleSet is a working component", "[retasubset]" )
{
    SECTION("Finding a tuple in the RetaTupleSet")
    {
        RetaTuple* t = new RetaTuple();
        t->randomize(3);
        s->add(t);
        REQUIRE(s->exists((const RetaTuple&)(*t)));
        RetaTuple* t2 = new RetaTuple();
        t2->randomize(3);
        // This new tuple is not added in the subset
        REQUIRE(!s->exists((const RetaTuple&)(*t2)));
        REQUIRE(s->exists((const RetaTuple&)(*t)));
    }


    SECTION("Initialization of a correct random RetaTupleSet")
    {
        unsigned int nbTuples = rand() % (int)1e5 + 1;
        unsigned int order = 3;
        for (unsigned int i=0; i < nbTuples; ++i)
        {
            RetaTuple* nt = new RetaTuple();
            nt->randomize(order);
            s->add(nt);
        }
        cout << *s;
        cout << nbTuples << " tuples of order " << order << " have been added" << endl;
        REQUIRE(nbTuples == s->size());

    }

    SECTION("Initialization of an inconsistent random RetaTupleSet (different orders in the same subset)")
    {
        RetaTuple* t = new RetaTuple();
        t->randomize(3);
        s->add(t);
        RetaTuple* t2 = new RetaTuple();
        t2->randomize(4);
        REQUIRE_THROWS(s->add(t2));
    }

    SECTION("Testing accessing stored tuple")
    {
        RetaTuple* t = new RetaTuple();
        t->randomize(3);
        RetaTuple* copy = new RetaTuple(*t); //copy has not H set
        t->setH(0.123);
        s->add(t);
        // Finding the tuple
        REQUIRE(s->exists((const RetaTuple&)(*copy)));
        cout << "Searching: " << *copy << " in: " << endl;
        cout << *s << endl;
        RetaTuple* found = s->find((const RetaTuple&)(*copy));
        REQUIRE(found != NULL);
        cout << "H value is: " << found->getH() << endl;
        RetaTuple* t2 = new RetaTuple();
        t2->randomize(3);
        // The new tuple can not be accessed (not added in subset)
        REQUIRE(s->find((const RetaTuple&)(*t2)) == NULL);
    }
}


TEST_CASE_METHOD(RetaContainerFixture, "RetaContainer is a working component", "[retacontainer]" )
{
    SECTION("Testing Container initialization")
    {
        unsigned int nbTuples = rand() % (int)1e5 + 1;
        unsigned int maxOrder = 6;
        for (unsigned int i=0; i < nbTuples; ++i)
        {
            unsigned int order = rand() % maxOrder + 1;
            if (order < 2)
            {
                order = 2;
            }
            RetaTuple* nt = new RetaTuple();
            nt->randomize(order);
            c->add(nt);
        }
        REQUIRE(c->getMaxOrder() == maxOrder);
    }

    SECTION("Testing for searching of a tuple")
    {
        srand (time(NULL));
        unsigned int nbTuples = rand() % (int)1e5 + 1;
        unsigned int randomPick = rand() % nbTuples;
        RetaTuple* tPicked = new RetaTuple();
        for (unsigned int i=0; i < nbTuples; ++i)
        {
            RetaTuple* t = new RetaTuple();
            t->randomize(3);
            if (i==randomPick)
            {
                *tPicked = *t; // tPicked will not have H set
                t->setH(0.12345);
                cout << "New test - Picking tuple: " << *tPicked << " from " << nbTuples << " stored in RetaContainer." << endl;
            }
            c->add(t);
        }

        // Check if the random picked is found
        REQUIRE(c->exists((const RetaTuple&)(*tPicked)));

        // Check that it is possible to retrieve a tuple in the container
        // and get it H value based on a created tuple
        RetaTuple* found = c->find((const RetaTuple&)(*tPicked));
        // Check that the tuple found as the H value
        REQUIRE(found != NULL);
        REQUIRE(found->getH() == 0.12345f);
    }

    SECTION("Combinatorial testing")
    {
        srand (time(NULL));
        unsigned int nbTuples = rand() % (int)10 + 1;
        RetaTuple* tPicked = new RetaTuple();
        cout << "Generating " << nbTuples << " tuples in the container for combinations..." << endl;
        for (unsigned int i=0; i < nbTuples; ++i)
        {
            RetaTuple* t = new RetaTuple();
            t->randomize(4);
            c->add(t);
        }
        c->getGlobalCombinations(s, 4, 4);
        cout << "Result: " << *s << endl;
        s->clear();

        c->getGlobalCombinations(s, 4, 3);
        cout << "Result: " << *s << endl;
        s->clear();

        c->getGlobalCombinations(s, 4, 2);
        cout << "Result: " << *s << endl;
    }
}
