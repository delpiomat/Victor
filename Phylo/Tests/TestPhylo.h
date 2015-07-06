/*
 * TestPhylo.h
 *
 *  Created on: 02 07 2015
 *      Author: Matteo Del Pioluogo
 */


#include <iostream>
#include <cppunit/TestFixture.h>
#include <cppunit/TestAssert.h>
#include <cppunit/TestCaller.h>
#include <cppunit/TestSuite.h>
#include <cppunit/TestCase.h>
#include <PhyloSupport.h>
using namespace std;
using namespace Victor;


class TestPhylo : public CppUnit::TestFixture {
private:
private:

public:

    TestPhylo(){
    }


    static CppUnit::Test *suite() {
        CppUnit::TestSuite *suiteOfTests = new CppUnit::TestSuite("TestAlign");
        suiteOfTests->addTest(new CppUnit::TestCaller<TestPhylo>("Test1 - Calculate distance of 2 seq.",
                &TestPhylo::testPhylo_A));
        suiteOfTests->addTest(new CppUnit::TestCaller<TestPhylo>("Test2 - .",
                &TestPhylo::testPhylo_B));
        suiteOfTests->addTest(new CppUnit::TestCaller<TestPhylo>("Test3 - .",
                &TestPhylo::testPhylo_C));

        return suiteOfTests;
    }

    /// Setup method

    void setUp() {
    }

    /// Teardown method

    void tearDown() {
    }

protected:

    void testPhylo_A() {
        CPPUNIT_ASSERT(1==0);
    }

    void testPhylo_B() {


        CPPUNIT_ASSERT(7 == 1);
    }

    void testPhylo_C() {

    	CPPUNIT_ASSERT(1 == 1);
    }

};
