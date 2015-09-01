/*
 * TestPhylogen.cc
 *
 *  Created on: 02 07 2015
 *      Author: Matteo Del Pioluogo
 */

#include <iostream>
#include <cppunit/TestSuite.h>
#include <cppunit/ui/text/TestRunner.h>

#include <TestPhylo.h>

using namespace std;

int main(){

	string a="AACCAA";
	string b="AACCAA";
	string c="AACCAA";
	string o="CC";
	string t="CC";
	vector <string> num(2);
	vector <string> letter(3);
	num[0]=o;
	num[1]=t;
	letter[0]=a;
	letter[1]=b;
	letter[2]=c;
	vector <double> w1(2);
	vector <double> w2(3);
	w1[0]=1;
	w1[1]=1;
	w2[0]=1;
	w2[1]=1;
	w2[2]=1;
	vector <string> result=Phylo::PhyloSupport::AlingMultiSvsMultiS(num,letter,w1,w2,false,1);
	cout<<result[0][0]<<"asdasdasdsdsdd"<<endl;

	for(unsigned int i=0;i<result.size();i++) {
		for(unsigned int j=0;j<result[i].size();j++) {
			cout<<result[i][j];
		}
		cout<<endl;
	}

	CppUnit::TextUi::TestRunner runner;
	cout << "Creating Test Suites: OK" << endl;
        runner.addTest(TestPhylo::suite());
	cout<< "Running the unit tests."<<endl;
	runner.run();

	return 0;
}
























































