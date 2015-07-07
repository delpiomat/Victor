/*  This file is part of Victor.

    Victor is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Victor is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Victor.  If not, see <http://www.gnu.org/licenses/>.
 */


// Includes:
#include <stdio.h>
#include <iostream>
#include <string>
#include <vector>
#include "PhyloSupport.h"
#include <SubMatrix.h>
#include <Alignment.h>
#include <AlignmentData.h>
#include <SequenceData.h>
#include <NWAlign.h>
#include <AGPFunction.h>
#include <ScoringS2S.h>
#include <IoTools.h>
#include <math.h>

// Global constants, typedefs, etc. (to avoid):
using namespace Victor::Align2;
using namespace Victor;
using std::string;
using std::vector;
using std::ifstream;
using std::cout;
using std::endl;


namespace Victor { namespace Phylo{


    /** @brief
     *
     **/
	// CONSTRUCTORS/DESTRUCTOR:

	/**
	 *@Description Basic contructor
	 */
	PhyloSupport::PhyloSupport(){

	}


	/**
	 *@Description Basic destructor
	 */
	PhyloSupport::~PhyloSupport() {

	}


	vector<Alignment> PhyloSupport::calcAlignmentV(Alignment *aliSec, vector<vector<double> > &distance , bool ktuples,bool verbose){
		string seq1Name, seq2Name, seq1, seq2, sec1, sec2;

		Alignment newAli;
		string path = getenv("VICTOR_ROOT");
		if (path.length() < 3)
			cout << "Warning: environment variable VICTOR_ROOT is not set." << endl;

		string dataPath = path + "data/";

		//n=4? no idea. Default in APP/subali.cc
		int n=4;

		AlignmentData *ad;

		//Default matrix
		string matrixFileName="blosum62.dat";
		matrixFileName = dataPath + matrixFileName;
		ifstream matrixFile(matrixFileName.c_str());
		if (!matrixFile)
			ERROR("Error opening substitution matrix file.", exception);
		SubMatrix sub(matrixFile);

		//Default gap function
		double openGapPenalty=12;
		double extensionGapPenalty=3;
		GapFunction *gf = new AGPFunction(openGapPenalty, extensionGapPenalty);


		//Default Structure
		 Structure *str;
		 str = 0;

		//Default csep
		 double cSeq;
		 cSeq = 1.00;

		//Default ScoringScheme
		ScoringScheme *ss;

		//Global Align
		Align *a;
		cout << "Start Suboptimal Needleman-Wunsch alignments:" << endl;

		double suboptPenaltyMul=1;
		double suboptPenaltyAdd=1;

		unsigned int suboptNum=1;

		vector<Alignment> a2;

		Alignment a3;

		vector <Alignment> alignV(aliSec->size()+1);//vector of align

		double score=0;
		//size not count target. so +1
		for(unsigned int index=0;index<aliSec->size()+1;index++){
			if(verbose)
				cout<<" index="<<index;
			if(index==aliSec->size()){
				seq1 = Alignment::getPureSequence(aliSec->getTarget());
				seq1Name = aliSec->getTargetName();
			}
			else{
				seq1 = Alignment::getPureSequence(aliSec->getTemplate(index));
				seq1Name = aliSec->getTemplateName(index);
			}
			for(unsigned int j=0;j<aliSec->size()+1;j++){
				if(verbose)
					cout<<" j="<<index;
				if(j==aliSec->size()){
					seq2 = Alignment::getPureSequence(aliSec->getTarget());
					seq2Name = aliSec->getTargetName();
				}else{
					seq2 = Alignment::getPureSequence(aliSec->getTemplate(j));
					seq2Name = aliSec->getTemplateName(j);
				}

				ad = new SequenceData(n, seq1, seq2, seq1Name, seq2Name);

				ss = new ScoringS2S(&sub, ad, str, cSeq);

				a = new NWAlign(ad, gf, ss);

				a->setPenalties(suboptPenaltyMul, suboptPenaltyAdd);
				a2 = a->generateMultiMatch(suboptNum);
				if (a2.size() == 0)
					ERROR("No output alignments generated.", exception);

				a2[0].cutTemplate(1);
				if(ktuples==true)
				{
					score=PhyloSupport::distanceCalcTwoSeqktuples(a2[0].getTarget(),a2[0].getTemplate(0));
				}
				else
					score=PhyloSupport::distanceCalcTwoSeq(a2[0].getTarget(),a2[0].getTemplate(0));
				distance[index][j]=score;
				if(verbose)
					cout<<"SCORE="<<distance[index][j]<<" ";
				if(j==0)
					newAli = a2[0];
				else
					newAli.addAlignment(a2[0]);
			}//end for j
			if(verbose){
				cout<<" seq1= "<<seq1Name<<endl;
				cout<<" seq2= "<<seq2Name<<endl;
			}
			alignV[index]=newAli;
		}//end for index

		cout<<"NWAlign ok"<<endl;
		return alignV;
	}

    /**
     *  Calculate the distance for two sequence. distance = 1-identity(PID)
     *  Percent sequence identity (PID)=Calculates the identity as the number of identical positions,
     *  it gives a positive value when the two sequences have a gap in the same position,
     *  divided by the length of the alignment.
     *
     */
   double PhyloSupport::distanceCalcTwoSeq(string seq1,string seq2){
        if (seq1.length() != seq2.length())
            cout << "Warning: sequence lengths do not match:\n"
                << "seq1 = " << seq1.length() << "\n"
            << "seq2 = " << seq2.length() << "\n";

        double tmp = 0;
        for (unsigned int i = 0; i < seq1.length(); i++){
            if ( seq1[i] == seq2[i])
                tmp++;
        }

        return  (1- (tmp / seq1.length()));
    }


   /**
    *  Calculate the distance for two sequence by ktuples formula.
    *  return number of k-tuples (number of residues identical) minus
	*	constant penality for each gap.
    *
    */
  double PhyloSupport::distanceCalcTwoSeqktuples(string seq1,string seq2){
       if (seq1.length() != seq2.length())
           cout << "Warning: sequence lengths do not match:\n"
               << "seq1 = " << seq1.length() << "\n"
           << "seq2 = " << seq2.length() << "\n";

       double penality = 0;
       double d=0;
       for (unsigned int i = 0; i < seq1.length(); i++){
    	   if(seq1[i]=='-' || seq2[i]=='-')
    		   penality+=1;
    	   else if( seq1[i] == seq2[i] )
               d++;
       }
       penality=(penality/seq1.length())*0.9*d;

       return 1-((d-penality)/seq1.length());
   }


	std::vector<std::string> &PhyloSupport::split(const std::string &s, char delim, std::vector<std::string> &elems) {
		    std::stringstream ss(s);
		    std::string item;
		    while (std::getline(ss, item, delim)) {
		        elems.push_back(item);
		    }
		    return elems;
		}

		std::vector<std::string> PhyloSupport::split(const std::string &s, char delim) {
		    std::vector<std::string> elems;
		    split(s, delim, elems);
		    return elems;
		}


		double PhyloSupport::calcDivR(vector<double> vDist){
			double tmpR=0;
			for(unsigned int i=0;i<vDist.size();i++){
				tmpR+=vDist[i];
			}
			return tmpR;
		}
	    /**
	     * Create a string of matrix value;
	     *
	     */
		void PhyloSupport::printMatrix(vector<vector<double> > &distance){
			string sMatrix="";
			//std::setprecision(3)//for minus precision
			for(unsigned int i=0; i<distance.size();i++){
				for(unsigned int j=0; j<distance.size();j++){
					cout<<distance[i][j]<<" ";
				}
				cout<<endl;
			}
		}



}} // namespace
