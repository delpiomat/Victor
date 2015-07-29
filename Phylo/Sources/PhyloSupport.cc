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

    // ATTRIBUTES:
	double PhyloSupport::openGapPenalty=12;
	double PhyloSupport::extensionGapPenalty=3;
	double PhyloSupport::cSeq=0.8;
	double PhyloSupport::downs=999.9;
	double PhyloSupport::downa=999.9;
	double PhyloSupport::ups=999.9;
	double PhyloSupport::upa=999.9;
	unsigned int PhyloSupport::weightingScheme=0;

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

		GapFunction *gf = new AGPFunction(openGapPenalty, extensionGapPenalty);


		//Default Structure
		 Structure *str;
		 str = 0;

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


	vector<string> PhyloSupport::AlingSvsS(string seq1,string seq2,bool verbose){

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

		GapFunction *gf = new AGPFunction(openGapPenalty, extensionGapPenalty);


		//Default Structure
		 Structure *str;
		 str = 0;

		//Default ScoringScheme
		ScoringScheme *ss;

		//Global Align
		Align *a;
		cout << "Start Suboptimal Needleman-Wunsch alignments:" << endl;

		double suboptPenaltyMul=1;
		double suboptPenaltyAdd=1;

		unsigned int suboptNum=1;

		vector<Alignment> a2;


		ad = new SequenceData(n, seq1, seq2, "seq1", "seq2");

		ss = new ScoringS2S(&sub, ad, str, cSeq);

		a = new NWAlign(ad, gf, ss);

		a->setPenalties(suboptPenaltyMul, suboptPenaltyAdd);
		a2 = a->generateMultiMatch(suboptNum);
		if (a2.size() == 0)
			ERROR("No output alignments generated.", exception);
		//a2[0].cutTemplate(1);
		cout<<"NWAlign SvsS OK-------------------"<<endl;
		vector<string> result(2);

		cout<<"dimension of a2 "<<a2.size()<<" a2[0].getTarget() "<<a2[0].getTarget()<<endl<<" a2[0].getTemplate() "<<a2[0].getTemplate()<<endl<<" a2[0].getScore() "<<a2[0].getScore()<<endl<<""<<""<<""<<""<<endl;
		result[0]=a2[0].getTarget();
		result[1]=a2[0].getTemplate();
		cout<<"result "<<a2.size()<<" result[1] "<<result[1]<<endl<<" result[0] "<<result[0]<<endl;
		cout<<"NWAlign SvsS Oexit-------------"<<endl;
		return result;
}


	vector<string> PhyloSupport::AlingMultiSvsMultiS(vector<string> seq1,vector<string> seq2,vector <double> vWeigth1,vector <double> vWeigth2,bool verbose){


			string path = getenv("VICTOR_ROOT");
			if (path.length() < 3)
				cout << "Warning: environment variable VICTOR_ROOT is not set." << endl;

			string dataPath = path + "data/";

			//Default matrix
			string matrixFileName="blosum62.dat";
			matrixFileName = dataPath + matrixFileName;
			ifstream matrixFile(matrixFileName.c_str());
			if (!matrixFile)
				ERROR("Error opening substitution matrix file.", exception);
			SubMatrix sub(matrixFile);


			//cout<<"le sequenze sono lunghe"<<seq1[0].size()<<" "<<seq2[0].size()<<endl;
			vector< vector<double> > scoreM(seq1[0].size(), vector<double>(seq2[0].size()));

			for(unsigned int i=0;i<scoreM.size();i++)
			{
				for(unsigned int j=0;j<scoreM[0].size();j++)
				{
					scoreM[i][j]=0;
					for(unsigned int k=0;k<seq1.size();k++)
					{
						for(unsigned int w=0;w<seq2.size();w++)
						{
							//cout<<"score("<<i<<"-"<<j<<")= "<<seq1[k][i]<<" "<<seq2[w][j]<<" usando blosum da "<<sub.score[seq1[k][i]][seq2[w][j]]<<" peso "<<vWeigth1[k]<<" "<<vWeigth2[w]<<endl;
							scoreM[i][j]+=(sub.score[seq1[k][i]][seq2[w][j]])*(vWeigth1[k]*vWeigth2[w]);
						}
					}
					if(i!=j){
						scoreM[i][j]-=scoreM[i][j]/(2*(abs(i-j)));
					}
					//cout<<scoreM[i][j]<<"\t ";
				}
				//cout<<endl<<endl;
			}


			unsigned int i=scoreM.size()-1;
			unsigned int j=scoreM[0].size()-1;
			cout<<"\t \t \t ---------start multi Align -----------------score i j max"<<scoreM[i][j]<<endl;
			unsigned int tempI=i;
			unsigned int tempJ=j;
			for(unsigned int w=i-1;w<scoreM.size()-1;w++)
			{
				if(w!=0)
				{
					if(scoreM[w][j]>scoreM[tempI][j])
						tempI=w;
				}
			}
			for(unsigned int w=j-1;w<scoreM[0].size()-1;w++)
			{
				if(w!=0)
				{
					if(scoreM[i][w]>scoreM[i][tempJ])
						tempJ=w;
				}
			}
			bool up=false;
			bool vertical=false;
			if(scoreM[tempI][j]>scoreM[i][tempJ]){
				i=tempI;
				up=true;
			}
			else{
				vertical=true;
				j=tempJ;
			}

			if(up)
			{
				for(unsigned int h=0;h<(scoreM.size()-1)-i;h++){
					for(unsigned int w=0; w<seq1.size();w++ )
					{
						seq1[w]=PhyloSupport::insertGapPosition(seq1[w],seq1[w].size());
					}
				}
			}else
			{
				for(unsigned int h=0;h<(scoreM[0].size()-1)-j;h++){
					for(unsigned int w=0; w<seq2.size();w++ )
					{
						seq2[w]=PhyloSupport::insertGapPosition(seq2[w],seq2[w].size());
					}
				}
			}


			while(i!=0 && j!=0){

				//cout<<"max{"<<scoreM[i-1][j]<<" "<<scoreM[i-1][j-1]<<" "<<scoreM[i][j-1]<<"}  "<<endl;

				if(scoreM[i-1][j-1]>=scoreM[i-1][j] && scoreM[i-1][j-1]>=scoreM[i][j-1])
				{
						i--;
						j--;
						//cout<<"diagonal "<<scoreM[i][j]<<endl;
				}
				else if(scoreM[i-1][j]>scoreM[i-1][j-1] && scoreM[i-1][j]>scoreM[i][j-1])
				{
					for(unsigned int w=0; w<seq1.size();w++ )
					{
						seq1[w]=PhyloSupport::insertGapPosition(seq1[w],i);
					}
					i--;
					//cout<<"up "<<scoreM[i][j]<<endl;
				}
				else{

					for(unsigned int w=0; w<seq2.size();w++ )
					{
						//cout<<seq2[w]<<endl;
						seq2[w]=PhyloSupport::insertGapPosition(seq2[w],j);
					//	cout<<"insert GAP in position "<<j<<endl;
						//cout<<seq2[w]<<endl;
					}
					j--;
				//	cout<<"left "<<scoreM[i][j]<<endl;
				}
				//cout<<endl;

			}
			cout<<"find the good way"<<endl;

			vector<string> seqFinal(seq1.size()+seq2.size());

			cout<<"seq1 "<<seq1.size()<<" seq2 "<<seq2.size()<<endl;
			cout<<"seqFinal size "<<seqFinal.size()<<endl;

			cout<<"seq1[0] size "<<seq1[0].size()<<" seq2[0] size"<<seq2[0].size()<<endl;


			bool onlyGap=true;
			while(seq1[0].size()<seq2[0].size() && onlyGap){

				for(unsigned int g=0;g<seq2.size();g++){
					//cout<<seq2[g]<<endl<<endl;
					if(seq2[g][0]!='-'){
						onlyGap=false;
					}
				}

				if(!onlyGap){
					for(unsigned int w=0; w<seq1.size();w++)
					{
						//cout<<"uno dei seq1 "<<seq1[w]<<endl<<endl;
						seq1[w]=PhyloSupport::insertGapPosition(seq1[w],0);
					}
				}
				else{
					for(unsigned int g=0;g<seq2.size();g++){
						//cout<<"uno dei seq1 "<<seq1[w]<<endl<<endl;
						seq2[g]=seq2[g].substr(1);
					}
				}
				onlyGap=true;
			}
			cout<<"seq1[0] "<<seq1[0].size()<<" seq2[0] "<<seq2[0].size()<<endl;

			onlyGap=true;
			while(seq2[0].size()<seq1[0].size() && onlyGap){
				for(unsigned int g=0;g<seq1.size();g++){
					//cout<<seq1[g]<<endl<<endl;
					if(seq1[g][0]!='-'){
						onlyGap=false;
					}
				}

				if(!onlyGap){
					for(unsigned int w=0; w<seq2.size();w++)
					{
						seq2[w]=PhyloSupport::insertGapPosition(seq2[w],0);
					}
				}
				else{
					for(unsigned int g=0;g<seq1.size();g++){
						seq1[g]=seq1[g].substr(1);
					}
				}
				onlyGap=true;
			}

			//cout<<"seq1"<<endl;
			for(unsigned int w=0; w<seq1.size();w++ )
			{
				//cout<<"della seq size= "<<seq1[w].size()<<endl;
				//cout<<seq1[w]<<endl<<endl;
				seqFinal[w]=seq1[w];
			}
			//cout<<endl;
			//cout<<"seq2"<<endl;
			for(unsigned int w=0; w<seq2.size();w++ )
			{
				//cout<<"della seq size= "<<seq2[w].size()<<endl;
				//cout<<seq2[w]<<endl<<endl;

				seqFinal[w+seq1.size()]=seq2[w];
			}
			cout<<"\t \t \t ---------finish multi Align prima del return----------------"<<endl;
			return seqFinal;
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


  	  string PhyloSupport::insertGapPosition(string seq, int position){
  		  string part1="";
  		  string part2="";
  		  part1=seq.substr(0,position);
  		  part2=seq.substr(position);
  		  seq=part1+"-"+part2;
  		  return seq;
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


		string PhyloSupport::intToString ( int num )
		  {
		     ostringstream ss;
		     ss << num;
		     return ss.str();
		  }


}} // namespace
