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

// --*- C++ -*------x-----------------------------------------------------------
//
//
// Description:
//
// -----------------x-----------------------------------------------------------
#include <Alignment.h>
#include <AlignmentData.h>
#include <SequenceData.h>
#include <SubMatrix.h>
#include <NWAlign.h>
#include <AGPFunction.h>
#include <ScoringS2S.h>

#include <GetArg.h>

#include <PhyloTree.h>
#include <PhyloTreeUPGMA.h>
#include <PhyloTreeNJ.h>

#include <NewickTree.h>
#include <PhyloSupport.h>

#include <string>
#include <vector>
#include <sstream>
#include <ctime>

using namespace Victor::Align2;
using namespace Victor::Phylo;
using namespace Victor;
using namespace std;

vector<Alignment> calcAlignmentV(Alignment *aliSec, vector<vector<double> > &distance){
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
	cout << "\nSuboptimal Needleman-Wunsch alignments:\n" << endl;

	double suboptPenaltyMul=1;
	double suboptPenaltyAdd=1;

	unsigned int suboptNum=1;

	vector<Alignment> a2;

	Alignment a3;

	vector <Alignment> alignV(aliSec->size()+1);//vector of align

	double score=0;
	//size not count target. so +1
	for(unsigned int index=0;index<aliSec->size()+1;index++){
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
			score=PhyloTree::distanceCalcTwoSeq(a2[0].getTarget(),a2[0].getTemplate(0));
			distance[index][j]=score;
			cout<<"SCORE="<<distance[index][j]<<" ";
			if(j==0)
				newAli = a2[0];
			else
				newAli.addAlignment(a2[0]);
		}//end for j
		cout<<"ultima seq1 "<<seq1Name<<endl;
		cout<<"corrisponde al target "<<newAli.getTargetName()<<endl;
		cout<<"ultima seq2 "<<seq2Name<<endl;
		cout<<"fine j"<<"AlignV ha size "<<alignV.size()<<endl;

		cout<<"ongi align vale?"<<newAli.size()<<endl;

		alignV[index]=newAli;
	}//end for index

	cout<<"NWAlign ok"<<endl;
	return alignV;
}


/// Show command line options and help text.

void
sShowHelp() {
    cout << "\nP Phylogenetic Tree GENERATOR" //This program read FASTA file and calculate phylogenetic tree. Use NWALIGN for align.
            << "\nThis program read FASTA file and calculate phylogenetic tree. Use NWALIGN for align."
            << "\nOptions:"
            << "\n"
            << "\n * [--in <name>]     \t Name of input FASTA file"
            << "\n   [--out <name>]    \t Name of output file (default = to screen)"
            << "\n   [-m <name>]       \t Name of substitution matrix file (default = blosum62.dat)"
			<< "\n"
            << "\n   [--gf <0|1|2>]    \t Gap function (default = 0, i.e. AGP)"
            << "\n                     \t --gf=0: AGP (Affine Gap Penalty) function (default)."
            << "\n                     \t --gf=1: VGP (Variable Gap Penalty) function."
            << "\n                     \t --gf=2: VGP2 (Variable Gap Penalty) function."
            << "\n   [-o <double>]     \t Open gap penalty (default = 12.00)"
            << "\n   [-e <double>]     \t Extension gap penalty (default = 3.00)"
            << "\n"
            << "\n   [--cluster <0|1>] \t CLuster function (default = 0, i.e. UPGMA)"
            << "\n                     \t --gf=0: UPGMA (Affine Gap Penalty) Cluster (default)."
            << "\n                     \t --gf=1: NJ (Variable Gap Penalty) Cluster."
            << "\n   [--cSeq <double>] \t Coefficient for sequence alignment (default = 0.80)"
            << "\n   [--cStr <double>] \t Coefficient for structural alignment (default = 0.20)"
            << "\n"
            << "\n   [--verbose]       \t Verbose mode"
            << "\n" << endl;
}


/// Show command line options and help text.

int main(int argc, char **argv) {

    string inputFileName, outputFileName, matrixFileName, matrixStrFileName;
    double openGapPenalty, extensionGapPenalty;
    unsigned int  gapFunction, cluster;
    double cSeq, cStr;
    bool verbose;
    struct tm* newtime;
    time_t t;

    // --------------------------------------------------
    // 0. Treat options
    // --------------------------------------------------
    if (getArg("h", argc, argv)) {
        sShowHelp();
        return 1;
    }


    getArg("-in", inputFileName, argc, argv, "cyc1_input.fasta");
	getArg("-out", outputFileName, argc, argv, "!");

	getArg("m", matrixFileName, argc, argv, "blosum62.dat");

	getArg("-cluster", cluster, argc, argv, 0);

	getArg("-gf", gapFunction, argc, argv, 0);
	getArg("o", openGapPenalty, argc, argv, 12.00);
	getArg("e", extensionGapPenalty, argc, argv, 3.00);

	getArg("-cSeq", cSeq, argc, argv, 0.80);
	getArg("-cStr", cStr, argc, argv, 0.20);

	verbose = getArg("-verbose", argc, argv);


    // --------------------------------------------------
    // 1. Load FASTA
    // --------------------------------------------------
	Alignment aliSec;
		string path = getenv("VICTOR_ROOT");
		if (path.length() < 3)
			cout << "Warning: environment variable VICTOR_ROOT is not set." << endl;
		string dataPath = path + "data/";
		string finalPath=dataPath+inputFileName;

		ifstream secFile(finalPath.c_str());
				if (!secFile)
					ERROR("Error opening secondary structure FASTA file.", exception);
		aliSec.loadFasta(secFile);
		if (aliSec.size() < 1)
			ERROR("Secondary structure FASTA file must contain two sequences.", exception);



    // --------------------------------------------------
    // 3. Load data
    // --------------------------------------------------

	NewickTree tree;
	string out="Error Tree";
	if(cluster==0){
		cout<<"##########################################Starting PhyloTreeUPGMA#####################################"<<endl;
		tree.upgma(aliSec,true);
		cout<<"##########################################Tree create with PhyloTreeUPGMA#############################"<<endl;
		out=tree.printNewickTree();
	}
	else if(cluster==1){
		cout<<"##########################################Starting PhyloTreeNJ########################################"<<endl;
		tree.neighborJoining(aliSec);
		cout<<"##########################################Tree create with PhyloTreeNJ################################"<<endl;
		out=tree.printNewickTree();
	}
	else
	{
		cout<<"invalid cluster, select [--cluster <0|1>]";
	}

    // --------------------------------------------------
    // 2. Output Tree
    // --------------------------------------------------
	if(outputFileName=="!"){
		cout<<"stampiamo l'albero creato "<<endl;
		cout<<out<<endl;
	}
	else{
		cout<<"stampiamo l'albero creato nel file data/"<<outputFileName<<endl;
		cout<<out<<endl;
	}

	string sequ1="AAAA--BAAAAAACAAA--------AAAAACAAA";
	string sequ2="A-AA-ABAAAAACCAAAAA--AABAAAABBBBBB";
	string sequ3="AAAABAAAAAACAAAAAAAACAAA";
	string sequ4="AAAACCAAAAAAABAAAABBBBBB";
	cout<<"Kimura distance of "<<sequ1<<" and "<<sequ2<<" ="<<PhyloSupport::distanceCalcTwoSeqKimura(sequ1,sequ2)<<endl;
	cout<<"Kimura distance of "<<sequ2<<" and "<<sequ1<<" ="<<PhyloSupport::distanceCalcTwoSeqKimura(sequ2,sequ1)<<endl;
	cout<<"Phylo distance of "<<sequ1<<" and "<<sequ2<<" ="<<PhyloSupport::distanceCalcTwoSeq(sequ1,sequ2)<<endl;
	cout<<"Phylo distance of "<<sequ2<<" and "<<sequ1<<" ="<<PhyloSupport::distanceCalcTwoSeq(sequ2,sequ1)<<endl;
	cout<<"Align2 Victor distance of "<<sequ1<<" and "<<sequ2<<" ="<<1-aliSec.calculatePairwiseIdentity(sequ1,sequ2)<<endl;
	cout<<"Align2 Victor distance of "<<sequ2<<" and "<<sequ1<<" ="<<1-aliSec.calculatePairwiseIdentity(sequ2,sequ1)<<endl;

	cout<<"Kimura distance of "<<sequ3<<" and "<<sequ4<<" ="<<PhyloSupport::distanceCalcTwoSeqKimura(sequ3,sequ4)<<endl;
	cout<<"Kimura distance of "<<sequ4<<" and "<<sequ3<<" ="<<PhyloSupport::distanceCalcTwoSeqKimura(sequ4,sequ3)<<endl;
	cout<<"Phylo distance of "<<sequ3<<" and "<<sequ4<<" ="<<PhyloSupport::distanceCalcTwoSeq(sequ3,sequ4)<<endl;
	cout<<"Phylo distance of "<<sequ4<<" and "<<sequ3<<" ="<<PhyloSupport::distanceCalcTwoSeq(sequ4,sequ3)<<endl;
	cout<<"Align2 Victor distance of "<<sequ3<<" and "<<sequ4<<" ="<<1-aliSec.calculatePairwiseIdentity(sequ3,sequ4)<<endl;
	cout<<"Align2 Victor distance of "<<sequ4<<" and "<<sequ3<<" ="<<1-aliSec.calculatePairwiseIdentity(sequ4,sequ3)<<endl;


	return 0;

}
