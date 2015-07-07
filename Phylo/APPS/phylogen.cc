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
            << "\n                     \t --gf=0: UPGMA (default)."
            << "\n                     \t --gf=1: NJ."
            << "\n   [--cSeq <double>] \t Coefficient for sequence alignment (default = 0.80)"
            << "\n   [--cStr <double>] \t Coefficient for structural alignment (default = 0.20)"
            << "\n"
			<< "\n   [--kimura]     	\t use kimura formula for calculate distance of pairwise sequence"
            << "\n   [--verbose]       	\t Verbose mode"
            << "\n" << endl;
}


/// Show command line options and help text.

int main(int argc, char **argv) {

    string inputFileName, outputFileName, matrixFileName, matrixStrFileName;
    double openGapPenalty, extensionGapPenalty;
    unsigned int  gapFunction, cluster;
    double cSeq, cStr;
    bool verbose=false;
    bool kimura=false;
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

	kimura = getArg("-kimura", argc, argv);
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
		tree.upgma(aliSec,kimura,verbose);
		cout<<"##########################################Tree create with PhyloTreeUPGMA#############################"<<endl;
		out=tree.printNewickTree();
	}
	else if(cluster==1){
		cout<<"##########################################Starting PhyloTreeNJ########################################"<<endl;
		tree.neighborJoining(aliSec,kimura,verbose);
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

	cout<<"------------------------End-------------------------"<<endl;

	return 0;

}
