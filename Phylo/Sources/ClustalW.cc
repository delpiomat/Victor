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
#include <SubMatrix.h>
#include <Alignment.h>
#include <AlignmentData.h>
#include <SequenceData.h>
#include <NWAlign.h>
#include <AGPFunction.h>
#include <ScoringS2S.h>

#include <Profile.h>
#include <HenikoffProfile.h>
#include <PSICProfile.h>
#include <SeqDivergenceProfile.h>

#include <NewickTree.h>
#include <PhyloSupport.h>
#include <ClustalW.h>

#include <sstream>
#include <string>
#include <vector>

// Global constants, typedefs, etc. (to avoid):
using namespace Victor::Align2;
using namespace Victor;
using namespace std;

namespace Victor { namespace Phylo{


	ClustalW::ClustalW(){
	}

	ClustalW::ClustalW(NewickTree gT) {

		cout<<"constructor"<<endl;
		guideTree=gT;
        progressiveAlign();
	}


	ClustalW::~ClustalW() {
	}


	void ClustalW::progressiveAlign(){

		cout<<" Progressive alignament Start "<<endl;
		vector<string> tmpV(2);
		vector<double> tmpWeigth(2);
		int tokenSize=60;


		vector <iNode*> nodeTree(guideTree.getNumberOfLeaf());
		for(unsigned int i=0;i<guideTree.getNumberOfLeaf();i++){
			nodeTree[i]=guideTree.getLeafInPosition(i).getRoot();
		}

		vector <string> seqV(1);
		vector <double> weigthV(1);
		for(unsigned int j=1;j<guideTree.getNumberOfLeaf();j++){
			//cout<<"\t --------- NEW START ---------------"<<j<<endl;

			for(unsigned int i=0;i<nodeTree.size();i++){

				if(j==guideTree.getNumberOfLeaf() && i==0)
				{
					if(nodeTree[0]->allignSeq[0].size()<nodeTree[1]->allignSeq[0].size())
						tmpV=PhyloSupport::AlingMultiSvsMultiS2(nodeTree[0]->allignSeq,nodeTree[1]->allignSeq,nodeTree[0]->weigthV,nodeTree[1]->weigthV,false,tokenSize);
					else
						tmpV=PhyloSupport::AlingMultiSvsMultiS2(nodeTree[1]->allignSeq,nodeTree[0]->allignSeq,nodeTree[1]->weigthV,nodeTree[0]->weigthV,false,tokenSize);
					nodeTree[1]->ClustalW=true;
					nodeTree[0]->allignSeq=tmpV;
				}
				else if(j==1){
					seqV[0]=nodeTree[i]->seq;
					weigthV[0]=nodeTree[i]->weigth;
					nodeTree[i]->allignSeq=seqV;
					nodeTree[i]->weigthV=weigthV;
				}
				else if(j==2 && nodeTree[i]->parent->numberOfChildLeaf==2){
					tmpV=PhyloSupport::AlingSvsS(nodeTree[i]->seq,nodeTree[i]->seq);
					weigthV[0]=nodeTree[i]->weigth;
					if(weigthV.size()==1)
						weigthV.push_back(nodeTree[i+1]->weigth);
					else
						weigthV[1]=nodeTree[i+1]->weigth;

					nodeTree.erase(nodeTree.begin()+i+1);
					nodeTree[i]->parent->allignSeq=tmpV;
					nodeTree[i]->parent->weigthV=weigthV;
					nodeTree[i]=nodeTree[i]->parent;
				}
				else if(j>2 && nodeTree[i]->parent->numberOfChildLeaf==j && !nodeTree[i]->ClustalW){

					if(nodeTree[i]->isLeft){
						if(nodeTree[i]->allignSeq[0].size()<nodeTree[i]->parent->right->allignSeq[0].size())
							tmpV=PhyloSupport::AlingMultiSvsMultiS2(nodeTree[i]->allignSeq,nodeTree[i]->parent->right->allignSeq,nodeTree[i]->weigthV,nodeTree[i]->parent->right->weigthV,false,tokenSize);
						else
							tmpV=PhyloSupport::AlingMultiSvsMultiS2(nodeTree[i]->parent->right->allignSeq,nodeTree[i]->allignSeq,nodeTree[i]->parent->right->weigthV,nodeTree[i]->weigthV,false,tokenSize);
						nodeTree[i]->parent->right->ClustalW=true;
						vector <double> weigthTMP(nodeTree[i]->weigthV.size()+nodeTree[i]->parent->right->weigthV.size());
						for(unsigned int y=0;y<nodeTree[i]->weigthV.size();y++){
							weigthTMP[y]=nodeTree[i]->weigthV[y];
						}
						for(unsigned int y=0;y<nodeTree[i]->parent->right->weigthV.size();y++){
							weigthTMP[y+nodeTree[i]->weigthV.size()]=nodeTree[i]->parent->right->weigthV[y];
						}
						nodeTree[i]->parent->weigthV=weigthTMP;
					}
					else{//is right child
						if(nodeTree[i]->allignSeq[0].size()<nodeTree[i]->parent->left->allignSeq[0].size())
							tmpV=PhyloSupport::AlingMultiSvsMultiS2(nodeTree[i]->allignSeq,nodeTree[i]->parent->left->allignSeq,nodeTree[i]->weigthV,nodeTree[i]->parent->right->weigthV,false,tokenSize);
						else
							tmpV=PhyloSupport::AlingMultiSvsMultiS2(nodeTree[i]->parent->left->allignSeq,nodeTree[i]->allignSeq,nodeTree[i]->parent->right->weigthV,nodeTree[i]->weigthV,false,tokenSize);
						nodeTree[i]->parent->left->ClustalW=true;
						vector <double> weigthTMP(nodeTree[i]->weigthV.size()+nodeTree[i]->parent->left->weigthV.size());
						for(unsigned int y=0;y<nodeTree[i]->weigthV.size();y++){
							weigthTMP[y]=nodeTree[i]->weigthV[y];
						}
						for(unsigned int y=0;y<nodeTree[i]->parent->left->weigthV.size();y++){
							weigthTMP[y+nodeTree[i]->weigthV.size()]=nodeTree[i]->parent->left->weigthV[y];
						}
						nodeTree[i]->parent->weigthV=weigthTMP;
					}
					nodeTree[i]=nodeTree[i]->parent;
					nodeTree[i]->allignSeq=tmpV;
				}
				else if(nodeTree[i]->ClustalW){
					nodeTree.erase(nodeTree.begin()+i);
				}

			}
		}


		string outString=ClustalW::printClustalWFromat(nodeTree[0]->allignSeq);
		cout<<endl<<outString<<endl;

		//write clustalW
	    ofstream outFile("out.clustalw");
	    if(!outFile) {
	        cout<<"Error During creation File out.clustalw";
	    }
	    outFile<<outString<<endl;

	    outFile.close();
		cout<<"endl clustaW"<<endl;

	}



	string ClustalW::printClustalWFromat(vector<string> seq){
		string txt="";
		string tmp="";
		unsigned int j=0;
		unsigned int index=0;
		unsigned int dim=50;
		while(j<seq[0].size()){
			for(unsigned int i=0;i<seq.size();i++){
			txt+="seq"+PhyloSupport::intToString(i)+" \t \t ";
				for(j=dim*index;j<=dim+dim*index && j<seq[i].size();j++){
					txt+=seq[i][j];
				}
			txt+=" "+PhyloSupport::intToString(j);
			cout<<endl;
			txt+="\n";
			}
		txt+="\n\n";
		index++;
		}


		return txt;

	}



}} // namespace
