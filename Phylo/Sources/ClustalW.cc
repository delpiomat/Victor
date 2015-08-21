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

		/*cout<<" Progressive alignament Start "<<endl;
		vector<string> tmpV(2);
		vector<double> tmpWeigth(2);

		cout<<"leaf number "<<guideTree.getNumberOfLeaf()<<endl;

		vector <iNode*> nodeTree(guideTree.getNumberOfLeaf());
		for(unsigned int i=0;i<guideTree.getNumberOfLeaf();i++){
			nodeTree[i]=guideTree.getLeafInPosition(i).getRoot();
		}

		vector <string> seqV(1);
		vector <double> weigthV(1);
		for(unsigned int j=1;j<guideTree.getNumberOfLeaf()+1;j++){
			cout<<"\t --------- NEW START ---------------"<<j<<endl;

			for(unsigned int i=0;i<nodeTree.size();i++){
				cout<<"\t --------- i= "<<i<<" j= ---------------"<<j<<endl;

				if(j==guideTree.getNumberOfLeaf() && i==0)
				{
					tmpV=PhyloSupport::AlingMultiSvsMultiS(nodeTree[0]->allignSeq,nodeTree[1]->allignSeq,nodeTree[0]->weigthV,nodeTree[1]->weigthV);
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
					cout<<"temporaneo weigthV.size = "<<weigthV.size()<<endl;
					if(weigthV.size()==1)
						weigthV.push_back(nodeTree[i+1]->weigth);
					else
						weigthV[1]=nodeTree[i+1]->weigth;
					cout<<"temporaneo weigthV.size =  "<<weigthV.size()<<endl;
					cout<<"size di node tree prima di erase "<<nodeTree.size()<<endl;
					nodeTree.erase(nodeTree.begin()+i+1);
					cout<<"size di node tree dopo erase "<<nodeTree.size()<<endl;
					nodeTree[i]->parent->allignSeq=tmpV;
					nodeTree[i]->parent->weigthV=weigthV;
					nodeTree[i]=nodeTree[i]->parent;
					cout<<"size di node tree dopo cambio padre "<<nodeTree.size()<<endl;
				}
				else if(j>2 && nodeTree[i]->parent->numberOfChildLeaf==j && !nodeTree[i]->ClustalW){

					if(nodeTree[i]->isLeft){
						tmpV=PhyloSupport::AlingMultiSvsMultiS(nodeTree[i]->allignSeq,nodeTree[i]->parent->right->allignSeq,nodeTree[i]->weigthV,nodeTree[i]->parent->right->weigthV);
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
						tmpV=PhyloSupport::AlingMultiSvsMultiS(nodeTree[i]->allignSeq,nodeTree[i]->parent->left->allignSeq,nodeTree[i]->weigthV,nodeTree[i]->parent->right->weigthV);
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
					cout<<"size di node tree prima di erase dentro un caso gia FATTO "<<nodeTree.size()<<" con i= "<<i<<endl;
					nodeTree.erase(nodeTree.begin()+i);
					cout<<"size di node tree dopo erase "<<nodeTree.size()<<endl;
				}

			}
		}



		cout<<"\t nodeTree.size() = "<<nodeTree.size()<<endl;
		cout<<"size di [0] "<<nodeTree[0]->allignSeq.size()<<endl;
		cout<<"seq di [0] "<<endl<<nodeTree[0]->allignSeq[0]<<endl<<endl;

		cout<<"la5"<<endl<<nodeTree[0]->allignSeq[5]<<endl;

		cout<<endl<<ClustalW::printClustalWFromat(nodeTree[0]->allignSeq)<<endl;

		cout<<"endl clustaW"<<endl;*/

		/*string o="TDVIYQIVTDRFADGDRTNNPAGDAFSGDRSNLKLYFGGDWQGIID";
		string t="TDVIYQIVTDRFSDGNPGNNPSGAIFSQNCIDLHKYCGGDWQGIID";

		string a="DYFINTDYQVIYQIFTDRFSDGNPANNPTVIYQIVTDRFYQVIYQVTDVIYQIVTDRFADGDRTNNPAGDAFSGDRSNLKLYFGGDWQGIID";
		string b="QWYQVIYQIFTDRFSDGNPANNPTYPTHTSLKKYFGGDWQNNPTGDTDVIYQIVTDRFADGDRTNNPAGDAFSGDRSNLKLYFGGDWQGIID";
		string c="PYRVYQIVTDRFVDGNSANNPTGAAFSSDHSNLKLYFGGDWQGITNTDVIYQIVTDRFADGDRTNNPAGDAFSGDRSNLKLYFGGDWQGIID";*/

		/*string o="TDVIYQIVTDRFADGDRTNNPA";
		string t="NPVIYQIVTDRFSDGNPGNNPS";

		string a="YQVIYQIFTDRFSDGNPANNPT";
		string b="VIVIYQIVTDRFVDGNTSNNPT";
		string c="RFTVYQIVTDRFVDGNSANNPT";*/

		/*string o="ACACACACGCACAAAAACACACCC";
		string t="ACACACACGCACAAAAACACACCC";

		string a="ACACACTGGCACACACACACAAAA";
		string b="ACACACTGGCACACACACACAAAA";
		string c="ACACACTGGCACACACACACACAA";*/

		string o="AAAAAAAAAACCCCCCCCCCCCCC";
		string t="AAAAAAAAAACCCCCCCCCCCCCC";

		string a="ACGGGGGGCACACACACACCCCCAAAAAAAAAACCCCCCCCCCCCCC";
		string b="TACACTGGCACACACACACAAAAAAAAAAAAAACCCCCCCCCCCCCC";
		string c="ACGGGGGGCACACACACACACAAAAAAAAAAAACCCCCCCCCCCCCC";

		vector <string> num(2);
		num[0]=o;
		num[1]=t;

		vector <string> letter(3);
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
		PhyloSupport::AlingMultiSvsMultiS2(num,letter,w1,w2,false,1);
		PhyloSupport::AlingMultiSvsMultiS2(num,letter,w1,w2,false,2);
		PhyloSupport::AlingMultiSvsMultiS2(num,letter,w1,w2,false,3);
		PhyloSupport::AlingMultiSvsMultiS2(num,letter,w1,w2,false,4);
		PhyloSupport::AlingMultiSvsMultiS2(num,letter,w1,w2,false,5);
		PhyloSupport::AlingMultiSvsMultiS2(num,letter,w1,w2,false,6);
		PhyloSupport::AlingMultiSvsMultiS2(num,letter,w1,w2,false,num[0].size()/2);
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
					txt+=seq[i][j+dim*index];
				}
			txt+=" "+PhyloSupport::intToString(j);
			txt+="\n";
			}
		txt+="\n\n";
		index++;
		}


		return txt;

	}


}} // namespace
