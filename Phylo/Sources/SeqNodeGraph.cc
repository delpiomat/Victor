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

#include <NewickTree.h>
#include <PhyloSupport.h>
#include <SeqNodeGraph.h>

#include <sstream>
#include <string>
#include <vector>

// Global constants, typedefs, etc. (to avoid):
using namespace Victor::Align2;
using namespace Victor;
using namespace std;

namespace Victor { namespace Phylo{

	// CONSTRUCTORS:


	SeqNodeGraph::SeqNodeGraph(){
		indexStart=0;
		indexFinish=0;
		totNumSeq=0;
		vector <double> taxV(0);
		taxEdge=taxV;
		tokenSize=0;
		numSeq=0;
		averageTax=0;
		numTokenInSeq2=0;
	}


	//totNumSeq how much of seq in this node minimun 1
	//numSeqInS2 how much seq in the seq2 important for size of vector TaxEge;
	//Weigth depends on guide tree;
	SeqNodeGraph::SeqNodeGraph(unsigned int indexS,unsigned int indexF, unsigned int totNumS,vector <string> tokenS,unsigned int numSeqInS2){
		indexStart=indexS;
		indexFinish=indexF;
		totNumSeq=totNumS;
		vector <double> taxV(numSeqInS2);
		taxEdge=taxV;
		tokenSeq=tokenS;
		tokenSize=tokenSeq[0].size();
		numSeq=tokenSeq.size();
		averageTax=0;
		numTokenInSeq2=numSeqInS2;
	}


	/**
	 *@Description Basic destructor
	 */
	SeqNodeGraph::~SeqNodeGraph() {
		//destroy_tree();//to do
	}

	// PREDICATES:
	unsigned int SeqNodeGraph::returnBestEdge(){
		unsigned int best=0;
		for(unsigned int i=0; i<taxEdge.size();i++){
			//cout<<"taxEdge[best]<taxEdge["<<i<<"] "<<"best="<<best<<" "<<taxEdge[best]<<" < "<<taxEdge[i]<<" "<<(taxEdge[best]<taxEdge[i])<<endl;
			//cout<<"calcolo con zero e i= "<<i<<" tax= "<<taxEdge[i]<<endl;
			if(taxEdge[best]<taxEdge[i]){
				best=i;
			}
		}
		return best;
	}


	int SeqNodeGraph::returnBestEdgeAfterIndex(unsigned index){
		unsigned int best=index+1;
		best=index+1;
		for(unsigned int i=index+1; i<taxEdge.size();i++){
			if(taxEdge[best]<taxEdge[i]){
				//cout<<"taxEdge[best]<taxEdge[i] "<<taxEdge[best]<<" < "<<taxEdge[i]<<endl;
				best=i;
			}
		}
		if(best>=taxEdge.size())
			best=-1;
		return best;
	}

	string SeqNodeGraph::getTokenSeq(unsigned int index){
		return tokenSeq[index];
	}

    int SeqNodeGraph::getIndexStart(){
    	return indexStart;
    }
    int SeqNodeGraph::getIndexFinish(){
    	return indexFinish;
    }
    unsigned int SeqNodeGraph::getTotNumSeq(){
    	return totNumSeq;
    }
    double SeqNodeGraph::getTaxEdgeInPosition(unsigned int i){
    	return taxEdge[i];
    }
    char SeqNodeGraph::getCharOfTokenSeq(unsigned int seqNum, unsigned int pos){
    	return tokenSeq[seqNum][pos];
    }
    unsigned int SeqNodeGraph::getTokenSize(){
    	return tokenSize;
    }
    double SeqNodeGraph::getAverageTax(){
    	return averageTax;
    }




    // OPERATORS:
    void SeqNodeGraph::setIndexStart(unsigned int i){
    	indexStart=i;
    }
    void SeqNodeGraph::setIndexFinish(unsigned int i){
    	indexFinish=i;
    }
    void SeqNodeGraph::setTotNumSeq(int i){
    	numSeq=i;
    }
    void SeqNodeGraph::setTaxEdgeInPosition(unsigned int i, double tax){
    	taxEdge[i]=tax+taxEdge[i];
    }
    void SeqNodeGraph::setTokenSize(int size){
    	tokenSize=size;
    }

    void SeqNodeGraph::printTaxEdge(){
    	cout<<"tax vector"<<endl;
    	for(int i=0;i<taxEdge.size();i++){
    		cout<<taxEdge[i]<<" ";
    	}
    	cout<<endl;
    }

    void SeqNodeGraph::calculateAverageTax(){
    	//cout<<"average"<<endl;
    	averageTax=0;
    	for(unsigned int i=0;i<taxEdge.size();i++){
    		averageTax+=taxEdge[i];
    }
    	averageTax=averageTax/taxEdge.size();
    	//cout<<"end averange "<<averageTax<<endl;
    }

    void SeqNodeGraph::setNode(SeqNodeGraph* node, vector <SeqNodeGraph*> vNode){

    	//cout<<"dentro set node-----------------"<<endl;

    	//matrix config
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
		//end matrix config

    	for(unsigned int i=0; i<vNode.size( );i++)
    	{//how much long vNode, how many node[i] in vNode
    		for(unsigned int j=0; j<vNode[i]->getTotNumSeq() ;j++)
			{//how seq in  node[i] of vNode
    			for(unsigned int x=0; x<vNode[i]->getTokenSize();x++)
    			{//for all char in one string in node[i] of vNode
    				for(unsigned int y=0; y<node->getTotNumSeq();y++)
    				{//for all token of string in node
    					//if(i==0)
    						//cout<<"clacolo Score---------"<<endl;
    					for(unsigned int w=0; w<node->getTokenSize();w++)
    					{//for all char in each seq of node
    						node->setTaxEdgeInPosition(i,sub.score[node->getCharOfTokenSeq(y,w)][vNode[i]->getCharOfTokenSeq(j,x)]);
    						if(i==0){
    							//cout<<" "<<node->getCharOfTokenSeq(y,w)<<" & "<<vNode[i]->getCharOfTokenSeq(j,x)<<endl;
    							//cout<<"\t lo score in posizione i= "<<i<<" vale "<<node->getTaxEdgeInPosition(i)<<endl;
    						}
    					}//end for

    				}//end for
    			}//end for
    		}//end for
    	}//end final for
    	//cout<<"end of all FOR"<<endl;

    	node->calculateAverageTax();
    	//cout<<"near return"<<endl;
    }


    void SeqNodeGraph::setNodeWithWeigth(SeqNodeGraph* node, vector <SeqNodeGraph*> vNode,vector <double> vWeigth1,vector <double> vWeigth2){

        	//cout<<"dentro set node-----------------"<<endl;

        	//matrix config
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
    		//end matrix config

        	for(unsigned int i=0; i<vNode.size( );i++)
        	{//how much long vNode, how many node[i] in vNode
        		for(unsigned int j=0; j<vNode[i]->getTotNumSeq() ;j++)
    			{//how seq in  node[i] of vNode
        			for(unsigned int x=0; x<vNode[i]->getTokenSize();x++)
        			{//for all char in one string in node[i] of vNode
        				for(unsigned int y=0; y<node->getTotNumSeq();y++)
        				{//for all token of string in node
        					//if(i==0)
        						//cout<<"clacolo Score---------"<<endl;
        					for(unsigned int w=0; w<node->getTokenSize();w++)
        					{//for all char in each seq of node
        						node->setTaxEdgeInPosition(i,sub.score[node->getCharOfTokenSeq(y,w)][vNode[i]->getCharOfTokenSeq(j,x)]*vWeigth1[y]*vWeigth2[i]);
        						if(i==0){
        							//cout<<" "<<node->getCharOfTokenSeq(y,w)<<" & "<<vNode[i]->getCharOfTokenSeq(j,x)<<endl;
        							//cout<<"\t lo score in posizione i= "<<i<<" vale "<<node->getTaxEdgeInPosition(i)<<endl;
        						}
        					}//end for

        				}//end for
        			}//end for
        		}//end for
        	}//end final for
        	//cout<<"end of all FOR"<<endl;

        	node->calculateAverageTax();
        	//cout<<"near return"<<endl;
        }








}} // namespace
