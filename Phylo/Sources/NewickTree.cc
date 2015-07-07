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

#include <sstream>
#include <string>
#include <vector>

// Global constants, typedefs, etc. (to avoid):
using namespace Victor::Align2;
using namespace Victor;
using namespace std;

namespace Victor { namespace Phylo{

	// CONSTRUCTORS:
	/**
	 *
	 * @param ad
	 * @param gf
	 * @param ss
	 */

std::vector<double> distance;

	NewickTree::NewickTree() {
		root=new iNode;
		root->name="NULL";
		root->left=NULL;
		root->right=NULL;
		root->seq="noPURESeq";
		root->allignSeq="noSeq";
		root->branchLength=-1;
		root->divergenceR=0;
		root->nseq=-1;
		root->isLeaf=true;
	}

	NewickTree::NewickTree(int position,string name,string pureSeq,double divR){
		root=new iNode;
		root->left=NULL;
		root->right=NULL;
		root->allignSeq="noSeq";
		root->branchLength=-1;

		root->name=name;
		root->divergenceR=divR;
		root->seq=pureSeq;
		root->nseq=position;
		root->isLeaf=true;
	}


	NewickTree::NewickTree(int position,string name,string pureSeq){
		root=new iNode;
		root->left=NULL;
		root->right=NULL;
		root->allignSeq="noSeq";
		root->branchLength=-1;

		root->name=name;
		root->divergenceR=0;
		root->seq=pureSeq;
		root->nseq=position;
		root->isLeaf=true;
	}



	/**
	 *@Description Create new Tree. Used for union two different Tree.
	 */
	NewickTree::NewickTree( NewickTree *rTree, NewickTree *lTree) {
		root=new iNode;
		root->allignSeq="noSeq";
		root->branchLength=-1;
		root->isLeaf=false;//internal
		root->name="";
		root->seq="";
		root->left=rTree->getRoot();
		root->right=lTree->getRoot();
		root->divergenceR=0;
	}


	/**
	 *@Description Basic destructor
	 */
	NewickTree::~NewickTree() {
		//destroy_tree();//to do
	}


	// PREDICATES:
	iNode *NewickTree::getRoot(){
		return root;
	}

	iNode *NewickTree::getRightChild(){
		return root->right;
    }

	iNode*NewickTree::getLeftChild(){
		return root->left;
    }

    string NewickTree::getName(){
    	return root->name;
    }

    double NewickTree::getRdiv(){
    	return root->divergenceR;
    }


	// MODIFIERS:



    // OPERATORS:
    /**
     * */
    void NewickTree::neighborJoining(Align2::Alignment ali,bool kimuraD, bool verbose){

		int Nseq=ali.size()+1;//number of seq
		vector<vector<double> > distance(Nseq, std::vector<double>(Nseq, 0));

		cout<<"---------Creating Distance Matrix--------------/"<<endl;
		vector<Alignment>vet=PhyloSupport::calcAlignmentV(&ali,distance,kimuraD,verbose);

		cout<<"---------Print Distance Matrix--------------/"<<endl;
		for(unsigned int i=0; i<distance.size(); i++){
			for(unsigned int j=0; j<distance.size(); j++){
				cout<<distance[i][j]<<" ";
			}
			cout<<endl;
		}

		vector<NewickTree> trees(Nseq);

		//Create all node
		for(unsigned int i=0;i<trees.size();i++){
			if(i==trees.size()-1)
				trees[i] = NewickTree(i,PhyloSupport::split(ali.getTargetName(), ' ')[0],ali.getTarget(),PhyloSupport::calcDivR(distance[distance.size()-1]));

			else
				trees[i] = NewickTree(i,PhyloSupport::split(ali.getTemplateName(i), ' ')[0],ali.getTemplate(i),PhyloSupport::calcDivR(distance[i]));

		}

		//matrix Mij for calculate minimun distance
		vector<vector<double> > Mij(trees.size(), std::vector<double>(trees.size(), 0));

		//memorize the nearest seq
		unsigned int minj=0;
		unsigned int mini=1;
		//vector for new distance
		vector<double> newMediaD ( trees.size());

		//for creat new matrix;
		unsigned int newi=0;
		unsigned int newj=0;

		int count=0;
		while(trees.size()!=3){//start 27,  end 3
			count=1;
			cout<<"----------------WHile num= "<<trees.size()<<"---------------"<<endl;
			//calculate a new distanceMatrix using for each pair of seq
			//Example: M(ij)=d(ij) -[r(i) + r(j)] or in the case of the pair A,B:
			//Example: M(AB)=d(AB) -[(r(A) + r(B)]
			// r=divergenceR, d(AB)=distanceA-B, M(AB)=newDistnace

			for(unsigned int i=0; i<trees.size(); i++){
				for(unsigned int j=0; j<trees.size(); j++){
					if(i!=j){
						Mij[i][j]=distance[i][j] - ( trees[i].getRdiv()+trees[j].getRdiv() );
					}
					else
						Mij[i][j]=0;//diagonal useless
				}
			}

			cout<<endl<<"----------------------------Create Min Distance Matrix--------------------------------"<<endl;
			mini=0;
			minj=1;
			//cout<<"Find min distance"<<endl;
			for(unsigned int i=0; i<trees.size(); i++){
				//cout<<"seq "<<i<<endl;
				for(unsigned int j=0; j<trees.size(); j++){
					//cout<<trees[i].getDistance(j)<<" ";
					if( Mij[mini][minj]  > Mij[i][j] && i!=j)
					{
						minj=j;
						mini=i;
					}
				}
				//cout<<endl;
			}

			cout<<endl<<"----------------------------min i j find! "<<mini<<" "<<minj<<" --------------------------------"<<endl;


			//distance of new node U respect A B the children
			//S(AU) =(d(AB)/2) + ([r(A)-r(B)]/2)
			//S(BU) =(d(AB)/2) + ([r(B)-r(A)]/2)
			trees[mini].setBranchLength((distance[mini][minj]/2) +  ((trees[mini].getRdiv()-trees[minj].getRdiv() )/2) );
			trees[minj].setBranchLength((distance[mini][minj]/2) +  ((trees[minj].getRdiv()-trees[mini].getRdiv() )/2) );


			vector<vector<double> > tmpdistance(trees.size()-1, std::vector<double>(trees.size()-1, 0));
			newi=0;
			newj=0;
			cout<<"startcopy"<<endl;
			//copy old matrix in new matrix, but not mini (A) and minj (B) and not the new node (U)
			for(unsigned int i=0; i<distance.size();i++){
				for(unsigned int j=0; j<distance.size();j++){
					if(i!=mini && j!=minj && i!=minj && j!=mini){
						if(i<mini && i<minj){
							newi=i;
							cout<<"caso 0,";
						}
						else if(i>mini && i<minj){
							newi=i-1;
							cout<<"caso 1,";
						}
						else if(i>mini && i>minj){
							cout<<"caso 2,";
							newi=i-2;
						}
						if(j<mini && j<minj){
							newj=j;
							cout<<"caso 3 ";
						}
						else if(j>mini && j<minj){
							newj=j-1;
							cout<<"caso 4 ";
						}
						else if(j>mini && j>minj){
							newj=j-2;
							cout<<"caso 5 ";
						}
						cout<<"mini e minj "<<mini<<" "<<minj<<" i e j "<<i<<" "<<j<<" newi e new j "<<newi<<" "<<newj<<" "<<"valore "<<distance[i][j]<<endl;
						tmpdistance[newi][newj]=distance[i][j];
					}
				}//end for j
			}//end for i

			//insert the new distance
			/*	example:
				d(CU) = (d(AC) + d(BC) - d(AB) )/2
				d(DU) = (d(AD) + d(BD) - d(AB) )/2
				d(EU) = (d(AE) + d(BE) - d(AB) )/2
				d(FU) = (d(AF) + d(BF) - d(AB) )/2
			*/
			newi=0;
			for(unsigned int i=0; i<distance.size();i++){
				if(i!=mini && i!=minj){
					if(i<mini && i<minj){
						newi=i;
						cout<<"caso 0,";
					}
					else if(i>mini && i<minj){
						newi=i-1;
						cout<<"caso 1,";
					}
					else if(i>mini && i>minj){
						cout<<"caso 2,";
						newi=i-2;
					}
				tmpdistance[newi][tmpdistance.size()-1]=(distance[mini][i]+distance[minj][i]-distance[mini][minj])/2;
				cout<<distance[mini][i]<<" + "<<distance[minj][i]<<" - "<<distance[mini][minj]<<"/2 = "<<tmpdistance[newi][tmpdistance.size()-1]<<endl;
				tmpdistance[tmpdistance.size()-1][newi]=tmpdistance[newi][tmpdistance.size()-1];
				}
			}

			cout<<"3 Difference matrix---------------------------------"<<endl;
			PhyloSupport::printMatrix(distance);
			cout<<"stop first matrix---------------------------------"<<endl;
			PhyloSupport::printMatrix(tmpdistance);
			distance=tmpdistance;
			cout<<"end tmp matrix----------------------------------"<<endl;
			cout<<"THE NEW matrix---------------------------------"<<endl;
			//ccreate ne distance matrix dimension-1
			distance=tmpdistance;
			PhyloSupport::printMatrix(distance);
			cout<<"the distance of tree vector "<<trees.size()<<endl;
			cout<<endl<<"---------------Create distance of other node from new node!--------------------------"<<endl;

			//insert new tree
			trees.push_back(NewickTree(&trees[mini],&trees[minj]));

			//erase trees
			trees.erase(trees.begin()+mini);
			//minus 1, because  length was decreased
			trees.erase(trees.begin()+minj-1);
			cout<<"the distance of tree vector "<<trees.size()<<endl;

			//SET NEW r
			for(unsigned int i=0; i<trees.size();i++)
				trees[i].setRdiv(PhyloSupport::calcDivR(distance[i]));


		}//end while
		NewickTree *tmpT;

		//The root is placed by a "mid-point" method (15) at a position
		//where the means of the branch lengths on either side of the root
		//are equal
		//now use upgma for decide root
		if(distance[0][1]<distance[1][2] || distance[0][2]<distance[1][2]){
			if(distance[0][1]<distance[0][2]){//min distance 0,1
				trees[0].setBranchLength(distance[0][1]/2-trees[0].getValueOFChildBranchLength());
				trees[1].setBranchLength(distance[0][1]/2-trees[1].getValueOFChildBranchLength());
				tmpT=new NewickTree(&trees[0],&trees[1]);
				trees[2].setBranchLength((distance[1][2]/2+distance[0][2]/2)/2-trees[2].getValueOFChildBranchLength());
				tmpT->setBranchLength((distance[1][2]/2+distance[0][2]/2)/2-tmpT->getValueOFChildBranchLength());
				tmpT=new NewickTree(tmpT,&trees[2]);
			}
			else{//min distance 0,2
				trees[0].setBranchLength(distance[0][2]/2-trees[0].getValueOFChildBranchLength());
				trees[2].setBranchLength(distance[0][2]/2-trees[2].getValueOFChildBranchLength());
				tmpT=new NewickTree(&trees[0],&trees[2]);
				trees[1].setBranchLength((distance[1][2]/2+distance[1][0]/2)/2-trees[1].getValueOFChildBranchLength());
				tmpT->setBranchLength((distance[1][2]/2+distance[1][0]/2)/2-tmpT->getValueOFChildBranchLength());
				tmpT=new NewickTree(tmpT,&trees[1]);
			}

		}
		else{//min distance 1,2
			trees[1].setBranchLength(distance[1][2]/2-trees[1].getValueOFChildBranchLength());
			trees[2].setBranchLength(distance[1][2]/2-trees[2].getValueOFChildBranchLength());
			tmpT=new NewickTree(&trees[1],&trees[2]);
			trees[0].setBranchLength((distance[0][2]/2+distance[0][1]/2)/2-trees[0].getValueOFChildBranchLength());
			tmpT->setBranchLength((distance[0][1]/2+distance[0][2]/2)/2-tmpT->getValueOFChildBranchLength());
			tmpT=new NewickTree(tmpT,&trees[0]);
			}
		root=tmpT->getRoot();


	}//end method


    /**
     *
     * */
    void NewickTree::upgma(Align2::Alignment ali,bool kimuraD, bool verbose){

    		int Nseq=ali.size()+1;//number of seq
    		vector<vector<double> > distance(Nseq, std::vector<double>(Nseq, 0));

    		vector<Alignment>vet=PhyloSupport::calcAlignmentV(&ali,distance,kimuraD,verbose);
       		vector<NewickTree> trees(Nseq);


    		if(verbose){
    			cout<<"---------Print Distance Matrix--------------/"<<endl;
    			PhyloSupport::printMatrix(distance);
    		}


    		//Create all node
    		for(unsigned int i=0;i<trees.size();i++){
    			if(i==trees.size()-1)
    				trees[i] = NewickTree(i,PhyloSupport::split(ali.getTargetName(), ' ')[0],ali.getTarget());

    			else
    				trees[i] = NewickTree(i,PhyloSupport::split(ali.getTemplateName(i), ' ')[0],ali.getTemplate(i));

    		}

    		//memorize the nearest seq
    		unsigned int minj=0;
    		unsigned int mini=1;
    		//vector for new distance
    		vector<double> newMediaD ( trees.size());

    		//for creat new matrix;
    		unsigned int newi=0;
    		unsigned int newj=0;


    		while(trees.size()!=1){

    			mini=0;
    			minj=1;
    			//search nearest seq
    			for(unsigned int i=0; i<trees.size(); i++){
    				for(unsigned int j=0; j<trees.size(); j++){
    					if( distance[mini][minj]  > distance[i][j] && i!=j)
    					{
    						minj=j;
    						mini=i;
    					}
    				}
    			}

    			if(verbose)
    				cout<<endl<<"----------------------------min i j find! "<<mini<<" "<<minj<<" --------------------------------"<<endl;


    			//distance of new node U respect A B the children
    			trees[mini].setBranchLength((distance[mini][minj]/2) - trees[mini].getValueOFChildBranchLength());
    			trees[minj].setBranchLength((distance[mini][minj]/2) - trees[mini].getValueOFChildBranchLength());


    			vector<vector<double> > tmpdistance(trees.size()-1, std::vector<double>(trees.size()-1, 0));
    			newi=0;
    			newj=0;

    			//copy old matrix in new matrix, but not mini (A) and minj (B) and not the new node (U)
    			for(unsigned int i=0; i<distance.size();i++){
    				for(unsigned int j=0; j<distance.size();j++){
    					if(i!=mini && j!=minj && i!=minj && j!=mini){
    						if(i<mini && i<minj){
    							newi=i;
    						}
    						else if(i>mini && i<minj){
    							newi=i-1;
    						}
    						else if(i>mini && i>minj){
    							newi=i-2;
    						}
    						if(j<mini && j<minj){
    							newj=j;
    						}
    						else if(j>mini && j<minj){
    							newj=j-1;
    						}
    						else if(j>mini && j>minj){
    							newj=j-2;
    						}
    						if(verbose)
    							cout<<"mini e minj "<<mini<<" "<<minj<<" i e j "<<i<<" "<<j<<" newi e new j "<<newi<<" "<<newj<<" "<<"value "<<distance[i][j]<<endl;
    						tmpdistance[newi][newj]=distance[i][j];
    					}
    				}//end for j
    			}//end for i

    			//insert the new distance
    			/*	example:
    				d(CU) = (d(AC) + d(BC))/2
    				d(DU) = (d(AD) + d(BD))/2
    				d(EU) = (d(AE) + d(BE))/2
    				d(FU) = (d(AF) + d(BF))/2
    			*/
    			newi=0;
    			for(unsigned int i=0; i<distance.size();i++){
    				if(i!=mini && i!=minj){
    					if(i<mini && i<minj){
    						newi=i;
    					}
    					else if(i>mini && i<minj){
    						newi=i-1;
    					}
    					else if(i>mini && i>minj){
    						newi=i-2;
    					}
    				tmpdistance[newi][tmpdistance.size()-1]=(distance[mini][i]+distance[minj][i])/2;
    				if(verbose)
    					cout<<distance[mini][i]<<" + "<<distance[minj][i]<<"/2 = "<<tmpdistance[newi][tmpdistance.size()-1]<<endl;
    				tmpdistance[tmpdistance.size()-1][newi]=tmpdistance[newi][tmpdistance.size()-1];
    				}
    			}


    			distance=tmpdistance;
    			if(verbose){
    				cout<<"The new distance matrix------------------------"<<endl;
    				PhyloSupport::printMatrix(distance);
    				cout<<"the distance of tree vector "<<trees.size()<<endl;
    				cout<<endl<<"---------------Create distance of other node from new node!--------------------------"<<endl;
    			}
    			//insert new tree
    			trees.push_back(NewickTree(&trees[mini],&trees[minj]));

    			//erase trees
    			trees.erase(trees.begin()+mini);
    			//minus 1, because  length was decreased
    			trees.erase(trees.begin()+minj-1);
    			if(verbose)
    				cout<<"the distance of tree vector "<<trees.size()<<endl;


    		}//end while
    		root=trees[0].getRoot();

    		cout<<"------------------end UPGMA-------------------------------"<<endl;
    	}//end method


    // OPERATORS:
    /// Print in a string the tree using Newick format
    string NewickTree::printNewickTree()
    {
    	return printNewickTreeNode(*root)+";";
    }


    string NewickTree::printNewickTreeNode(iNode node)
    {
    	string printTree="";
    	std::string varAsString="";
    	std::ostringstream sstream;
    	if(node.right != NULL)
    		printTree+= "(" +printNewickTreeNode(*node.left);
    	if(node.left != NULL)
    		printTree+= "," + printNewickTreeNode(*node.right)+")";
    	if(node.branchLength<=0)
		{
    		return printTree+node.name;
		}
    	else
    		sstream << node.branchLength;
    		varAsString = sstream.str();
    		return printTree+node.name+":"+varAsString;

    }


    /**
     * Return max values of this node children branchLength.
     * Return zero if left child and right are NULL
     */

    double NewickTree::getValueOFChildBranchLength(){
    	double value1=0;
    	double value2=0;
    	if(root->left!=NULL)
    	{
    		value1=root->left->branchLength;
    	}
    	if(root->right!=NULL)
    	{
    		value2=root->right->branchLength;
    	}
    	if(value1>value2)
    		return value1;
    	return value2;
    }


	/**
	* Set the value of node if is -1 its escape during print
	* @param val
	*
	*/
	void NewickTree::setBranchLength(double val){
		if(val>0)
			root->branchLength=val;
		else
			root->branchLength=0;
	}

	/**
	* Set the value of divergence will use in NJ if is -1 its escape during print
	* @param val
	*
	*/
    void NewickTree::setRdiv(double val){
    	root->divergenceR=val;
    }


}} // namespace
