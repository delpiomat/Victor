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


#ifndef _PhyloSupport_h_
#define _PhyloSupport_h_

// Includes:

#include <string>
#include <vector>
#include <iostream>
#include <Alignment.h>

using namespace Victor::Align2;
using std::string;
using std::vector;


// Global constants, typedefs, etc. (to avoid):
namespace Victor { namespace Phylo {

    /**@brief  Methods to manages the global statistic data.
     * 
     *@Description  Use for create Tree during phylogeny
     * */
    class PhyloSupport {
    public:

        // CONSTRUCTORS/DESTRUCTOR:
    	PhyloSupport();


        ~PhyloSupport();


        //STATIC
        //Calc multi Align
        static vector<Alignment> calcAlignmentV(Alignment *aliSec, vector<vector<double> > &distance, bool kimura=false, bool verbose=false);
        //Calc distance from 2 seq of DNA or protein
        static double distanceCalcTwoSeq(string seq1,string seq2);
        static double distanceCalcTwoSeqKimura(string seq1,string seq2);
        static std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems);
        static std::vector<std::string> split(const std::string &s, char delim);
        static double calcDivR(vector<double> vDist);
        static void printMatrix(vector<vector<double> > &distance);



    protected:
    // HELPERS:

    // ATTRIBUTES:



    private:
        // HELPERS:

    };


}}
#endif

