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


#ifndef _ClustalW_h_
#define _ClustalW_h_

// Includes:

#include <string>
#include <vector>
#include <iostream>
#include <Alignment.h>


using std::string;
using std::vector;


// Global constants, typedefs, etc. (to avoid):
namespace Victor { namespace Phylo {

    /**@brief  Methods to manages the global statistic data.
     * 
     *@Description  Use for ClustalW
     * */
    class ClustalW {

    public:

        // CONSTRUCTORS/DESTRUCTOR:
    	ClustalW();
    	ClustalW(NewickTree gT);


        ~ClustalW();


    	// PREDICATES:



	    // OPERATORS:

        // HELPERS:
        static string printClustalWFromat(vector <string> seq);
        double scoreClustalW(vector<string> allignSeq);


    protected:
	    // OPERATORS:
        void progressiveAlign();

        // HELPERS:


        // ATTRIBUTES:
        NewickTree guideTree;



    private:
        // HELPERS:

        // ATTRIBUTES:


    };


}}
#endif

