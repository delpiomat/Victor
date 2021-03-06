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


#ifndef _AMINOACIDHYDROGEN_H_
#define	_AMINOACIDHYDROGEN_H_

#include <AminoAcid.h>
#include <AminoAcidCode.h>
#include <string.h>
#include <map>
#include <list>



namespace Victor { namespace Biopool { 

    /**@brief Implements aminoacid hydrogens.
     * 
     *  Includes methods that allow to add hydrogens to STANDARD aminoacids.
     * */
    class AminoAcidHydrogen {
    public:

        static void loadParam(string inputFile);
        static void setHydrogen(AminoAcid* aa, bool verbose);

    private:

        static map<AminoAcidCode, vector<vector<string> > > paramH;

    };


}} //namespace


#endif	/* _AMINOACIDHYDROGEN_H_ */

