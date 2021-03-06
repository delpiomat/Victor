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
// Description:     Calculate scores for profile to profile alignment using
//                  euclidean distance.
//
// -----------------x-----------------------------------------------------------

#include <EDistance.h>

namespace Victor { namespace Align2{

    // CONSTRUCTORS:
    /**
     * 
     * @param pro1
     * @param pro2
     */
    EDistance::EDistance(Profile *pro1, Profile *pro2) : ScoringFunction(),
    pro1(pro1), pro2(pro2), offset(1.00) {
    }
    /**
     *  
     * @param pro1
     * @param pro2
     * @param offset
     */
    EDistance::EDistance(Profile *pro1, Profile *pro2, double offset)
    : ScoringFunction(), pro1(pro1), pro2(pro2), offset(offset) {
    }
    /**
     *  
     * @param orig
     */
    EDistance::EDistance(const EDistance &orig) : ScoringFunction(orig) {
        copy(orig);
    }

    EDistance::~EDistance() {
    }


    // OPERATORS:
    /**
     *  
     * @param orig
     * @return 
     */
    EDistance&
            EDistance::operator =(const EDistance &orig) {
        if (&orig != this)
            copy(orig);
        POSTCOND((orig == *this), exception);
        return *this;
    }


    // PREDICATES:
    /**
     *  
     * @param i
     * @param j
     * @return 
     */
    double
    EDistance::scoringSeq(int i, int j) {
        const string residue_indices = "ARNDCQEGHILKMFPSTWYV";

        double freq1 = 0.00;
        double freq2 = 0.00;
        double s = 0.00;

        for (unsigned int k = 0; k < 20; k++) {
            freq1 = pro1->getAminoFrequency(residue_indices[k], (i - 1));
            freq2 = pro2->getAminoFrequency(residue_indices[k], (j - 1));

            s += ((freq1 - freq2) * (freq1 - freq2));
        }

        return offset - sqrt(s);
    }


    // MODIFIERS:
    /**
     *  
     * @param orig
     */
    void
    EDistance::copy(const EDistance &orig) {
        ScoringFunction::copy(orig);
        pro1 = orig.pro1->newCopy();
        pro2 = orig.pro2->newCopy();
        offset = orig.offset;
    }
    /**
     * 
     * @return 
     */
    EDistance*
    EDistance::newCopy() {
        EDistance *tmp = new EDistance(*this);
        return tmp;
    }

}} // namespace
