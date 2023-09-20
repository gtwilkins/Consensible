/*
 * Copyright (C) 2017 Glen T. Wilkins <glen.t.wilkins@gmail.com>
 * Written by Glen T. Wilkins
 * 
 * This file is part of the consensible software package <https://github.com/gtwilkins/Consensible>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef ALIGN_RESULT_H
#define ALIGN_RESULT_H

#include "types.h"

struct SNPs;
struct SNP;
struct Bubble;

struct AlignResult
{
    AlignResult();
    int getSeqLen( int i );
    void trimFromEnd( int sIndex, int trimLen, bool drxn );
    string s_[2];
    int start_, len_, lIgnore[2], rIgnore[2], score_;
};

struct SnpAlignResult : AlignResult
{
    struct BubbleAlignCoords
    {
        static void reverse( vector<BubbleAlignCoords>& bubbles, int base );
        bool isDeletion();
        Bubble* bubble_;
        SNPs* snps_;
        SNP* snp_;
        int start_, end_, coord_[2];
        vector<BubbleAlignCoords> bubbles_;
    };
    void reverse();
    vector<BubbleAlignCoords> bubbles_;
//    vector< pair<SNPs*, int> > snps_;
};

#endif /* ALIGN_RESULT_H */

