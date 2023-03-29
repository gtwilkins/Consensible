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

#ifndef GLIN_ALIGNMENT_H
#define GLIN_ALIGNMENT_H

#include "types.h"
#include "bubble.h"
#include "align_result.h"

struct SNPs;

class Alignment
{
    void score();
    void debug();
    bool isFreeStart( int i, int j );
protected:
    string a_, b_;
    vector< vector<int> > m_, p_;
    int hit_, miss_, gapOpen_, gapExt_, iMax_, jMax_;
    bool freeEnds_[2], scored_;
public:
    Alignment( string& a, string& b, int aStart, int aLen, int bStart, int bLen );
    AlignResult align( int hit, int miss, int gapOpen, int gapExt, bool lAnchored, bool rAnchored );
    AlignResult align( bool lAnchored, bool rAnchored );
    static AlignResult alignBySeed( string& a, string& b, int aStart, int bStart, int len );
};

class SnpAlignment : Alignment
{
    void score();
    vector<SNPs*> snps_;
    vector< vector< pair<int, int> > > snp_;
    int coord_[2];
public:
    SnpAlignment( string& a, string& b, int aStart, int aLen, int bStart, int bLen, vector<SNPs*>& snps );
    SnpAlignResult align( bool lAnchored, bool rAnchored );
};


#endif /* GLIN_ALIGNMENT_H */

