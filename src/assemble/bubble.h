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


#ifndef BUBBLE_H
#define BUBBLE_H

#include "types.h"
#include "align_result.h"
#include "consensus_map.h"

struct Match;

struct SNP
{
    vector< pair<Match*, int> > matches_;
    string seq_;
};

struct SNPs
{
    void addSnp( string seq, Match* m, int coord );
    string resolve( string base );
    vector<SNP> snps_;
    int start_, len_;
};

struct Bubble
{
    struct BubbleCoord
    {
        int bubble_, match_, len_;
    };
    struct BubbleMap
    {
        BubbleMap( Match* match, vector< pair<int,int> >& coords, int lanchor, int ranchor );
        void addCoord( int match, int bubble );
        Match* match_;
        vector<BubbleCoord> coords_;
        int anchors_[2];
    };
    Bubble(): type_( 2 ){};
    Bubble( AlignResult& result, ConMap* a, ConMap* b, bool drxn );
    void addMatch( AlignResult& result, ConMap* a, ConMap* b, bool drxn );
    vector<BubbleMap> maps_;
    string template_, consensus_;
    int start_, len_, type_;
};


#endif /* BUBBLE_H */

