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

#ifndef GLIN_MATCH_H
#define GLIN_MATCH_H

#include "types.h"
#include "align_result.h"

class Target;
struct MappedRead;

struct Match
{
    Match( Target* tar, MappedRead* read, AlignResult& ar, int tarCoord );
//    Match( Target* tar, MappedRead* read, string t, string r, vector<pair<int,int>> &anchors, int tarCoord );
    void anchorCoords( Match* match, vector< pair<int,int>>& coords );
    int getAnchor( int drxn );
    vector<pair<int,int>> getGaps();
    static bool isBridge( Match* l, Match* r, AlignResult& result, int competingCons );
    bool isInsert( int coord );
    int size();
    void updateCoord( int i, int coord );
    void updateCoords( Match* match, vector<int> coords, int start );
    void updateCoords( Match* match, vector< pair<int,int>>& lcoords, vector< pair<int,int>>& rcoords );
    Target* tar_;
    MappedRead* read_;
    vector<int> coords_;
    int anchors_[2], score_;
};


#endif /* GLIN_MATCH_H */

