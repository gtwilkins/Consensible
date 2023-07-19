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

#ifndef CONSENSUS_MAP_H
#define CONSENSUS_MAP_H

#include "types.h"
#include "align_result.h"

struct Match;

struct ConMap
{
    static vector<ConMap*> getFoldableEnds( vector<ConMap*>& maps, vector<pair<int,int>>& remapped, bool d );
    int unmapped( int d );
    int uncoorded( int d );
    bool match( ConMap* cm, vector<pair<ConMap*, AlignResult>>& hits, int drxn );
    void setFinalFoldCoords( int (&base)[2], int (&read)[2], bool drxn );
    void unsetMappedEnd( int cutoff, vector<Bubble*>* bubbles, vector<SNPs*>* snps, bool drxn );
    void updateCoords( ConMap* donor, AlignResult& result, int dIndex, bool drxn );
    void updateCoords( SnpAlignResult& result, bool drxn );
    void updateCoord( int coord, int range );
    void updateMapped( int mapped );
    Match* node_;
    int coord_[2]/*of consensus*/, range_[2],/*of read*/ mapped_[2]/*of read, includes bubbles etc.*/;
};


#endif /* CONSENSUS_MAP_H */

