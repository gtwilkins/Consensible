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
    int excess( int d );
    bool match( ConMap* cm, vector<pair<ConMap*, AlignResult>>& hits, int drxn );
    void updateCoords( ConMap* donor, AlignResult& result, int dIndex, bool drxn );
    void updateCoords( SnpAlignResult& result, bool drxn );
    Match* node_;
    int coord_[2]/*of consensus*/, range_[2],/*of read*/ mapped_[2]/*of read, includes bubbles etc.*/;
};


#endif /* CONSENSUS_MAP_H */

