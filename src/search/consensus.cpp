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

#include "consensus.h"
#include "target.h"
#include <cassert>
#include <algorithm>

Consensus::Consensus( vector<Match*>& matches, Target* tar )
: tar_( tar )
{
    coord_[0] = matches[0]->coords_[matches[0]->anchors_[0]];
    coord_[1] = matches[0]->coords_[matches[0]->anchors_[1]];
    for ( Match* gm : matches )
    {
        ConMap cm;
        cm.node_ = gm;
        cm.coord_[0] = gm->coords_[gm->anchors_[0]];
        cm.coord_[1] = gm->coords_[gm->anchors_[1]];
        cm.range_[0] = gm->anchors_[0];
        cm.range_[1] = gm->anchors_[1];
        coord_[0] = min( coord_[0], cm.coord_[0] );
        coord_[1] = max( coord_[1], cm.coord_[1] );
        maps_.push_back( cm );
    }
    sort( maps_.begin(), maps_.end(), []( ConMap& a, ConMap& b ){ return a.coord_[0] == b.coord_[0] ? a.coord_[1] > b.coord_[1] : a.coord_[0] < b.coord_[0]; } );
}
