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

#include "match.h"
#include "read.h"
#include <cassert>

Match::Match( Target* tar, MappedRead* read, string t, string r, vector<pair<int,int>> &anchors, int tarCoord )
: tar_( tar ), read_( read ), score_( 0 )
{
    assert( t.size() == r.size() );
    assert( !anchors.empty() );
    anchors_[0] = anchors[0].first;
    anchors_[1] = anchors.back().second;
    for ( int i = 0; i < t.size(); i++ )
    {
        if ( i == anchors[0].first ) anchors_[0] = coords_.size();
        if ( i == anchors.back().second ) anchors_[1] = coords_.size();
        if ( i >= anchors[0].first && i <= anchors.back().second ) score_ += t[i] == r[i] ? 1 : -2;
        if ( r[i] != '-' ) coords_.push_back( tarCoord );
        if ( t[i] != '-' ) tarCoord++;
    }
    assert( coords_.size() == read->seq_.size() );
    read_->matches_.push_back( this );
}

int Match::getAnchor( int drxn )
{
    return coords_[anchors_[drxn]];
}

vector<pair<int, int>> Match::getGaps()
{
    vector<pair<int, int>> gaps;
    for ( int i = 1; i < coords_.size(); i++ ) if ( coords_[i] == coords_[i-1] )
    {
        int gap = 1;
        for ( ; i+1 < coords_.size() && coords_[i] == coords_[i+1]; i++ ) gap++;
        gaps.push_back( make_pair( coords_[i], gap ) );
    }
    
    return gaps;
}