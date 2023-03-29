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

#include "target.h"
#include "local_alignment.h"
#include "consensus.h"
#include "alignment.h"
#include <algorithm>
#include <cassert>
#include <fstream>

bool Target::addMatch( MappedRead* read, int coord )
{
    int readLen = read->seq_.length();
    int startCoord = max( 0, coord-readLen );
    int cutLen = coord-startCoord + ( 2*readLen );
    string t = seq_.substr( startCoord, cutLen ), align[2];
    LocalAlignment la( t, read->seq_, true, false );
    la.realign( align[0], align[1], false );
    vector<pair<int,int>> anchors = LocalAlignment::degapByAnchoring( align[0], align[1], true );
    if ( anchors.empty() ) assert( false );
    if ( anchors.empty() ) return false;
    matches_.push_back( new Match( this, read, align[0], align[1], anchors, startCoord ) );
    
    return true;
}

vector<pair<int, int>> Target::getGaps()
{
    vector<pair<int, int>> gaps;
    for ( Match* gm : matches_ ) for ( pair<int, int> gap : gm->getGaps() )
    {
        int i = gaps.size();
        while ( i > 0 && gaps[i-1].first >= gap.first ) i--;
        if ( i < gaps.size() && gaps[i].first == gap.first ) gaps[i].second = max ( gaps[i].second, gap.second );
        else gaps.insert( gaps.begin()+i, gap );
    }
    return gaps;
}

void Target::sortMatches()
{
    sort( matches_.begin(), matches_.end(), []( Match* a, Match* b ){
        return !a->coords_.empty() && !b->coords_.empty() ? a->coords_[0] < b->coords_[0] : a->coords_.empty() && !b->coords_.empty();
    });
}

vector<Consensus*> Target::assemble()
{
    vector<Consensus*> consensus;
    for ( int i = 0; i < matches_.size(); i++ )
    {
        int anchors[2]{ matches_[i]->getAnchor( 0 ), matches_[i]->getAnchor( 1 ) };
        vector<Match*> overlapping{ matches_[i] };
        for ( ; i+1 < matches_.size() && matches_[i+1]->getAnchor( 0 ) < anchors[1]; i++ )
        {
            anchors[1] = max( anchors[1], matches_[i+1]->getAnchor( 1 ) );
            overlapping.push_back( matches_[i+1] );
        }
        consensus.push_back( new Consensus( overlapping, this ) );
    }
    for ( int i = 0; i+1 < consensus.size(); i++ ) if ( Consensus::bridge( consensus[i], consensus[i+1] ) )
    {
        delete consensus[i+1];
        consensus.erase( consensus.begin() + 1 + i-- );
    }
    return consensus;
}

void Target::print( ofstream& ofs )
{
    sortMatches();
    vector<pair<int, int>> gaps = getGaps();
    int gp = 0;
    ofs << ">" + header_ + "\n";
    for ( int i = 0; i < seq_.size(); i++ )
    {
        if ( gp < gaps.size() && gaps[gp].first == i ) ofs << string( gaps[gp++].second, '-' );
        if ( i < seq_.size() ) ofs << seq_[i];
    }
    ofs << ( gp < gaps.size() ? string( gaps[gp].second, '-' ) : "" ) + "\n";
    for ( Match* gm : matches_ ) if ( !gm->coords_.empty() )
    {
        gp = 0;
        ofs << ">" + to_string( gm->read_->id_ ) + "\n" + string( gm->coords_[0], '-' );
        while ( gp < gaps.size() && gaps[gp].first < gm->coords_[0] ) ofs << string( gaps[gp++].second, '-' );
        for ( int i = 0; i < gm->coords_.size(); i++ )
        {
            int gap = gp < gaps.size() && gaps[gp].first == gm->coords_[i] ? gaps[gp++].second : 0;
            string ins = "";
            while ( i+1 < gm->coords_.size() && gm->coords_[i] == gm->coords_[i+1] ) ins += gm->read_->seq_[i++];
            assert( ins.size() <= gap );
            ofs << ( !gm->coords_[i] ? string( gap-ins.size(), '-' ) : "" ) 
                    + ins 
                    + ( gm->coords_[i] && gm->coords_[i] < seq_.size() ? string( gap-ins.size(), '-' ) : "" ) 
                    + ( i  && gm->coords_[i] - gm->coords_[i-1] > 1 ? string( gm->coords_[i] - gm->coords_[i-1] - 1, '-' ) : "" )
                    + gm->read_->seq_[i];
        }
        ofs << "\n";
    }
}