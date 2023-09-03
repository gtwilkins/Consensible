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

#include "consensus_resolution.h"
#include <algorithm>
#include <cassert>

void ConsensusResolution::branch( vector<Bubble*>& branches, vector<ConMap*>& maps, vector<ConsensusResolution>& resolves, bool drxn )
{
    vector<ConsensusResolution> branched;
    for ( Bubble* b : branches )
    {
        assert( !b->len_ );
        unordered_set<ConMap*> mapped;
        for ( int i = 0; i < maps.size(); i++ ) if ( drxn ? b->start_ <= maps[i]->coord_[1] : maps[i]->coord_[0] < b->start_ ) mapped.insert( maps[i] );
        if ( b->score( mapped ) > 0 ) branched.push_back( b );
    }
    sort( branched.begin(), branched.end(), []( ConsensusResolution& a, ConsensusResolution& b ){ return a.score_ > b.score_; } );
    if ( branched.empty() ) return;
    for ( int i = 0; i < resolves.size(); i++ ) if ( drxn ? branched[0].end_ < resolves[i].end_ : resolves[i].start_ < branched[0].start_ )
    {
        resolves.erase( resolves.begin() + i-- );
    }
    resolves.insert( drxn ? resolves.end() : resolves.begin(), branched[0] );
}

vector<ConsensusResolution> ConsensusResolution::create( vector<Bubble*>& bubs, vector<SNPs*>& snps, vector<ConMap*>& maps )
{
    vector<ConsensusResolution> resolves;
    sort( maps.begin(), maps.end(), []( ConMap* a, ConMap* b ){ return a->coord_[0] == b->coord_[0] ? a->coord_[1] > b->coord_[1] : a->coord_[0] < b->coord_[0]; } );
    int i = 0;
    for ( Bubble* b : bubs )
    {
        unordered_set<ConMap*> mapped;
        while ( i < maps.size() && maps[i]->coord_[1] <= b->start_ ) i++;
        for ( int j = i; j < maps.size() && ( maps[j]->coord_[0] < b->start_ + b->len_ ); j++ ) if ( maps[j]->isConsensus( b->start_, b->start_+b->len_ ) ) mapped.insert( maps[j] );
        if ( b->score( mapped ) > 0 ) resolves.push_back( ConsensusResolution( b ) );
    }
    i = 0;
    for ( SNPs* s : snps )
    {
        unordered_set<Match*> mapped;
        while ( i < maps.size() && maps[i]->coord_[1] <= s->start_ ) i++;
        for ( int j = i; j < maps.size() && ( maps[j]->coord_[0] < s->start_ + s->len_ ); j++ ) if ( maps[j]->isConsensus( s->start_, s->start_+s->len_ ) ) mapped.insert( maps[j]->node_ );
        if ( s->score( mapped ) > 0 ) resolves.push_back( ConsensusResolution( s ) );
    }
    sort( resolves.begin(), resolves.end(), []( ConsensusResolution& a, ConsensusResolution& b ){ return a.start_ == b.start_ ? a.end_ > b.end_ : a.start_ < b.start_; } );
    return resolves;
}

string ConsensusResolution::getConsensus()
{
    if ( bub_ ) return bub_->getConsensus();
    if ( snps_ ) return snps_->getConsensus();
    assert( false );
    return "";
}

bool ConsensusResolution::isClash( ConsensusResolution& rhs )
{
    if ( rhs.start_ == start_ && rhs.end_ == end_ ) return true;
    return start_ < rhs.end_ && rhs.start_ < end_;
}

vector<ConsensusResolution> ConsensusResolution::resolve( vector<Bubble*> branches[2], vector<Bubble*>& bubs, vector<SNPs*>& snps, vector<ConMap*>& maps )
{
    vector<ConsensusResolution> resolves = create( bubs, snps, maps );
    vector<int> doms;
    for ( int i = 0; i < resolves.size(); )
    {
        int len = 1, best = i;
        while ( i+len < resolves.size() && resolves[i].isClash( resolves[i+len] ) ) len++;
        for ( int j = i+1; j < i+len; j++ ) if ( resolves[j].score_ > resolves[best].score_ ) best = j;
        if ( best == i && resolves[i].score_ > 0 ) doms.push_back( i );
        if ( best == i ) i += len;
        else i = best;
    }
    
    for ( int i = 1; i < doms.size(); i++ )
    {
        int best = doms[i];
        for ( int j = doms[i-1]+1; j < doms[i]; j++ ) if ( !resolves[j].isClash( resolves[doms[i-1]] ) && !resolves[j].isClash( resolves[doms[i]] ) )
        {
            if ( best == doms[i] || resolves[j].score_ > resolves[best].score_ ) best = j;
        }
        if ( best != doms[i] && resolves[best].score_ > 0 ) doms.insert( doms.begin() + i--, best );
    }
    vector<ConsensusResolution> result;
    for ( int i : doms ) result.push_back( resolves[i] );
    branch( branches[0], maps, result, 0 );
    branch( branches[1], maps, result, 1 );
    sort( result.begin(), result.end(), []( ConsensusResolution& a, ConsensusResolution& b ){ return a.start_ == b.start_ ? a.end_ > b.end_ : a.start_ < b.start_; } );
    return result;
}
