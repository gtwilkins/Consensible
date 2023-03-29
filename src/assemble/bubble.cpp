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

#include "bubble.h"
#include "match.h"
#include "read.h"
#include <cassert>

void SNPs::addSnp( string seq, Match* m, int coord )
{
    for ( SNP& snp : snps_ ) if ( seq == snp.seq_)
    {
        snp.matches_.push_back( make_pair( m, coord ) );
        return;
    }
    SNP snp;
    snp.seq_ = seq;
    snp.matches_.push_back( make_pair( m, coord ) );
    snps_.push_back( snp );
}

string SNPs::resolve( string base )
{
    string seq = base;
    int best = 0;
    for ( SNP& snp : snps_ ) if ( snp.matches_.size() > best || ( snp.matches_.size() == best && snp.seq_ == base ) )
    {
        best = snp.matches_.size();
        seq = snp.seq_;
    }
    return seq;
}

Bubble::BubbleMap::BubbleMap( Match* match, vector< pair<int,int> >& coords, int lanchor, int ranchor )
: match_( match )
{
    for ( pair<int, int> coord : coords ) addCoord( coord.first, coord.second );
    anchors_[0] = lanchor;
    anchors_[1] = ranchor;
}

void Bubble::BubbleMap::addCoord( int match, int bubble )
{
    if ( !coords_.empty() && coords_.back().bubble_+coords_.back().len_ == bubble && coords_.back().match_+coords_.back().len_ == match )
    {
        coords_.back().len_++;
        return;
    }
    coords_.push_back( Bubble::BubbleCoord() );
    coords_.back().bubble_ = bubble;
    coords_.back().match_ = match;
    coords_.back().len_ = 1;
}

Bubble::Bubble( AlignResult& result, ConMap* a, ConMap* b, bool drxn )
: len_( 0 ), type_( drxn )
{
    int seqLen = result.getSeqLen( 0 );
    start_ = a->node_->coords_[ drxn ? a->node_->size()-seqLen-1 : seqLen ];
    int start = drxn ? a->node_->size()-seqLen : result.lIgnore[0];
    int len = seqLen - ( drxn ? result.rIgnore[0] : result.lIgnore[0] );
    template_ = a->node_->read_->seq_.substr( start, len );
    vector< pair<int,int> > coords[2];
    int counts[3]{ drxn ? a->node_->size()-seqLen : 0, drxn ? b->node_->size()-result.getSeqLen( 1 ) : 0, 0 };
    for ( int i = 0; i < result.start_+result.len_; i++ )
    {
        if ( i >= result.start_ )
        {
            for ( int s : { 0, 1 } ) if ( result.s_[s][i] != '-' ) coords[s].push_back( make_pair( counts[s], counts[2] ) );
            if ( result.s_[0][i] != '-' ) counts[2]++;
        }
        for ( int s : { 0, 1 } ) if ( result.s_[s][i] != '-' ) counts[s]++;
    }
    if ( !coords[0].empty() && coords[0][0].first < a->mapped_[0] ) a->mapped_[0] = coords[0][0].first;
    if ( !coords[0].empty() && coords[0].back().first > a->mapped_[1] ) a->mapped_[1] = coords[0].back().first;
    if ( !coords[1].empty() && coords[1][0].first < b->mapped_[0] ) b->mapped_[0] = coords[1][0].first;
    if ( !coords[1].empty() && coords[1].back().first > b->mapped_[1] ) b->mapped_[1] = coords[1].back().first;
    maps_.push_back( BubbleMap( a->node_, coords[0], 0, template_.size() ) );
    maps_.push_back( BubbleMap( b->node_, coords[1], 0, template_.size() ) );
}

void Bubble::addMatch( AlignResult& result, ConMap* a, ConMap* b, bool drxn )
{
    for ( BubbleMap& bm : maps_ ) if ( bm.match_ == a->node_ )
    {
        vector<pair<int,int>> coords;
        int coord[2]{ drxn ? a->node_->size()-result.getSeqLen( 0 ) : result.lIgnore[0], drxn ? b->node_->size()-result.getSeqLen( 1 ) : result.lIgnore[1] };
        int i = result.start_;
        for ( BubbleCoord bc : bm.coords_ )
        {
            for ( ; i < result.start_+result.len_ && coord[0] < bc.match_+bc.len_; i++ )
            {
                if ( bc.match_ <= coord[0] && result.s_[0][i] != '-' && result.s_[1][i] != '-' )
                {
                    b->mapped_[drxn] = drxn ? max( b->mapped_[1], coord[1] ) : min( b->mapped_[0], coord[1] );
                    coords.push_back( make_pair( coord[1], coord[0]-bc.match_ ) );
                }
                for ( int s : { 0, 1 } ) if ( result.s_[s][i] != '-' ) coord[s]++;
            }
        }
        if ( !coords.empty() ) maps_.push_back( BubbleMap( b->node_, coords, coords[0].second, coords.back().second ) );
        return;
    }
    assert( false );
}
