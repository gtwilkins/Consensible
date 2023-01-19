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

string Consensus::SNPs::resolve( string base )
{
    string seq = base;
    int best = 0;
    for ( SNP& snp : snps_ ) if ( snp.matches_.size() > best || ( snp.matches_.size() == best && snp.seq_ == base ) )
    {
        int best = snp.matches_.size();
        seq = snp.seq_;
    }
    return seq;
}


Consensus::Consensus( vector<Match*>& matches, Target* tar )
: tar_( tar ), template_( tar->seq_ )
{
    coord_[0] = matches[0]->coords_[matches[0]->anchors_[0]];
    coord_[1] = matches[0]->coords_[matches[0]->anchors_[1]];
    for ( Match* m : matches ) addMatch( m );
    sort( maps_.begin(), maps_.end(), []( ConMap& a, ConMap& b ){ return a.coord_[0] == b.coord_[0] ? a.coord_[1] > b.coord_[1] : a.coord_[0] < b.coord_[0]; } );
}

void Consensus::addMatch( Match* m )
{
    ConMap cm;
    cm.node_ = m;
    cm.coord_[0] = m->coords_[m->anchors_[0]];
    cm.coord_[1] = m->coords_[m->anchors_[1]];
    cm.range_[0] = m->anchors_[0];
    cm.range_[1] = m->anchors_[1];
    coord_[0] = min( coord_[0], cm.coord_[0] );
    coord_[1] = max( coord_[1], cm.coord_[1] );
    maps_.push_back( cm );
    
    int prv = m->anchors_[1], cur = m->anchors_[1];
    int minLen = 1;
    for ( int i = m->anchors_[0]; i <= m->anchors_[1]; i++ )
    {
        if ( ( i == m->anchors_[1] || i == cur + minLen ) && prv < cur ) addMismatch( m, prv, cur );
        if ( i >= m->anchors_[1] ) break;
        bool ins = i+1 < m->coords_.size() && m->coords_[i] == m->coords_[i+1];
        bool match = !ins && m->read_->seq_[i] == m->tar_->seq_[m->coords_[i]];
        bool del = i && ( m->coords_[i-1]+1 < m->coords_[i] );
        
        if ( ( !match || del ) && cur + minLen <= i ) prv = i-1;
        if ( !match ) cur = m->anchors_[1];
        if ( match && ( i < cur || del ) ) cur = i;
    }
}

void Consensus::addMismatch( Match* match, int lGood, int rGood )
{
    int matchLen = rGood - lGood - 1;
    int tarLen = lGood < 0 || rGood == match->anchors_[1] ? matchLen : match->coords_[rGood] - match->coords_[lGood] - 1;
    if ( matchLen < 2 && tarLen < 2 )
    {
        int start = match->coords_[lGood+1];
        int i = 0;
        while ( i < snps_.size() && ( snps_[i].start_ == start ? snps_[i].len_ < tarLen : snps_[i].start_ < start ) ) i++;
        if ( i == snps_.size() || snps_[i].start_ != start || snps_[i].len_ != tarLen  )
        {
            snps_.insert( snps_.begin() + i, SNPs() );
            snps_[i].start_ = start;
            snps_[i].len_ = tarLen;
        }
        SNPs::SNP snp;
        snp.seq_ = matchLen ? match->read_->seq_.substr( lGood+1, 1 ) : "";
        int j = 0;
        while ( j < snps_[i].snps_.size() && snps_[i].snps_[j].seq_ != snp.seq_ ) j++;
        if ( j == snps_[i].snps_.size() ) snps_[i].snps_.insert( snps_[i].snps_.begin() + j, snp );
        snps_[i].snps_[j].matches_.push_back( make_pair( match, lGood+1 ) );
    }
    else
    {
        Bubble bubble;
        bubble.start_ = match->coords_[lGood+1];
        bubble.len_ = tarLen;
        bubble.template_ = match->read_->seq_.substr( lGood+1, matchLen );
        bubble_.push_back( bubble );
    }
}

void Consensus::resolveBubbles()
{
    
}

void Consensus::resolveBranches()
{
    vector<ConMap> map( maps_.begin(), maps_.end() );
    for ( int d : { 0, 1 } )
    {
        sort( map.begin(), map.end(), [&]( ConMap& a, ConMap& b ){ return d ? a.coord_[0] < b.coord_[0] : a.coord_[1] < b.coord_[1]; } );
    }
    for ( int d : { 0, 1 } ) for ( Branch& b : branch_[d] )
    {
        
    }
}

void Consensus::resolve()
{
//    resolveBranches();
//    resolveBubbles();
    consensus_ = "";
    int i = 0, j = 0, coord = coord_[0];
    while ( i < bubble_.size() || j < snps_.size() )
    {
        for ( ; i < bubble_.size() && ( j >= snps_.size() || bubble_[i].start_ < snps_[j].start_+snps_[j].len_ ) ; i++ )
        {
            
        }
        
        for ( ; j < snps_.size() && ( i >= bubble_.size() || snps_[j].start_+snps_[j].len_ < bubble_[i].start_ ) ; j++ )
        {
            consensus_ += template_.substr( coord, snps_[j].start_-coord ) + snps_[j].resolve( template_.substr( snps_[j].start_, snps_[j].len_ ) );
            coord = snps_[j].start_ + snps_[j].len_;
        }
    }
    assert( false );
}
