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

void SNPs::addSnp( vector<SNPs*>& snps, string seq, ConMap* cm, int baseCoord, int matchCoord, int len )
{
    SNPs* cur = NULL;
    for ( SNPs* ss : snps ) if ( ss->start_ == baseCoord && ss->len_ == len ) cur = ss;
    if ( !cur )
    {
        cur = new SNPs();
        cur->start_ = baseCoord;
        cur->len_ = len;
        snps.push_back( cur );
    }
    cur->addSnp( seq, cm->node_, matchCoord );
}

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

void SNPs::addMatch( SnpAlignResult& result, SnpAlignResult::BubbleAlignCoords& bac, ConMap* cm, int& i, int& read )
{
    assert( bac.snp_->seq_.empty() || result.s_[0][i] == result.s_[1][i] );
    for ( pair<Match*, int>& match : bac.snp_->matches_ ) assert( match.first != cm->node_ );
    bac.snp_->matches_.push_back( make_pair( cm->node_, read ) );
    if ( bac.snp_->seq_.size() )
    {
        cm->updateMapped( read++ );
        i++;
    }
}

bool SNPs::isFoldTarget( int (&limits)[2] )
{
    return limits[0] <= start_ && start_+len_ <= limits[1];
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

void SNPs::setRemapped( vector<pair<int,int>>& remapped, unordered_map<SNPs*,int>& oldSnps, vector<SNPs*>& curSnp )
{
    for ( SNPs* s : curSnp )
    {
        auto it = oldSnps.find( s );
        if ( it == oldSnps.end() || it->second < s->snps_.size() ) remapped.push_back( make_pair( s->start_, s->start_+s->len_ ) );
    }
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

Bubble::Bubble( ConMap* cm, bool drxn )
: start_( drxn ? cm->coord_[1]+1 : cm->coord_[0] ), len_( 0 ), type_( drxn )
{
    assert( cm->range_[drxn] == cm->mapped_[drxn] );
    template_ = cm->node_->read_->seq_.substr( drxn ? cm->mapped_[1] : 0, cm->unmapped( drxn ) );
    BubbleMap bm( cm );
    BubbleCoord bc;
    bm.anchors_[0] = bc.bubble_ = 0;
    bm.anchors_[1] = bc.len_ = template_.size();
    bc.match_ = drxn ? cm->mapped_[1] : 0;
    bm.coords_.push_back( bc );
    maps_.push_back( bm );
    cm->mapped_[drxn] = drxn ? cm->node_->size() : 0;
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

Bubble::~Bubble()
{
    for ( IntraBubble* ib : bubbles_ ) delete ib;
}

void Bubble::abort()
{
    for ( BubbleMap& bm : maps_ ) if ( !bm.coords_.empty() )
    {
        if ( bm.coords_[0].match_ <= bm.map_->mapped_[0] )
        {
            bm.map_->mapped_[0] = min( bm.map_->range_[0], bm.coords_.back().match_+bm.coords_.back().len_ );
        }
        if ( bm.coords_.back().match_ + bm.coords_.back().len_ >= bm.map_->mapped_[1] )
        {
            bm.map_->mapped_[1] = max( bm.map_->range_[1], bm.coords_[0].match_ );
        }
    }
    delete this;
}

void Bubble::addBubble( vector<Bubble*>& bubbles, string seq, ConMap* cm, int baseCoord, int matchCoord, int len )
{
    Bubble* bub = NULL;
    for ( Bubble* b : bubbles ) if ( b->start_ == baseCoord && b->len_ == len && b->template_ == seq ) bub = b;
    if ( !bub )
    {
        bub = new Bubble();
        bub->template_ = seq;
        bub->start_ = baseCoord;
        bub->len_ = len;
        bubbles.push_back( bub );
    }
    BubbleMap map( cm );
    BubbleCoord bc( 0, matchCoord, bub->template_.size() );
    map.anchors_[0] = 0;
    map.anchors_[1] = bub->template_.size();
    map.coords_.push_back( bc );
    bub->maps_.push_back( map );
}

void Bubble::addBubble( pair<int,int> bubble, pair<int,int> read, ConMap* cm )
{
    int len = bubble.second-bubble.first;
    string seq = cm->node_->read_->seq_.substr( read.first, read.second-read.first );
    
    Bubble* bub = NULL;
    for ( Bubble* b : bubs_ ) if ( b->start_ == bubble.first && b->len_ == len && b->template_ == seq )
    {
        bub = b;
        assert( false );
    }
    if ( !bub )
    {
        bub = new Bubble();
        bub->template_ = cm->node_->read_->seq_.substr( read.first, read.second-read.first );
        bub->start_ = bubble.first;
        bub->len_ = bubble.second-bubble.first;
        bubs_.push_back( bub );
    }
    BubbleMap map( cm );
    BubbleCoord bc( 0, read.first, bub->template_.size() );
    map.anchors_[0] = 0;
    map.anchors_[1] = bub->template_.size();
    map.coords_.push_back( bc );
    bub->maps_.push_back( map );
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

void Bubble::addMatch( SnpAlignResult& result, SnpAlignResult::BubbleAlignCoords& bac, ConMap* cm, int& i, int& read )
{
    for ( BubbleMap& bm : maps_ ) assert( bm.map_ != cm );
    BubbleMap map( cm );
    BubbleCoord bc( 0, read, 0 );
    map.anchors_[0] = i <= result.start_ ? template_.size() : 0;
    map.anchors_[1] = result.start_+result.len_ < bac.end_ ? 0 : template_.size();
    
    int coord = 0;
    vector<SnpAlignResult::BubbleAlignCoords>::iterator it = bac.bubbles_.begin();
    for ( ; i < min( bac.end_, result.start_+result.len_ ); i++ )
    {
        for ( ; it != bac.bubbles_.end() && i == it->start_; it++ )
        {
            // Update anchors
            if ( result.start_ < it->end_ ) map.anchors_[0] = min( map.anchors_[0], it->start_ );
            if ( result.start_+result.len_ > it->start_ ) map.anchors_[1] = it->end_;
            
            // Terminate and add most recent match run
            if ( bc.len_ ) map.coords_.push_back( bc );
            
            // Add bubble if currently on a mismatch run
            if ( bc.bubble_+bc.len_ < coord || bc.match_ + bc.len_ < read ) addBubble( make_pair( bc.bubble_+bc.len_, coord ), make_pair( bc.match_+bc.len_, read ), cm );
            
            // Ascend bubble tree
            if ( it->bubble_ ) it->bubble_->addMatch( result, *it, cm, i, read );
            if ( it->snps_ ) it->snps_->addMatch( result, *it, cm, i, read );
            
            // Reset match run
            coord = it->coord_[1];
            bc = BubbleCoord( coord, read, 0 );
        }
        if ( i >= min( bac.end_, result.start_+result.len_ ) ) break;
        
        // Update mapped
        if ( i >= result.start_ && i < result.start_ + result.len_ && result.s_[1][i] != '-' ) cm->updateMapped( read );
        
        // Update anchors
        if ( i == result.start_ ) map.anchors_[0] = min( map.anchors_[0], coord );
        if ( i == result.start_+result.len_ ) map.anchors_[1] = coord;
        
        // Start initial run
        if ( i == result.start_ ) bc = BubbleCoord( coord, read, 0 );
        
        if ( i < result.start_ || result.start_ + result.len_ <= i );
        else if ( result.s_[0][i] == result.s_[1][i] )
        {
            if ( bc.bubble_+bc.len_ < coord || bc.match_ + bc.len_ < read )
            {
                if ( bc.len_ ) map.coords_.push_back( bc );
                addBubble( make_pair( bc.bubble_+bc.len_, coord ), make_pair( bc.match_+bc.len_, read ), cm );
                bc = BubbleCoord( coord, read, 0 );
            }
            bc.len_++;
        }
        if ( result.s_[0][i] != '-' ) coord++;
        if ( result.s_[1][i] != '-' ) read++;
    }
    
    // Clean up
    if ( bc.len_ ) map.coords_.push_back( bc );
    if ( bc.bubble_+bc.len_ < coord || bc.match_ + bc.len_ < read ) addBubble( make_pair( bc.bubble_+bc.len_, coord ), make_pair( bc.match_+bc.len_, read ), cm );
    for ( ; i < bac.end_; i++ ) if ( result.s_[1][i] != '-' ) read++;
    
    maps_.push_back( map );
}

bool Bubble::isFoldTarget( ConMap* cm, int (&limits)[2], bool finalised, bool drxn )
{
    int coord[2]{ !type_ ? start_-(int)template_.size() : start_, type_ == 1 ? start_+(int)template_.size() : start_+len_ };
    
    // Must not be entirely outside limits
    if ( limits[1] < coord[0] || coord[1] < limits[0] ) return false;
    
    // Do not allow moving of limits if finalised
    if ( finalised && ( drxn ? coord[0] < limits[0] : limits[1] < coord[1] ) ) return false;
    
    // Do not allow bubble entirely within conmap if limits have not been finalised
    if ( !finalised && ( drxn ? coord[1] <= cm->coord_[1]+1 : cm->coord_[0] <= coord[0] ) ) return false;
    
    // Do not allow moving of limits if too far away from conmap end
    if ( !finalised && ( drxn ? coord[0] <= max( cm->coord_[0], cm->coord_[1]-9 ) : min( cm->coord_[0]+10, cm->coord_[1] ) <= coord[1] ) ) return false;
    
    // Move limits if necessary
    limits[!drxn] = drxn ? min( limits[0], coord[0] ) : max( limits[1], coord[1] );
    
    return true;
}

void Bubble::setRemapped( vector<pair<int,int>>& remapped, unordered_map<Bubble*,int>& oldBubbles, vector<Bubble*>& curBubbles )
{
    for ( Bubble* b : curBubbles )
    {
        auto it = oldBubbles.find( b );
        if ( it == oldBubbles.end() || it->second < b->maps_.size() ) remapped.push_back( make_pair( b->start_, b->start_+b->len_ ) );
    }
}
