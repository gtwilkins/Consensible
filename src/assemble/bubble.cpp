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
#include "alignment.h"
#include <cassert>
#include <algorithm>
#include <limits>
#include <iostream>

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
    for ( SNP& snp : snps_ ) if ( seq == snp.seq_ )
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
    bool finished = read == cm->node_->size();
    assert( bac.snp_->seq_.empty() || result.s_[0][i] == result.s_[1][i] || !len_ );
    int len = bac.snp_->seq_.size();
    if ( !bac.snp_->seq_.empty() && result.s_[0][i] != result.s_[1][i] )
    {
        string seq = result.s_[1].substr( i, bac.end_-bac.start_ );
        len = seq.size();
        assert( seq.size() == 1 && seq != "-" );
        if ( !finished ) addSnp( seq, cm->node_, read );
    }
    else if ( !finished )
    {
        for ( pair<Match*, int>& match : bac.snp_->matches_ ) assert( match.first != cm->node_ );
        bac.snp_->matches_.push_back( make_pair( cm->node_, read ) );
    }
    if ( len )
    {
        cm->updateMapped( read++ );
        i++;
    }
}

string SNPs::getConsensus()
{
    assert( !snps_.empty() || score_ > 0 );
    return  snps_[0].seq_;
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

int SNPs::score( unordered_set<Match*>& matches )
{
    for ( SNP& snp : snps_ ) for ( pair<Match*, int>& match: snp.matches_ )
    {
        assert( matches.find( match.first ) == matches.end() );
    }
    sort( snps_.begin(), snps_.end(), []( SNP& a, SNP& b ){ return a.matches_.size() > b.matches_.size(); } );
    return score_ = ( snps_.empty() ? 0 : snps_[0].matches_.size() ) - matches.size();
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
    anchors_[0] = min( anchors_[0], bubble );
    anchors_[1] = max( anchors_[1], bubble+1 );
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

void Bubble::BubbleMap::addCoord( int match, int bubble, int len )
{
    int i = 0;
    for ( ; i < coords_.size(); i++ )
    {
        if ( coords_[i].bubble_ == bubble+len && coords_[i].match_ == match+len )
        {
            len += coords_[i].len_;
            coords_.erase( coords_.begin() + i );
            addCoord( match, bubble, len );
            return;
        }
        if ( coords_[i].bubble_ + coords_[i].len_ == bubble && coords_[i].match_ + coords_[i].len_ == match )
        {
            coords_[i].len_ += len;
            return;
        }
        if ( bubble < coords_[i].bubble_ ) break;
    }
    coords_.insert( coords_.begin() + i, BubbleCoord( bubble, match, len ) );
}

Bubble::BubbleMapDetails::BubbleMapDetails( Bubble* bub, BubbleMap& map )
: map_( map.map_ )
{
    start_ = map.map_->node_->size();
    end_ = 0;
    bub->getMapRange( map.map_, start_, end_ );
    assert( start_ < end_ );
//    assert( map.coords_.size() );
//    BubbleCoord* coords[2] = { &map.coords_[0], &map.coords_.back() };
//    start_ = coords[0]->match_;
//    end_ = coords[1]->match_ + coords[1]->len_;
    off_ = 0;
    bool complete[2]{ !map.anchors_[0], map.anchors_[1] == bub->template_.size() };
//    if ( bub->type_ == 2 ) assert( complete[0] && complete[1] );
    if ( bub->type_ == 1 ) assert( complete[0] );
    if ( bub->type_ == 0 ) assert( complete[1] );
    if ( complete[0] && complete[1] ) type_ = 2;
    else if ( complete[0] ) type_ = 1;
    else if ( complete[1] ) type_ = 0;
    else assert( false );
    seq_ = map.map_->node_->read_->seq_.substr( start_, end_-start_ );
}

int Bubble::BubbleMapDetails::setEnd( Bubble* bub, BubbleMap& map )
{
    assert( map.coords_.size() );
    start_ = min( start_, map.coords_[0].match_ );
    end_ = max( end_, map.coords_.back().match_+map.coords_.back().len_ );
    for ( Bubble* b : bub->bubs_ ) for ( BubbleMap& bm : b->maps_ ) if ( bm.map_ == map.map_ ) setEnd( b, bm );
}

bool Bubble::BubbleAltPath::add( BubbleMapDetails& map )
{
    bool added = map.seq_ == seq_;
    if ( !added && map.type_ != 2 )
    {
        size_t l = seq_.find( map.seq_ ), r = seq_.rfind( map.seq_ );
        if ( map.type_ == 1 && type_ != 0 && l != string::npos && l == 0 )
        {
            added = true;
        }
        if ( map.type_ == 0 && type_ != 1 && r != string::npos && r == seq_.size()-map.seq_.size() )
        {
            added = true;
        }
    }
    if ( !added ) return false;
    map.off_ = !map.type_ ? seq_.size()-map.seq_.size() : 0;
    map_.push_back( map );
    return true;
}

void Bubble::BubbleAltPath::assign( int* coords )
{
    cut_[0] = cut_[1] = -1;
    bool asc  = true;
    for ( int i = 0; i < path_.size(); i++ )
    {
        if ( path_[i].second.first <= coords[0] && coords[0] <= path_[i].second.second && ( !i || !asc || path_[i].second.first < coords[0] ) )
        {
            cut_[0] = i;
            coord_[0] = asc ? coords[0]-path_[i].second.first : path_[i].first->template_.size()-( path_[i].second.second - coords[0] );
        }
        if ( i == path_.size()/2 ) asc = false;
        if ( path_[i].second.first <= coords[1] && coords[1] <= path_[i].second.second && ( i+1 == path_.size() || asc || coords[1] < path_[i].second.second ) )
        {
            cut_[1] = i;
            coord_[1] = asc ? coords[1]-path_[i].second.first : path_[i].first->template_.size()-( path_[i].second.second - coords[1] );
        }
    }
    assert( cut_[0] >= 0 && cut_[1] >= 0 );
}

bool Bubble::BubbleAltPath::bridge( Bubble* b, BubbleAltPath& rhs, vector<BubbleAltPath>& paths )
{
    for ( int i = 0; i < map_.size(); i++ ) for ( int j = 0; j < rhs.map_.size(); j++ )
    {
        BubbleMapDetails* l = &map_[i],* r = &rhs.map_[j];
        int off;
        if ( !b->getOffset( l->map_, r->map_, off ) ) continue;
        off -= l->start_;
        int len = 0, miss = 0;
        for ( int i = 0; i+off < l->seq_.size() && i < r->seq_.size(); i++ )
        {
            if ( l->seq_[i+off] != r->seq_[i] ) miss++;
            len++;
        }
        assert( len && !miss );
        BubbleAltPath path( *l );
        path.seq_ += r->seq_.substr( len );
        path.type_ = 2;
        bool dupe = false;
        for ( BubbleAltPath& alt : paths ) if ( alt.seq_ == path.seq_ ) dupe = true;
        if ( dupe ) continue;
        while ( ++i < map_.size() ) assert( path.add( map_[i] ) );
        while ( j < rhs.map_.size() ) assert( path.add( rhs.map_[j++] ) );
        paths.push_back( path );
        return true;
    }
    return false;
}

vector<Bubble::BubbleAltPath> Bubble::BubbleAltPath::build( Bubble* b )
{
    vector<BubbleMapDetails> maps;
    for ( BubbleMap& bm : b->maps_ ) maps.push_back( BubbleMapDetails( b, bm ) );
    sort( maps.begin(), maps.end(), []( BubbleMapDetails& a, BubbleMapDetails& b ){ return a.seq_.size() > b.seq_.size(); } );
    vector<BubbleAltPath> paths;
    vector<pair<BubbleAltPath, bool>> branch[2];
    for ( BubbleMapDetails& bmd : maps )
    {
        bool added = false;
        for ( BubbleAltPath& path : paths ) if ( path.add( bmd ) ) added = true;
        if ( !added ) paths.push_back( BubbleAltPath( bmd ) );
    }
    for ( int i = 0; i < paths.size(); i++ ) if ( paths[i].type_ != 2 )
    {
        branch[paths[i].type_].push_back( make_pair( paths[i], false ) );
        paths.erase( paths.begin() + i-- );
    }
    
    for ( pair<BubbleAltPath, bool>& l : branch[1] ) for ( pair<BubbleAltPath, bool>& r : branch[0] )
    {
        if ( r.second && l.second ) continue;
        if ( l.first.bridge( b, r.first, paths ) ) l.second = r.second = true;
    }
    return paths;
}

vector<Bubble*> Bubble::BubbleAltPath::getStack( int d )
{
    vector<Bubble*> stack;
    int i = cut_[0], mid = path_.size()/2;
    if ( d != 2 ) i = cut_[d];
    else if ( mid <= cut_[0] ) i = cut_[0];
    else if ( cut_[1] <= mid ) i = cut_[1];
    else i = mid;
    bool asc = i <= mid;
    for ( int j = asc ? 0 : path_.size(); asc ? j <= i : --j >= i; j+=( asc ? 1 : 0 ) ) stack.push_back( path_[j].first );
    return stack;
}

void Bubble::BubbleAltPath::remap( SnpAlignResult::BubbleAlignCoords& bac, int& coord )
{
    assert( bac.bubble_ || bac.end_-bac.start_ == 1 );
    if ( bac.snps_ )
    {
        string seq = result_.s_[1][bac.start_] == '-' ? "" : string( 1, result_.s_[1][bac.start_] );
        for ( BubbleMapDetails& map : map_ ) if ( map.off_ <= coord && coord+seq.size() <= map.seq_.size()+map.off_ ) bac.snps_->addSnp( seq, map.map_->node_, map.start_+coord-map.off_ );
        assert( bac.bubbles_.empty() );
        coord += seq.size();
        return;
    }
    if ( bac.bubble_->template_.empty() )
    {
        for ( BubbleMapDetails& map : map_ ) if ( map.off_ <= coord && coord <= map.seq_.size()+map.off_ )
        {
            if ( !map.type_ && ( coord-map.off_ ) == 0 ) continue;
            if ( map.type_ == 1 && ( coord-map.off_ ) == map.seq_.size() ) continue;
            BubbleMap bm( map.map_ );
            bm.anchors_[0] = bm.anchors_[1] = 0;
            bac.bubble_->maps_.push_back( bm );
        }
        return;
    }
    vector<pair<BubbleMap, BubbleMapDetails*>> maps;
    for ( BubbleMapDetails& map : map_ )
    {
        BubbleMap bm( map.map_ );
        bm.anchors_[0] = bac.bubble_->template_.size();
        bm.anchors_[1] = 0;
        maps.push_back( make_pair( bm, &map ) );
    }
    
    int b = 0, j = 0, nxt = bac.bubbles_.empty() ? bac.end_ : bac.bubbles_[0].start_, coordStart = coord;
    for ( int i = bac.start_; i < bac.end_; )
    {
        if ( i == nxt )
        {
            for ( pair<BubbleMap, BubbleMapDetails*>& map : maps ) if ( coord < map.second->seq_.size() ) map.first.anchors_[1] = max( map.first.anchors_[1], bac.bubbles_[j].bubble_->start_+bac.bubbles_[j].bubble_->len_ );
            remap( bac.bubbles_[j], coord );
            for ( pair<BubbleMap, BubbleMapDetails*>& map : maps ) if ( map.second->off_ < coord ) map.first.anchors_[0] = min( map.first.anchors_[0], bac.bubbles_[j].bubble_->start_ );
            b = bac.bubbles_[j].bubble_->start_ + bac.bubbles_[j].bubble_->len_;
            i = bac.bubbles_[j++].end_;
            nxt = j < bac.bubbles_.size() ? bac.bubbles_[j].start_ : bac.end_;
        }
        else if ( result_.s_[0][i] != result_.s_[1][i] )
        {
            int start = b, len = 0;
            string seq = "";
            for ( ; i < nxt && result_.s_[0][i] != result_.s_[1][i]; i++ )
            {
                if ( result_.s_[0][i] != '-' ) len++;
                if ( result_.s_[1][i] != '-' ) seq += result_.s_[1][i];
            }
            Bubble* bub = Bubble::addBubble( seq, bac.bubble_->bubs_, start, len );
            for ( pair<BubbleMap, BubbleMapDetails*>& map : maps ) if ( map.second->off_ < coord+seq.size() && coord < map.second->seq_.size()+map.second->off_ )
            {
                BubbleMap bm( map.second->map_ );
                bm.anchors_[0] = max( 0, map.second->off_-coord );
                bm.anchors_[1] = min( seq.size(), map.second->seq_.size()+map.second->off_-coord );
                assert( ( seq.empty() ? bm.anchors_[0] == bm.anchors_[1] : bm.anchors_[0] < bm.anchors_[1] ) && bm.anchors_[1]-bm.anchors_[0] <= seq.size() );
                map.first.anchors_[0] = min( bub->start_, map.first.anchors_[0] );
                map.first.anchors_[1] = max( bub->start_ + bub->len_, map.first.anchors_[1] );
                BubbleCoord bc( bm.anchors_[0], map.second->start_ + max( 0, coord-map.second->off_ ), bm.anchors_[1]-bm.anchors_[0] );
                assert( bub->template_.substr( bc.bubble_, bc.len_ ) == map.second->map_->node_->read_->seq_.substr( bc.match_, bc.len_ ) );
                if ( bc.len_ ) bm.coords_.push_back( bc );
                bub->maps_.push_back( bm );
            }
            coord += seq.size();
            b += len;
        }
        else
        {
            assert( result_.s_[0][i] != '-' );
            for ( pair<BubbleMap, BubbleMapDetails*>& map : maps ) if ( map.second->off_ <= coord && coord < map.second->seq_.size()+map.second->off_ )
            {
                map.first.addCoord( map.second->start_+coord-map.second->off_, b );
                map.first.anchors_[0] = min( map.first.anchors_[0], b );
                map.first.anchors_[1] = max( map.first.anchors_[1], b+1 );
            }
            
            coord++;
            b++;
            i++;
        }
    }
    bool added = false;
    for ( pair<BubbleMap, BubbleMapDetails*>& map : maps ) if ( map.second->off_ < coord && coordStart < map.second->seq_.size()+map.second->off_ )
    {
        assert( map.first.anchors_[0] < map.first.anchors_[1] || ( bac.bubble_->template_.empty() && map.first.anchors_[0] == map.first.anchors_[1] ) );
        bac.bubble_->maps_.push_back( map.first );
        added = true;
    }
    assert( added );
}

bool Bubble::BubbleAltPath::set( SnpAlignResult result )
{
    score_ = -1;
    result_ = result;
    int type = map_[0].type_, range = 0;
    for ( SnpAlignResult::BubbleAlignCoords& bac : result_.bubbles_ )
    {
        int scores[3] = { 0, 0, min( 2, (int)( bac.bubble_ ? bac.bubble_->template_ : bac.snp_->seq_ ).size() ) };
        vector<int> flankGaps[2];
        for ( int i = 0; i < result_.s_[0].size(); i++ )
        {
            int gap = 0;
            for ( int s : { 0, 1 } ) if ( result_.s_[s][i] == '-' ) 
            {
                for ( ; i < result_.s_[0].size() && result_.s_[s][i] == '-'; i++ ) gap++;
                if ( !s ) range += gap;
                break;
            }
            if ( gap )
            {
                int j = i <= bac.start_ ? 0 : ( i-gap < bac.end_ ? 2 : 1 );
                if ( i > bac.start_ && i-gap < bac.end_ ) j = 2;
                if ( ( type != 0 || range ) && ( type != 1 || range < seq_.size() ) ) scores[j] -= 3 + gap;
                if ( j != 2 ) flankGaps[j].push_back( gap );
                i--;
            }
            else
            {
                int j = i < bac.start_ ? 0 : ( i < bac.end_ ? 2 : 1 );
                scores[j] += result_.s_[0][i] == result_.s_[1][i] ? 1 : -1;
                range++;
            }
        }
        for ( int d : { 0, 1 } ) for ( int j = 1; j < flankGaps[d].size(); j++ ) scores[!d] -= 3 + flankGaps[d][j];
        int score = scores[2] + min( 0, max( scores[0], scores[1] ) );
        if ( score <= score_ ) continue;
        score_ = score;
        bestAlign_ = bac;
    }
    return score_ >= 0;
}

bool Bubble::BubbleAltPath::compatible( Bubble* b, SnpAlignResult::BubbleAlignCoords& bac )
{
    if ( bac.bubble_ && bac.bubble_->start_ < b->start_ ) return false;
    if ( bac.bubble_ && b->start_ + b->len_ < bac.bubble_->start_ + bac.bubble_->len_ ) return false;
    if ( bac.snps_ && bac.snps_->start_ < b->start_ ) return false;
    if ( bac.snps_ && b->start_ + b->len_ < bac.snps_->start_ + bac.snps_->len_ ) return false;
    return true;
}

bool Bubble::BubbleAltPath::presplit( Bubble* base, SnpAlignResult::BubbleAlignCoords& bac )
{
    if ( bac.bubble_ && bac.bubble_->start_ < base->start_ )
    {
        return false;
        int i = 0;
        for ( ; i < result_.s_[0].size() && result_.s_[1][i] == '-'; i++ ) assert( result_.s_[0][i] != '-');
        assert( !bac.start_ );
        int coord[2]{ i-1, i };
        vector<Bubble*> stacks[2]{ { bac.bubble_ }, { bac.bubble_ } };
        SnpAlignResult::BubbleAlignCoords* pBac = &bac;
        for ( int d : { 0, 1 } ) for ( int added = 1; added-- > 0; ) for ( SnpAlignResult::BubbleAlignCoords& sbac : pBac->bubbles_ ) if ( sbac.start_ <= coord[d] && coord[d] < sbac.end_ )
        {
            coord[d] = coord[d] - sbac.start_;
            stacks[d].push_back( sbac.bubble_ );
            pBac = &sbac;
            added = 1;
            assert( false );
            break;
        }
        
        bac.bubble_->split( stacks, coord, base->start_ );
        assert( false );
        return true;
    }
    if ( bac.bubble_ && base->start_ + base->len_ < bac.bubble_->start_ + bac.bubble_->len_ )
    {
        return false;
        assert( false );
        return true;
    }
    return false;
}

vector<Bubble*> Bubble::BubbleAltPath::split( Bubble* base, SnpAlignResult::BubbleAlignCoords& bac )
{
    vector<Bubble*> created;
    int len[2]{0};
    for ( int i = 0; i < result_.s_[0].size(); i++ )
    {
        if ( i < bac.start_ && result_.s_[1][i] != '-' ) len[0]++;
        if ( i < bac.end_ && result_.s_[1][i] != '-' ) len[1]++;
    }
    if ( bac.start_ && ( map_[0].type_ != 0 || len[0] ) )
    {
        created.push_back( new Bubble( base, *this, bac, len[0], 0 ) );
    }
    if ( bac.end_ < result_.s_[0].size() && ( map_[0].type_ != 1 || len[1] ) )
    {
        created.push_back( new Bubble( base, *this, bac, len[1], 1 ) );
    }
//    if ( bac.coord_[0] )
//    {
//        created.push_back( new Bubble( base, *this, bac, bac.coord_[0], 0 ) );
//    }
//    if ( bac.coord_[1] < seq_.size() ) created.push_back( new Bubble( base, *this, bac, bac.coord_[1], 1 ) );
    
    int coord = 0;
    for ( int i = 0; i < bac.start_; i++ ) if ( result_.s_[1][i] != '-' ) coord++;
    remap( bac, coord );
    
    for ( BubbleMapDetails& map : this->map_ ) base->remove( map.map_ );
    return created;
}

Bubble::Bubble( Bubble* base, BubbleAltPath& path, SnpAlignResult::BubbleAlignCoords& bac, int cut, bool drxn )
: template_( path.seq_.substr( drxn ? cut : 0, drxn ? path.seq_.size()-cut : cut ) ), type_( base->type_ != 2 && base->type_ == drxn ? base->type_ : 2 ), score_( 0 )
{
    int bacstart = bac.bubble_ ? bac.bubble_->start_ : bac.snps_->start_;
    int baclen = bac.bubble_ ? bac.bubble_->len_ : bac.snps_->len_;
    start_ = drxn ? bacstart+baclen : base->start_;
    len_ = drxn ? base->start_+base->len_-start_ : bacstart-start_;
    int coord[2]{ drxn ? cut : 0, drxn ? (int)path.seq_.size() : cut };
    int cutoff[2]{ coord[0] - int( drxn && template_.empty() ), coord[1] + int( !drxn && template_.empty() ) };
    for ( BubbleMapDetails& map : path.map_ ) if ( map.off_ < cutoff[1] && cutoff[0] < map.seq_.size()+map.off_ )
    {
        BubbleMap bm( map.map_ );
        bm.anchors_[0] = template_.size();
        bm.anchors_[1] = 0;
        for ( int i = max( map.off_, coord[0] ); i < min( map.off_+(int)map.seq_.size(), coord[1] ); i++ ) assert( map.map_->node_->read_->seq_[map.start_+i-map.off_] == template_[i-coord[0]] );
        for ( int i = max( map.off_, coord[0] ); i < min( map.off_+(int)map.seq_.size(), coord[1] ); i++ ) bm.addCoord( map.start_+i-map.off_, i-coord[0] );
        maps_.push_back( bm );
    }
}

Bubble::Bubble( ConMap* cm, bool drxn )
: start_( drxn ? cm->coord_[1]+1 : cm->coord_[0] ), len_( 0 ), type_( drxn ), score_( 0 )
{
    assert( cm->range_[drxn] == cm->mapped_[drxn] );
    template_ = cm->node_->read_->seq_.substr( drxn ? cm->mapped_[1]+1 : 0, cm->unmapped( drxn ) );
    BubbleMap bm( cm );
    BubbleCoord bc;
    bm.anchors_[0] = bc.bubble_ = 0;
    bm.anchors_[1] = bc.len_ = template_.size();
    bc.match_ = drxn ? cm->mapped_[1]+1 : 0;
    bm.coords_.push_back( bc );
    maps_.push_back( bm );
    cm->mapped_[drxn] = drxn ? cm->node_->size() : 0;
}

Bubble::Bubble( AlignResult& result, ConMap* a, ConMap* b, bool drxn )
: len_( 0 ), type_( drxn ), score_( 0 )
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

Bubble::Bubble( AlignResult& result, ConMap* l, ConMap* r )
: type_( 2 ), score_( 0 )
{
    // This function creates a bridging bubble from the overlap of two reads
    int count[2]{0}, anchor[2]{ 0, r->node_->size() };
    for ( int i = 0; i < result.s_[0].size() && anchor[1] == r->node_->size(); i++ )
    {
        bool anchored[2]{ l->node_->anchors_[0]+(result.s_[0][i] != '-'?0:1) <= count[0] && count[0] <= l->node_->anchors_[1]
                        , r->node_->anchors_[0]+(result.s_[1][i] != '-'?0:1) <= count[1] && count[1] <= r->node_->anchors_[1] };
//        if ( count[0] <= l->node_->anchors_[1] && count[1] < r->node_->anchors_[0] && result.s_[0][i] != '-' ) anchor[0] = count[0];
//        if ( r->node_->anchors_[0] <= count[1] && l->node_->anchors_[1] < count[0] && result.s_[1][i] != '-' ) anchor[1] = count[1];
        if ( anchored[0] && !anchored[1] && result.s_[0][i] != '-' ) anchor[0] = count[0];
        if ( anchored[1] && !anchored[0] && result.s_[1][i] != '-' ) anchor[1] = count[1];
        for ( int d : { 0, 1 } ) if ( result.s_[d][i] != '-' ) count[d]++;
    }
    assert( anchor[0] && anchor[1] < r->node_->size() ); // This would mean that no anchors were found
    count[0] = count[1] = 0;
    start_ = l->node_->coords_[anchor[0]] + 1;
    len_ = r->node_->coords_[anchor[1]] - start_;
    
    int pSeq = 0, excess[2]{ (int)l->node_->size()-1-anchor[0]-result.rIgnore[0], -result.lIgnore[1] };
    for ( int i = 0; i < result.s_[0].size() && ( count[1] < anchor[1] || result.s_[1][i] == '-' ); i++ )
    {
        if ( !pSeq && ( ( result.s_[0] == result.s_[1] && excess[0] <= excess[1] ) || excess[0] <= 0 ) ) pSeq = 1;
        
        if ( anchor[0] < count[0] && result.s_[pSeq][i] != '-' ) template_ += result.s_[pSeq][i];
        
        if ( anchor[0] < count[0] && result.s_[0][i] != '-' ) excess[0]--;
        if ( count[1] < anchor[1] && result.s_[1][i] != '-' ) excess[1]++;
        for ( int d : { 0, 1 } ) if ( result.s_[d][i] != '-' ) count[d]++;
    }
}

Bubble::~Bubble()
{
    for ( Bubble* b : bubs_ ) delete b;
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
    Bubble* bub = addBubble( seq, bubbles, baseCoord, len );
    BubbleMap map( cm );
    BubbleCoord bc( 0, matchCoord, bub->template_.size() );
    map.anchors_[0] = 0;
    map.anchors_[1] = bub->template_.size();
    if ( bc.len_ ) map.coords_.push_back( bc );
    bub->maps_.push_back( map );
}

void Bubble::addBubble( pair<int,int> bubble, pair<int,int> read, ConMap* cm )
{
    assert( read.first < cm->node_->size() );
    int len = bubble.second-bubble.first;
    string seq = cm->node_->read_->seq_.substr( read.first, read.second-read.first );
    
    Bubble* bub = addBubble( cm->node_->read_->seq_.substr( read.first, read.second-read.first ), bubs_, bubble.first, len );
    BubbleMap map( cm );
    BubbleCoord bc( 0, read.first, bub->template_.size() );
    map.anchors_[0] = 0;
    map.anchors_[1] = bub->template_.size();
    if ( bc.len_ ) map.coords_.push_back( bc );
    bub->maps_.push_back( map );
}

Bubble* Bubble::addBubble( string seq, vector<Bubble*>& bubbles, int start, int len )
{
    for ( Bubble* b : bubbles ) if ( b->start_ == start && b->len_ == len && b->template_ == seq )
    {
        assert( false );
        return b;
    }
    Bubble* bub = new Bubble();
    bub->template_ = seq;
    bub->start_ = start;
    bub->len_ = len;
    bubbles.push_back( bub );
    return bub;
}

void Bubble::addMatch( ConMap* cm, int start )
{
    BubbleMap map( cm );
    BubbleCoord bc( 0, start, template_.size() );
    map.anchors_[0] = 0;
    map.anchors_[1] = template_.size();
    map.coords_.push_back( bc );
    maps_.push_back( map );
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
    bool finished = read == cm->node_->size();
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
//            if ( result.start_ < it->end_ ) map.anchors_[0] = min( map.anchors_[0], it->start_ );
//            if ( result.start_+result.len_ > it->start_ ) map.anchors_[1] = it->end_;
            if ( result.start_ < it->end_ ) map.anchors_[0] = min( map.anchors_[0], it->bubble_->start_ );
            if ( result.start_+result.len_ > it->start_ ) map.anchors_[1] = max( map.anchors_[1], it->bubble_->start_+it->bubble_->len_ );
            
            // Terminate and add most recent match run
            if ( bc.len_ ) map.coords_.push_back( bc );
            
            // Add bubble if currently on a mismatch run
            if ( i <= result.start_ || result.start_+ result.len_ < i );
            else if ( bc.bubble_+bc.len_ < coord || bc.match_ + bc.len_ < read ) addBubble( make_pair( bc.bubble_+bc.len_, coord ), make_pair( bc.match_+bc.len_, read ), cm );
            
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
    if ( i == result.start_+result.len_ ) map.anchors_[1] = coord;
    if ( bc.match_ + bc.len_ == cm->node_->size() ) coord = bc.bubble_+bc.len_; // Prevent adding bogus deletions if the read has simply come to its end
    
    // Clean up
    if ( bc.len_ ) map.coords_.push_back( bc );
    if ( bc.bubble_+bc.len_ < coord || bc.match_ + bc.len_ < read ) addBubble( make_pair( bc.bubble_+bc.len_, coord ), make_pair( bc.match_+bc.len_, read ), cm );
    for ( ; i < bac.end_; i++ ) if ( result.s_[1][i] != '-' ) read++;
    
    if ( read && !finished ) maps_.push_back( map );
}

bool Bubble::cleanup()
{
    bool cleaned = maps_.empty();
    for ( int i = 0; i < bubs_.size(); i++ ) if ( bubs_[i]->cleanup() &&  bubs_[i]->maps_.empty() )
    {
        delete bubs_[i];
        bubs_.erase( bubs_.begin() + i-- );
        cleaned = true;
    }
    return cleaned;
}

bool Bubble::consolidate( vector<SNPs*>& snps, vector<Bubble*>& bubs, vector<Bubble*>& baseBubs, string& baseSeq )
{
    if ( template_.empty() ) return false;
    
    vector<Bubble::BubbleAltPath> matches;
    for ( Bubble::BubbleAltPath& bap : Bubble::BubbleAltPath::build( this ) )
    {
        SnpAlignment align( baseSeq, bap.seq_, start_, len_, 0, bap.seq_.size(), snps, bubs );
        if ( bap.set( align.align( true, true ) ) && bap.compatible( this, bap.bestAlign_  ) ) matches.push_back( bap );
    }
    if ( matches.empty() ) return false;
    
    sort( matches.begin(), matches.end(), []( Bubble::BubbleAltPath& a, Bubble::BubbleAltPath& b ){ return a.score_ > b.score_; } );
    Bubble::test( baseBubs );
    matches[0].presplit( this, matches[0].bestAlign_ );
    
    for ( Bubble* b : matches[0].split( this, matches[0].bestAlign_ ) ) baseBubs.push_back( b );
    vector<Bubble*> tests;
    for ( Bubble* b : baseBubs ) if ( b != this ) tests.push_back( b );
    Bubble::test( tests );
    return true;
    
//    for ( Bubble::BubbleAltPath& bap : Bubble::BubbleAltPath::build( this ) )
//    for ( Bubble::BubbleAltPath& bap : getAlignableSeqs() )
//    for ( pair<string,vector<Bubble*>>& as : getAlignableSeqs() )
//    {
//        SnpAlignment align( baseSeq, bap.seq_, start_, len_, 0, bap.seq_.size(), snps, bubs );
//        bap.result_ = align.align( true, true );
//        for ( int k = 0; k < bap.result_.bubbles_.size(); k++ )
////        for ( SnpAlignResult::BubbleAlignCoords& bac : result.bubbles_ )
//        {
//            SnpAlignResult::BubbleAlignCoords bac = bap.result_.bubbles_[k];
//            int scores[3] = { 0, 0, 2 };
//            for ( int i = 0; i < bap.result_.s_[0].size(); i++ )
//            {
//                int gap = 0;
//                for ( int s : { 0, 1 } ) if ( bap.result_.s_[s][i] == '-' ) while ( i < bap.result_.s_[0].size() && bap.result_.s_[s][i++] == '-' ) gap++;
//                int j = i < bac.start_ ? 0 : ( i < bac.end_ ? 2 : 1 );
//                if ( gap )
//                {
//                    if ( i > bac.start_ && i-gap < bac.end_ ) j = 2;
//                    scores[j] += 3 + gap;
//                    i--;
//                }
//                else scores[j] += bap.result_.s_[0][i] == bap.result_.s_[1][i] ? 1 : -1;
//            }
//            int score = scores[2] + max( 0, min( scores[0], scores[1] ) );
//            for ( Bubble* b : bap.split( this, bac, baseBubs ) ) baseBubs.push_back( b );
            
//            bap.assign( bac.coord_ );
//            vector<Bubble*> stacks[3]{ bap.getStack( 0 ), bap.getStack( 1 ), bap.getStack( 2 ) };
            
//            int cuts[2]{-1,-1};
////            as.second.insert( as.second.begin(), this );
////            vector<Bubble*> stacks[3];
////            vector<int> coords( as.second.size(), 0 );
////            int coord = 0;
//            bool ascending  = true;
//            for ( int i = 0; i < bap.path_.size(); i++ )
//            {
//                if ( bap.path_[i].second.first <= bac.coord_[0] && bac.coord_[0] <= bap.path_[i].second.second )
//                {
//                    if ( !i || !ascending || bap.path_[i].second.first < bac.coord_[0] ) cuts[0] = i;
//                }
//                if ( i == bap.path_.size()/2 ) ascending = false;
//                if ( bap.path_[i].second.first <= bac.coord_[1] && bac.coord_[1] <= bap.path_[i].second.second )
//                {
//                    if ( i+1 == bap.path_.size() || ascending || bac.coord_[1] < bap.path_[i].second.second ) cuts[1] = i;
//                }
//            }
//            cut( bap.getStack( 0 ), bap.getStack( 2 ), bap.coord_[0], 0 );
//            cut( bap.getStack( 1 ), bap.getStack( 2 ), bap.coord_[1], 1 );
//            assert( false );
//        }
//        assert( false );
//    }
//    assert( false );
//    return false;
}

void Bubble::cut( vector<Bubble*> cutStack, vector<Bubble*> keepStack, int cut, bool drxn )
{
    assert( false );
}

//vector<Bubble::BubbleAltPath> Bubble::getAlignableSeqs()
//{
//    vector<Bubble::BubbleAltPath> paths( 1, Bubble::BubbleAltPath() );
//    paths[0].seq_ = template_;
//    paths[0].path_.push_back( make_pair( this, make_pair( 0, template_.size() ) ) );
//    for ( Bubble* b : bubs_ ) for ( Bubble::BubbleAltPath bap : b->getAlignableSeqs() )
//    {
//        bap.seq_ = template_.substr( 0, b->start_ ) + bap.seq_ + template_.substr( b->start_ + b->len_ );
//        for ( pair<Bubble*,pair<int,int>>& path : bap.path_ ) path.second = make_pair( path.second.first+b->start_, path.second.second+b->start_);
//        bap.path_.insert( bap.path_.begin(), make_pair( this, make_pair( 0, b->start_ ) ) );
//        bap.path_.insert( bap.path_.end(), make_pair( this, make_pair( bap.path_.back().second.second, bap.seq_.size() ) ) );
//    }
//    return paths;
//}

//vector<pair<string,vector<Bubble*>>> Bubble::getAlignableSeqs()
//{
//    vector<pair<string,vector<Bubble*>>> seqs{ make_pair( template_, vector<Bubble*>{} ) };
//    for ( Bubble* b : bubs_ ) for ( pair<string,vector<Bubble*>> as : b->getAlignableSeqs() )
//    {
//        string seq = template_.substr( 0, b->start_ ) + as.first + template_.substr( b->start_ + b->len_ );
//        as.second.insert( as.second.begin(), b );
//        seqs.push_back( make_pair( seq, as.second ) );
//    }
//    return seqs;
//}

string Bubble::getConsensus()
{
    sort( bubs_.begin(), bubs_.end(), []( Bubble* a, Bubble* b ){ return a->start_ == b->start_ ? a->start_+a->len_ > b->start_+b->len_ : a->start_ < b->start_; } );
    vector<int> bubs;
    // Find dominant bubbles
    for ( int i = 0; i < bubs_.size(); )
    {
        int len = 1, best = i;
        while ( i+len < bubs_.size() && bubs_[i]->isClash( bubs_[i+len] ) ) len++;
        for ( int j = i+1; j < i+len; j++ ) if ( bubs_[j]->score_ > bubs_[best]->score_ ) best = j;
        if ( best == i && bubs_[i]->score_ > 0 ) bubs.push_back( i );
        if ( best == i ) i += len;
        else i = best;
    }
    
    // Fill in any positive bubbles that do not overlap with dominant bubbles
    for ( int i = 1; i < bubs.size(); i++ )
    {
        int best = bubs[i];
        for ( int j = bubs[i-1]+1; j < bubs[i]; j++ ) if ( !bubs_[j]->isClash( bubs_[bubs[i-1]] ) && !bubs_[j]->isClash( bubs_[bubs[i]] ) )
        {
            if ( best == bubs[i] || bubs_[j]->score_ > bubs_[best]->score_ ) best = j;
        }
        if ( best != bubs[i] && bubs_[best]->score_ > 0 ) bubs.insert( bubs.begin() + i--, best );
    }
    
    // Get consensus
    string seq;
    int coord = 0;
    for ( int b : bubs )
    {
        seq += template_.substr( coord, bubs_[b]->start_-coord ) + bubs_[b]->getConsensus();
        coord = bubs_[b]->start_ + bubs_[b]->len_;
    }
    if ( coord < template_.size() ) seq += template_.substr( coord );
    return seq;
}

vector<ConMap*> Bubble::getMapped()
{
    vector<ConMap*> mapped;
    for ( BubbleMap& bm : maps_ ) mapped.push_back( bm.map_ );
    return mapped;
}


void Bubble::getMapRange( ConMap* cm, int& start, int& end )
{
    for ( BubbleMap& bm : maps_ ) if ( bm.map_ == cm )
    {
        if ( !bm.coords_.empty() ) start = min( start, bm.coords_[0].match_ );
        if ( !bm.coords_.empty() ) end = max( end, bm.coords_.back().match_+bm.coords_.back().len_ );
        break;
    }
    for ( Bubble* b : bubs_ ) b->getMapRange( cm, start, end );
}

bool Bubble::getOffset( ConMap* l, ConMap* r, int& off )
{
    BubbleMap* maps[2]{ NULL, NULL };
    for ( int i = 0; i < maps_.size(); i++ ) for ( int d : { 0, 1 } ) if ( maps_[i].map_ == ( d ? r : l ) ) maps[d] = &maps_[i];
    if ( !maps[0] || !maps[1] ) return false;
    
    for ( BubbleCoord& lbc : maps[0]->coords_ ) for ( BubbleCoord& rbc : maps[1]->coords_ ) if ( lbc.bubble_ < rbc.bubble_+rbc.len_ && rbc.bubble_ < lbc.bubble_+lbc.len_ )
    {
        int diff = rbc.bubble_ - lbc.bubble_;
        off = lbc.match_+diff-rbc.match_;
        return true;
    }
    
    for ( Bubble* b : bubs_ ) if ( b->getOffset( l, r, off ) ) return true;
    return false;
}

bool Bubble::isAnchored()
{
    sort( maps_.begin(), maps_.end(), []( BubbleMap& a, BubbleMap& b ) { return a.anchors_[0] == b.anchors_[0] ? a.anchors_ [1] > b.anchors_[1] : a.anchors_[0] < b.anchors_[0]; });
    int anchor = 0;
    for ( BubbleMap& bm : maps_ ) if ( !bm.anchors_[0] || bm.anchors_[0] < anchor )
    {
        anchor = max( bm.anchors_[1], anchor );
        if ( anchor == template_.size() ) return true;
    }
    return false;
}

bool Bubble::isClash( Bubble* b )
{
    if ( b->start_ == start_ && b->len_ == len_ ) return true;
    return start_ < b->start_ + b->len_ && b->start_ < start_ + len_;
}

bool Bubble::isClash( BubbleCoord& bc )
{
    if ( bc.match_ == start_ && bc.len_ == len_ ) return true;
    return start_ < bc.match_ + bc.len_ && bc.match_ < start_ + len_;
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

bool Bubble::remove( ConMap* cm, bool doErase )
{
    bool removed = false;
    for ( int i = 0; i < maps_.size(); i++ ) if ( maps_[i].map_ == cm ) 
    {
        maps_.erase( maps_.begin() + i-- );
        removed = true;
        break;
    }
    for ( int i = 0; i < bubs_.size(); i++ )  if ( bubs_[i]->remove( cm ) && bubs_[i]->maps_.empty() && doErase )
    {
        delete bubs_[i];
        bubs_.erase( bubs_.begin() + i-- );
    }
    return removed;
}

void Bubble::remove( unordered_set<ConMap*>& removes )
{
    for ( int i = 0; i < maps_.size(); i++ ) if ( removes.find( maps_[i].map_ ) != removes.end() ) maps_.erase( maps_.begin() + i-- );
    for ( int i = 0; i < bubs_.size(); i++ )
    {
        bubs_[i]->remove( removes );
        if ( !bubs_[i]->maps_.empty() ) continue;
        delete bubs_[i];
        bubs_.erase( bubs_.begin() + i-- );
    }
    
}

bool Bubble::retract( vector<ConMap*>& orphans )
{
    if ( isAnchored() ) return false;
    for ( BubbleMap& bm : maps_ ) 
    {
        int coord[2]{ bm.map_->node_->size(), 0 };
        getMapRange( bm.map_, coord[0], coord[1] );
        if ( find( orphans.begin(), orphans.end(), bm.map_ ) == orphans.end() ) orphans.push_back( bm.map_ );
        if ( coord[1] <= coord[0] ) continue;
        if ( bm.map_->range_[1] < coord[1] ) bm.map_->mapped_[1] = min( coord[0]-1, bm.map_->range_[1] );
        if ( coord[0] < bm.map_->range_[0] ) bm.map_->mapped_[0] = max( coord[1], bm.map_->range_[0] );
    }
    return true;
}

int Bubble::score( unordered_set<ConMap*>& maps )
{
    for ( BubbleMap& bm : maps_ ) maps.erase( bm.map_ );
    score_ = maps_.size()-maps.size();
    for ( Bubble* b : bubs_ )
    {
        unordered_set<ConMap*> bmaps;
        for ( BubbleMap& bm : maps_ ) if ( bm.anchors_[0] <= b->start_ && b->start_+b->len_ <= bm.anchors_[1] ) bmaps.insert( bm.map_ );
        for ( Bubble* bb : bubs_ ) if ( ( b->start_ == bb->start_ && b->len_ == bb->len_ ) || ( b->start_ < bb->start_+bb->len_ && bb->start_ < b->start_+b->len_ ) ) for ( BubbleMap& bm : bb->maps_ ) bmaps.erase( bm.map_ );
        b->score( bmaps );
    }
    return score_;
}

void Bubble::setRemapped( vector<pair<int,int>>& remapped, unordered_map<Bubble*,int>& oldBubbles, vector<Bubble*>& curBubbles )
{
    for ( Bubble* b : curBubbles )
    {
        auto it = oldBubbles.find( b );
        if ( it == oldBubbles.end() || it->second < b->maps_.size() ) remapped.push_back( make_pair( b->start_, b->start_+b->len_ ) );
    }
}

vector<Bubble*> Bubble::split( vector<Bubble*> stacks[2], int coord[2], int anchor )
{
    unordered_set<ConMap*> missed, taken;
    for ( int d : { 0, 1 } ) for ( Bubble* b : stacks[d].back()->bubs_ ) if ( b->start_ <= coord[d] && coord[d] < b->start_+b->len_ ) for ( BubbleMap& bm : b->maps_ ) missed.insert( bm.map_ );
    for ( int d : { 0, 1 } ) for ( int i = stacks[d].size(); --i > 0; )
    {
        unordered_set<ConMap*> keep;
        for ( BubbleMap& bm : stacks[d][i]->maps_ ) keep.insert( bm.map_ );
        for ( Bubble* b : stacks[d][i-1]->bubs_ ) if ( b->isClash( stacks[d][i] ) ) for ( BubbleMap& bm : b->maps_ ) if ( keep.find( bm.map_ ) == keep.end() ) missed.insert( bm.map_ );
        for ( BubbleMap& bm : stacks[d][i-1]->maps_ ) for ( BubbleCoord& bc : bm.coords_ ) if ( stacks[d][i]->isClash( bc ) ) missed.insert( bm.map_ );
    }
    
    Bubble* bubs[2]{ new Bubble(), new Bubble() };
    bubs[0]->type_ = bubs[1]->type_ = type_;
    bubs[0]->start_ = start_;
    bubs[0]->len_ = anchor-start_;
    bubs[1]->start_ = anchor;
    bubs[1]->len_ = start_+len_-anchor;
    
    for ( int d : { 0, 1 } ) for ( int i = d ? stacks[1].size()-1 : 0; i >= 0 && i < stacks[d].size(); d ? i-- : i++ )
    {
        Bubble* up = i+1 < stacks[d].size() ? stacks[d][i+1] : NULL;
        int cutoff = up ? ( d ? up->start_+up->len_ : up->start_ ) : ( d ? coord[1] : coord[0]+1);
        int offset = d ? -cutoff : bubs[0]->template_.size();
        int len = d ? stacks[d][i]->template_.size()-cutoff : cutoff;
        // Take mapped reads
        for ( BubbleMap& bm : stacks[0][i]->maps_ ) if ( missed.find( bm.map_ ) == missed.end() )
        {
            bubs[d]->addCloned( bm, bubs[d]->template_.size(), bubs[d]->template_.size()+len, offset );
        }
        // Take sub-bubbles
        for ( Bubble* b : stacks[d][i]->bubs_ ) if ( b != up && ( d ? cutoff <= b->start_ : b->start_+b->len_ <= cutoff ) )
        {
            bubs[d]->addCloned( b->stealBubble( missed, offset ) );
        }
        // Add sequence
        bubs[d]->template_ += stacks[d][i]->template_.substr( d ? cutoff : 0, len );
    }
    
    for ( int d : { 0, 1 } ) for ( BubbleMap& bm : bubs[d]->maps_ ) taken.insert( bm.map_ );
    this->remove( taken );
    
//    for ( int i = 0; i < stacks[0].size(); i++ )
//    {
//        Bubble* nxt = i+1 < stacks[0].size() ? stacks[0][i+1] : NULL;
//        int cutoff = nxt ? nxt->start_ : coord[0]+1, offset = bubs[0]->template_.size();
//        for ( BubbleMap& bm : stacks[0][i]->maps_ ) if ( missed.find( bm.map_ ) == missed.end() )
//        {
//            bubs[0]->addCloned( bm, 0, offset+cutoff, offset );
//        }
//        for ( Bubble* b : stacks[0][i]->bubs_ ) if ( b != nxt && b->start_+b->len_ <= cutoff )
//        {
//            assert( false );
//            if ( clone ) bubs[0]->addCloned( b->stealBubble( missed, offset ) );
//        }
//        bubs[0]->template_ += stacks[0][i]->template_.substr( 0, cutoff );
//    }
//    for ( int i = stacks[1].size(); i-- > 0; )
//    {
//        Bubble* prv = i+1 < stacks[1].size() ? stacks[1][i+1] : NULL;
//        int cutoff = prv ? prv->start_+prv->len_ : coord[1];
//        for ( Bubble* b : stacks[1][i]->bubs_ ) if ( b != prv && cutoff <= b->start_+b->len_ )
//        {
//            Bubble* clone = b->stealBubble( missed, bubs[1]->template_.size()-cutoff );
//            assert( false );
//            if ( clone ) bubs[1]->bubs_.push_back( clone );
//        }
//        bubs[1]->template_ += stacks[1][i]->template_.substr( cutoff );
//    }
    assert( false );
}

void Bubble::addCloned( Bubble* b )
{
    if ( !b ) return;
    bubs_.push_back( b );
    for ( BubbleMap& bm : b->maps_ )
    {
        bool added = false;
        for ( BubbleMap& obm : maps_ ) if ( added = obm.map_ == bm.map_ )
        {
            bm.anchors_[0] = min( bm.anchors_[0], b->start_ );
            bm.anchors_[1] = max( bm.anchors_[1], b->start_+b->len_ );
            break;
        }
        if ( added ) continue;
        BubbleMap nbm( bm.map_ );
        nbm.anchors_[0] = b->start_;
        nbm.anchors_[1] = b->start_+b->len_;
        maps_.push_back( nbm );
    }
}

void Bubble::addCloned( BubbleMap bm, int start, int end, int offset )
{
    for ( int i = 0; i < bm.coords_.size(); i++ )
    {
        bm.coords_[i].bubble_ -= offset;
        int trim[2]{ max( start - bm.coords_[i].bubble_, 0 ), max( bm.coords_[i].bubble_ + bm.coords_[i].len_ - end, 0 ) };
        bm.coords_[i].bubble_ += trim[0];
        bm.coords_[i].match_ += trim[0];
        bm.coords_[i].len_ -= trim[0] + trim[1];
        if ( bm.coords_[i].len_ <= 0 ) bm.coords_.erase( bm.coords_.begin() + i-- );
    }
    if ( bm.coords_.empty() ) return;
    bm.anchors_[0] = bm.coords_[0].bubble_;
    bm.anchors_[1] = bm.coords_.back().bubble_ + bm.coords_.back().len_;
    for ( BubbleMap& obm : maps_ ) if ( obm.map_ == bm.map_ )
    {
        for ( BubbleCoord& bc : bm.coords_ ) obm.addCoord( bc.match_, bc.bubble_, bc.len_ );
        return;
    }
    maps_.push_back( bm );
}

Bubble* Bubble::stealBubble( unordered_set<ConMap*>& ignores, int startOffset )
{
    Bubble* b = new Bubble();
    b->template_ = template_;
    b->type_ = type_;
    b->start_ = start_ + startOffset;
    b->len_ = len_;
    for ( int i = 0; i < b->maps_.size(); i++ ) if ( ignores.find( b->maps_[i].map_ ) == ignores.end() )
    {
        b->maps_.push_back( b->maps_[i] );
        b->maps_.erase( b->maps_.begin() + i-- );
    }
    if ( b->maps_.empty() )
    {
        assert( false );
        delete b;
        return NULL;
    }
    for ( Bubble* sb : bubs_ )
    {
        Bubble* clone = sb->stealBubble( ignores, 0 );
        if ( clone ) b->bubs_.push_back( sb );
    }
    return b;
}

void Bubble::test()
{
    for ( BubbleMap& bm : maps_ ) for ( BubbleCoord& bc : bm.coords_ ) assert( bc.len_ );
    for ( Bubble* b : bubs_ ) b->test();
}

void Bubble::test( vector<Bubble*>& bubs )
{
    for ( Bubble* b : bubs ) 
    {
        b->test();
        assert( b->isAnchored() );
        int anchors[2]{ 0, (int)b->template_.size() };
        for ( BubbleMap& bm : b->maps_ )
        {
            assert( bm.anchors_[0] >= 0 && bm.anchors_[1] <= b->template_.size() );
            if ( bm.anchors_[0] == 0 ) anchors[0] = max( anchors[0], bm.anchors_[1] );
            if ( bm.anchors_[1] == b->template_.size() ) anchors[1] = min( anchors[1], bm.anchors_[0] );
        }
        assert( b->template_.empty() || anchors[1] < anchors[0] );
    }
}
