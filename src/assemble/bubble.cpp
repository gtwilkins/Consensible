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
    assert( bac.snp_->seq_.empty() || result.s_[0][i] == result.s_[1][i] || !len_ );
    int len = bac.snp_->seq_.size();
    if ( !bac.snp_->seq_.empty() && result.s_[0][i] != result.s_[1][i] )
    {
        string seq = result.s_[1].substr( i, bac.end_-bac.start_ );
        len = seq.size();
        assert( seq.size() == 1 && seq != "-" );
        addSnp( seq, cm->node_, read );
    }
    else
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

//void Bubble::BubbleMap::getSeq()
//{
//    
//}
//
//void Bubble::BubbleMap::getType()
//{
//    
//}

Bubble::BubbleMapDetails::BubbleMapDetails( Bubble* bub, BubbleMap& map )
: map_( map.map_ )
{
    assert( map.coords_.size() );
    BubbleCoord* coords[2] = { &map.coords_[0], &map.coords_.back() };
    start_ = coords[0]->match_;
    end_ = coords[1]->match_ + coords[1]->len_;
    off_ = 0;
    bool complete[2]{ !map.anchors_[0], map.anchors_[1] == bub->template_.size() };
    if ( bub->type_ == 2 ) assert( complete[0] && complete[1] );
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
        size_t l = seq_.find_first_of( map.seq_ ), r = seq_.find_last_of( map.seq_ );
        if ( map.type_ == 1 && map_[0].type_ != 0 && l != string::npos && l == 0 ) added = true;
        if ( map.type_ == 0 && map_[0].type_ != 1 && r != string::npos && r == seq_.size()-map.seq_.size() ) added = true;
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

vector<Bubble::BubbleAltPath> Bubble::BubbleAltPath::build( Bubble* b )
{
    vector<BubbleMapDetails> maps;
    for ( BubbleMap& bm : b->maps_ ) maps.push_back( BubbleMapDetails( b, bm ) );
    sort( maps.begin(), maps.end(), []( BubbleMapDetails& a, BubbleMapDetails& b ){ return a.seq_.size() > b.seq_.size(); } );
    vector<BubbleAltPath> paths;
    for ( BubbleMapDetails& bmd : maps )
    {
        bool added = false;
        for ( BubbleAltPath& path : paths ) if ( path.add( bmd ) ) added = true;
        if ( !added ) paths.push_back( BubbleAltPath( bmd ) );
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
        for ( BubbleMapDetails& map : map_ ) if ( map.off_ <= coord && coord+seq.size() <= map.seq_.size() ) bac.snps_->addSnp( seq, map.map_->node_, map.start_+coord-map.off_ );
        assert( bac.bubbles_.empty() );
        coord += seq.size();
        return;
    }
    vector<pair<BubbleMap, BubbleMapDetails*>> maps;
    for ( BubbleMapDetails& map : map_ )
    {
        BubbleMap bm( map.map_ );
        bm.anchors_[0] = bac.bubble_->template_.size();
        bm.anchors_[1] = 0;
        maps.push_back( make_pair( bm, &map ) );
//        int lExcess = coord - map.off_ + map.start_;
//        bm.anchors_[0] = max( 0, map.off_-coord );
//        bm.anchors_[1];
    }
    int b = 0, j = 0, nxt = bac.bubbles_.empty() ? bac.end_ : bac.bubbles_[0].start_, coordStart = coord;
    for ( int i = bac.start_; i < bac.end_; )
    {
        if ( i == nxt )
        {
            assert( false );
            for ( pair<BubbleMap, BubbleMapDetails*>& map : maps ) if ( coord < map.second->seq_.size() ) map.first.anchors_[1] = max( map.first.anchors_[1], bac.bubbles_[j].bubble_->start_+bac.bubbles_[j].bubble_->len_ );
            remap( bac.bubbles_[j], coord );
            for ( pair<BubbleMap, BubbleMapDetails*>& map : maps ) if ( map.second->off_ < coord ) map.first.anchors_[0] = min( map.first.anchors_[0], bac.bubbles_[j].bubble_->start_+bac.bubbles_[j].bubble_->len_ );
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
            assert( false );
            Bubble* bub = Bubble::addBubble( seq, bac.bubble_->bubs_, start, len );
            for ( pair<BubbleMap, BubbleMapDetails*>& map : maps ) if ( map.second->off_ < coord+seq.size() && coord-map.second->off_ < map.second->seq_.size() )
            {
                BubbleMap bm( map.second->map_ );
                bm.anchors_[0] = max( 0, map.second->off_-coord );
                bm.anchors_[1] = min( seq.size(), map.second->seq_.size()-coord-map.second->off_ );
                map.first.anchors_[0] = min( bub->start_, map.first.anchors_[0] );
                map.first.anchors_[1] = max( bub->start_ + bub->len_, map.first.anchors_[1] );
                assert( bm.anchors_[1] < bm.anchors_[0] && bm.anchors_[1]-bm.anchors_[0] <= seq.size() );
                BubbleCoord bc( bm.anchors_[0], bm.anchors_[0] + max( 0, coord-map.second->off_ ), bm.anchors_[1]-bm.anchors_[0] );
                assert( bub->template_.substr( bc.bubble_, bc.len_ ) == map.second->map_->node_->read_->seq_.substr( bc.match_, bc.len_ ) );
                bm.coords_.push_back( bc );
                bub->maps_.push_back( bm );
            }
//            for ( BubbleMapDetails& map : map_ ) if ( map.off_ < coord+seq.size() && coord < map.seq_.size() )
//            {
//                BubbleMap bm( map.map_ );
//                bm.anchors_[0] = max( 0, map.off_-coord );
//                bm.anchors_[1] = min( seq.size(), map.seq_.size()-coord-map.off_ );
//                assert( bm.anchors_[1] < bm.anchors_[0] && bm.anchors_[1]-bm.anchors_[0] <= seq.size() );
//                BubbleCoord bc( bm.anchors_[0], bm.anchors_[0] + max( 0, coord-map.off_ ), bm.anchors_[1]-bm.anchors_[0] );
//                assert( bub->template_.substr( bc.bubble_, bc.len_ ) == map.map_->node_->read_->seq_.substr( bc.match_, bc.len_ ) );
//                bm.coords_.push_back( bc );
//                bub->maps_.push_back( bm );
//            }
            coord += seq.size();
        }
        else
        {
            assert( result_.s_[0][i] != '-' );
            for ( pair<BubbleMap, BubbleMapDetails*>& map : maps ) if ( map.second->off_ <= coord && coord-map.second->off_ < map.second->seq_.size() )
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
    for ( pair<BubbleMap, BubbleMapDetails*>& map : maps ) if ( map.second->off_ < coord && coordStart < map.second->seq_.size() )
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
    for ( SnpAlignResult::BubbleAlignCoords& bac : result_.bubbles_ )
    {
        int scores[3] = { 0, 0, 2 };
        for ( int i = 0; i < result_.s_[0].size(); i++ )
        {
            int gap = 0;
            for ( int s : { 0, 1 } ) if ( result_.s_[s][i] == '-' ) while ( i < result_.s_[0].size() && result_.s_[s][i++] == '-' ) gap++;
            int j = i < bac.start_ ? 0 : ( i < bac.end_ ? 2 : 1 );
            if ( gap )
            {
                if ( i > bac.start_ && i-gap < bac.end_ ) j = 2;
                scores[j] += 3 + gap;
                i--;
            }
            else scores[j] += result_.s_[0][i] == result_.s_[1][i] ? 1 : -1;
        }
        int score = scores[2] + max( 0, min( scores[0], scores[1] ) );
        if ( score <= score_ ) continue;
        score_ = score;
        bestAlign_ = bac;
    }
    return score_ >= 0;
}

vector<Bubble*> Bubble::BubbleAltPath::split( Bubble* base, SnpAlignResult::BubbleAlignCoords& bac )
{
    vector<Bubble*> created;
    if ( bac.coord_[0] )
    {
        created.push_back( new Bubble( base, *this, bac, bac.coord_[0], 0 ) );
    }
    if ( bac.coord_[1] < seq_.size() ) created.push_back( new Bubble( base, *this, bac, bac.coord_[1], 1 ) );
    
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
    for ( BubbleMapDetails& map : path.map_ ) if ( map.off_ < coord[1] && coord[0] < map.seq_.size() )
    {
//        int off = !map.type_ ? path.seq_.size()-map.seq_.size() : 0;
//        if ( coord[1] <= cut || map.seq_.size()+off <= coord[0] ) assert( false );
//        if ( coord[1] <= cut || map.seq_.size()+off <= coord[0] ) continue;
        BubbleMap bm( map.map_ );
        bm.anchors_[0] = template_.size();
        bm.anchors_[1] = 0;
        for ( int i = max( map.off_, coord[0] ); i < min( map.off_+(int)map.seq_.size(), coord[1] ); i++ ) assert( map.map_->node_->read_->seq_[map.start_+i-map.off_] == template_[i-coord[0]] );
        for ( int i = max( map.off_, coord[0] ); i < min( map.off_+(int)map.seq_.size(), coord[1] ); i++ ) bm.addCoord( map.start_+i-map.off_, i-coord[0] );
        maps_.push_back( bm );
    } else assert( false );
}

Bubble::Bubble( ConMap* cm, bool drxn )
: start_( drxn ? cm->coord_[1]+1 : cm->coord_[0] ), len_( 0 ), type_( drxn ), score_( 0 )
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
    map.coords_.push_back( bc );
    bub->maps_.push_back( map );
}

void Bubble::addBubble( pair<int,int> bubble, pair<int,int> read, ConMap* cm )
{
    int len = bubble.second-bubble.first;
    string seq = cm->node_->read_->seq_.substr( read.first, read.second-read.first );
    
    Bubble* bub = addBubble( cm->node_->read_->seq_.substr( read.first, read.second-read.first ), bubs_, bubble.first, len );
    BubbleMap map( cm );
    BubbleCoord bc( 0, read.first, bub->template_.size() );
    map.anchors_[0] = 0;
    map.anchors_[1] = bub->template_.size();
    map.coords_.push_back( bc );
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
    BubbleCoord bc( 0, start, len_ );
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

bool Bubble::consolidate( vector<SNPs*>& snps, vector<Bubble*>& bubs, vector<Bubble*>& baseBubs, string& baseSeq )
{
    vector<Bubble::BubbleAltPath> matches;
    for ( Bubble::BubbleAltPath& bap : Bubble::BubbleAltPath::build( this ) )
    {
        SnpAlignment align( baseSeq, bap.seq_, start_, len_, 0, bap.seq_.size(), snps, bubs );
        if ( bap.set( align.align( true, true ) ) ) matches.push_back( bap );
    }
    sort( matches.begin(), matches.end(), []( Bubble::BubbleAltPath& a, Bubble::BubbleAltPath& b ){ return a.score_ > b.score_; } );
    if ( !matches.empty() )
    {
        for ( Bubble* b : matches[0].split( this, matches[0].bestAlign_ ) ) baseBubs.push_back( b );
        return true;
    }
    return false;
    
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

bool Bubble::isClash( Bubble* b )
{
    if ( b->start_ == start_ && b->len_ == len_ ) return true;
    return start_ < b->start_ + b->len_ && b->start_ < start_ + len_;
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

bool Bubble::remove( ConMap* cm )
{
    bool removed = false;
    for ( int i = 0; i < maps_.size(); i++ ) if ( maps_[i].map_ == cm ) 
    {
        maps_.erase( maps_.begin() + i-- );
        removed = true;
        break;
    }
    for ( int i = 0; i < bubs_.size(); i++ )  if ( bubs_[i]->remove( cm ) && bubs_[i]->maps_.empty() )
    {
        delete bubs_[i];
        bubs_.erase( bubs_.begin() + i-- );
    }
    return removed;
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
