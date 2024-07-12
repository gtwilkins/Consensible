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

#include "long_query.h"
#include "parameters.h"
#include "constants.h"
#include "alignment.h"
#include <algorithm>
#include <cassert>
#include <iostream>

extern Parameters params;


size_t ReadMatch::Align::gapLen( Align& rhs )
{
    assert( read_+len_ <= rhs.read_ && tar_+len_ <= rhs.tar_ );
    return max( rhs.read_ - ( read_+len_ ), rhs.tar_ - ( tar_+len_ ) );
}

size_t ReadMatch::Align::readEnd()
{
    return read_+len_;
}

size_t ReadMatch::Align::tarEnd()
{
    return tar_+len_;
}

size_t ReadMatch::Align::totalLen( Align& rhs )
{
    assert( read_+len_ <= rhs.read_ && tar_+len_ <= rhs.tar_ );
    return min( rhs.read_+len_ - read_, rhs.tar_+len_ - tar_ );
}

ReadMatch::ReadMatch( string seq, string& base, ReadId id, size_t coord, size_t len, bool drxn )
: seq_( seq ), id_( id )
{
    assert( len && len <= seq.size() );
    Align al;
    al.tar_ = drxn ? coord-len : coord+len-1;
    al.read_ = drxn ? seq.size()-len : len-1;
    al.len_ = 1;
    
    while ( al.read_ && al.tar_ && seq[al.read_-1] == base[al.tar_-1] && al.read_-- && al.tar_-- ) al.len_++;
    while ( al.read_+al.len_ < seq.size() && al.tar_+al.len_ < base.size() && seq[al.read_+al.len_] == base[al.tar_+al.len_] ) al.len_++;
    
    size_t readLens[2]{ al.read_, seq.size()-( al.read_+al.len_ ) };
    size_t readStarts[2]{ 0, al.read_+al.len_ };
    size_t baseLens[2]{ min( al.tar_, readLens[0]+15 ), min( base.size()-( al.tar_ +al.len_ ), readLens[1]+15 ) };
    size_t baseStarts[2]{ al.tar_-baseLens[0], al.tar_+al.len_ };
    
    for ( int d : { 0, 1 } )
    {
        Alignment align( base, seq, baseStarts[d], baseLens[d], readStarts[d], readLens[d] );
        AlignResult ar = align.align( 1, 2, 4, 4, d, !d );
        size_t coords[2]{ baseStarts[d], readStarts[d] };
        Align cur;
        for ( int i = 0; i < ar.s_[0].size(); i++ )
        {
            if ( ar.s_[0][i] == ar.s_[1][i] && i >= ar.start_ && i < ar.start_+ar.len_ )
            {
                assert( ar.s_[0][i] != '-' );
                if ( cur.len_ ) cur.len_++;
                else cur = Align( coords[0], coords[1], 1 );
            }
            else if ( cur.len_ )
            {
                align_.push_back( cur );
                cur = Align();
            }
            for ( int s : { 0, 1 } ) if ( ar.s_[s][i] != '-' ) coords[s]++;
        }
        if ( cur.len_ ) align_.push_back( cur );
        if ( !d ) align_.push_back( al );
    }
    
    read_[0] = align_[0].read_;
    read_[1] = align_.back().read_+align_.back().len_;
    tar_[0] = align_[0].tar_;
    tar_[1] = align_.back().tar_+align_.back().len_;
}

bool ReadMatch::assess( uint hitLen, uint errorRate )
{
    if ( align_.size() == 1 ) return align_[0].len_ >= min( (uint)seq_.size(), hitLen );
    
    vector<float> gaps;
    for ( int i = 1; i < align_.size(); i++ ) gaps.push_back( align_[i-1].gapLen( align_[i] ) );
    
    float cutoff = (float)errorRate / (float)100;
    int bestGap = 0, bestLen = 0, best[2]{ 0, 0 };
    for ( int i = 0; i < align_.size(); i++ )
    {
        int len = align_[i].len_, gap = 0;
        if ( len >= bestLen && ( bestLen = len ) ) best[0] = best[1] = i;
        for ( int j = i; j < gaps.size(); j++ )
        {
            gap += gaps[j];
            len += align_[j+1].len_;
            if ( ( (float)gap / (float)len ) > cutoff || bestLen > len || ( bestLen == len && bestGap < gap ) ) continue;
            best[0] = i;
            best[1] = j+1;
            bestGap = gap;
            bestLen = len;
        }
    }
    
    read_[0] = align_[best[0]].read_;
    read_[1] = align_[best[1]].read_+align_[best[1]].len_;
    tar_[0] = align_[best[0]].tar_;
    tar_[1] = align_[best[1]].tar_+align_[best[1]].len_;
    
    return bestLen >= hitLen;
}

int ReadMatch::score()
{
    int score = 0;
    bool scoring = false;
    for ( int i = 0; i < align_.size(); i++ ) if ( read_[0] <= align_[i].read_ && align_[i].read_+align_[i].len_ <= read_[1] )
    {
        score += align_[i].len_;
        if ( scoring ) score -= align_[i-1].gapLen( align_[i] );
        scoring = true;
    }
    return score;
}

char ReadMatch::seq( int i )
{
    assert( i >= 0 && i < seq_.size() );
    return seq_[i];
}

string& ReadMatch::seq()
{
    return seq_;
}

size_t ReadMatch::size()
{
    return seq_.size();
}

LongQuery::LongQuery( string seq, IndexReader* ir, uint hitLen, uint errorRate )
: seq_( seq ), ir_( ir ), hitLen_( hitLen ), maxLen_( params.readLen ), errorRate_( errorRate ), tarLen_( hitLen_ ), pBlock_( 0 ), p_( 0 ), flushLen_( 1000 )
{
    tarErrors_ = ( hitLen * errorRate / 100 ) + ( ( hitLen * errorRate ) % 100 ? 1 : 0 );
    if ( errorRate ) tarLen_ = ( tarErrors_ * 100 / errorRate ) + ( ( tarErrors_ * 100 ) % errorRate ? 1 : 0 );
    
    for ( uint tarLeft = tarLen_; tarLeft > 0; )
    {
        uint block = tarLeft / ( max( (uint)1, tarErrors_ ) + 2 - blocks_.size() );
        blocks_.push_back( block );
        assert( block && block <= tarLeft );
        tarLeft -= block;
    }
    
    coords_[0] = seq_.size();
    coords_[1] = 0;
}

void LongQuery::add( LQ& lq, LQState& lqe, CharId rank, CharId count )
{
    EndRank er;
    er.coords_[!lq.drxn_] = p_ + ( lq.drxn_ ? 0 : blocks_[pBlock_] );
    er.coords_[lq.drxn_] = er.coords_[!lq.drxn_] + ( lq.drxn_ ? lqe.len_ : -lqe.len_ );
    er.errors_ = lqe.errCount_;
    for ( CharId i = 0; i < count; i++ )
    {
        auto ins = endRanks_[lq.drxn_].insert( make_pair( rank+i, er ) );
        if ( !ins.second && ins.first->second.score() < er.score() ) ins.first->second = er;
    }
    
    coords_[0] = min( coords_[0], er.coords_[0] );
    coords_[1] = max( coords_[1], er.coords_[1] );
}

bool LongQuery::advance( LQ& lq, LQState& lqe, uint8_t c )
{
    if ( lqe.len_ < lq.t_.size() && lq.t_[lqe.len_] != c && ++lqe.errCount_ > lqe.errAllowed_ && !lqe.retrieve_ )
    {
        return false;
    }
    
    if ( lqe.pBlock_ < lq.blocks_.size() && lqe.len_+1 == lq.blocks_[lqe.pBlock_] )
    {
        if ( lqe.pBlock_ >= ( lq.blocks_.size()/2 ) && lqe.errCount_ <= lqe.pBlock_ ) lqe.retrieve_ = true;
        lqe.pBlock_++;
    }
    
    if ( lqe.pErr_ < lq.errBlocks_.size() && lqe.len_+1 == lq.errBlocks_[lqe.pErr_] )
    {
        lqe.errAllowed_++;
        lqe.pErr_++;
    }
    
    if ( lqe.len_ > lq.q_.size() ) assert( false );
    if ( lqe.len_ < lq.q_.size() ) lq.q_[lqe.len_] = c;
    else lq.q_.push_back( c );
    
    return true;
}

bool LongQuery::flush( QueryBinaries* qb, bool force )
{
    assert( result_.empty() );
    if ( !result_.empty() ) return true;
    if ( coords_[1] <= coords_[0] ) return false;
    if ( !force && coords_[1] < coords_[0]+flushLen_ && p_ < coords_[1]+maxLen_ ) return false;
    
    vector<ReadMatch*> matches = matches_;
    unordered_map<ReadId, ReadMatch*> matchMap;
    unordered_set<ReadMatch*> trashed;
    for ( ReadMatch* rm : matches ) matchMap[rm->id_] = rm;
    matches_.clear();
    
    for ( int d : { 0, 1 } )
    {
        vector<pair<CharId, EndRank> > ranks( endRanks_[d].begin(), endRanks_[d].end() );
        sort( ranks.begin(), ranks.end(), []( pair<CharId, EndRank>& a, pair<CharId, EndRank>& b ){ return a.first < b.first; } );
        for ( CharId i = 0; i < ranks.size(); )
        {
            CharId count = 1, j = i;
            while ( i+count < ranks.size() && ranks[i+count].first == ranks[i].first+count ) count++;
            for ( ReadId id : qb->getIds( ranks[i].first, count ) )
            {
                if ( d ) id = params.getRevId( id );
                ReadMatch* match = new ReadMatch( qb->getSequence( id ), seq_, id, ranks[j].second.coords_[d], ranks[j].second.coords_[1]-ranks[j].second.coords_[0], d );
                bool added = false;
                if ( match->assess( hitLen_, errorRate_ ) )
                {
                    auto ins = matchMap.insert( make_pair( match->id_, match ) );
                    if ( !( added = ins.second ) && ( added = ins.first->second->score() < match->score() ) )
                    {
                        trashed.insert( ins.first->second );
                        ins.first->second = match;
                    }
                }
                if ( added ) matches.push_back( match );
                else delete match;
                j++;
            }
            i += count;
        }
        endRanks_[d].clear();
    }
    
    for ( ReadMatch* rm : matches )
    {
        if ( trashed.find( rm ) == trashed.end() ) matches_.push_back( rm );
        else delete rm;
    }
    sort( matches_.begin(), matches_.end(), []( ReadMatch* a, ReadMatch* b ){ return a->tar_[0] == b->tar_[0] ? a->tar_[1] > b->tar_[1] : a->tar_[0] < b->tar_[0]; } );
    
    bool result = package( force );
    coords_[0] = matches_.empty() ? seq_.size() : matches_[0]->tar_[0];
    coords_[1] = 0;
    for ( ReadMatch* rm : matches_ ) coords_[1] = max( rm->tar_[1], coords_[1] );
    return result;
}

LongQuery::LQ LongQuery::initiate( CharId& rank, CharId& count, bool drxn )
{
    LQ lq( drxn );
    lq.coords_[0] = lq.coords_[1] = drxn ? p_ : p_ + blocks_[pBlock_];
    if ( drxn ) for ( ; lq.coords_[1] < seq_.size() && lq.t_.size() < maxLen_; ) lq.t_.push_back( 3 - charToInt[ seq_[ lq.coords_[1]++ ] ] );
    else for ( ; lq.coords_[0] > 0 && lq.t_.size() < maxLen_; ) lq.t_.push_back( charToInt[ seq_[ --lq.coords_[0] ] ] );
    lq.q_ = lq.t_;
    ir_->setBaseAll( lq.t_[0], lq.t_[1], rank, count );
    
    for ( uint p = pBlock_; lq.blocks_.size() <= tarErrors_; )
    {
        uint blockHalves[2]{ blocks_[p] / 2, blocks_[p] - ( blocks_[p] / 2 ) };
        if ( !lq.blocks_.empty() && lq.errBlocks_.size() < tarErrors_ ) lq.errBlocks_.push_back( lq.blocks_.back() );
        if ( !lq.blocks_.empty() && lq.errBlocks_.size() < tarErrors_ ) lq.errBlocks_.push_back( lq.blocks_.back() + blockHalves[drxn] );
        lq.blocks_.push_back( blocks_[p] + ( lq.blocks_.empty() ? 0 : lq.blocks_.back() ) );
        p = ( drxn ? ( p+1 < blocks_.size() ? p+1 : 0 ) : ( p ? p-1 : blocks_.size()-1 ) );
    }
    
    if ( lq.blocks_[tarErrors_ / 2] > lq.t_.size() ) count = 0;
    
    return lq;
}

bool LongQuery::package( bool force )
{
    if ( matches_.empty() ) return false;
    
    struct ResultStruct
    {
        ResultStruct( int i, int start, int end ): i( i ), j( i ), start( start ), end( end ){};
        int i, j, start, end, gap;
    };
    
    vector<ResultStruct> blocks{ ResultStruct( 0, matches_[0]->tar_[0], matches_[0]->tar_[1] ) };
    for ( int i = 1; i < matches_.size(); i++ )
    {
        if ( matches_[i]->tar_[0] < blocks.back().end )
        {
            blocks.back().j = i;
            blocks.back().end = max( blocks.back().end, (int)matches_[i]->tar_[1] );
        }
        else blocks.push_back( ResultStruct( i, matches_[i]->tar_[0], matches_[i]->tar_[1] ) );
    }
    
    assert( result_.empty() );
    for ( int i = 0; i < blocks.size(); i++ ) blocks[i].gap = i+1 < blocks.size() ? blocks[i+1].start - blocks[i].end : p_ - blocks[i].end;
    for ( int i = 0; i < blocks.size(); i++ ) if ( blocks[i].end+maxLen_ < p_ || force )
    {
        result_.insert( result_.end(), matches_.begin(), matches_.begin() + blocks[i].j+1 );
        matches_.erase( matches_.begin(), matches_.begin() + blocks[i].j+1 );
        flushLen_ = 1000;
        return true;
    }
    
    flushLen_ += 1000;
    
    return false;
}

bool LongQuery::query()
{
    if ( p_ + blocks_[pBlock_] > seq_.size() ) return false;
    
    for ( int d : { 0, 1 } )
    {
        CharId rank, count;
        LQ lq = initiate( rank, count, d );
        if ( count ) query( lq, LQState(), rank, count, lq.t_[1] );
    }
    
    p_ += blocks_[pBlock_];
    pBlock_ = pBlock_+1 < blocks_.size() ? pBlock_+1 : 0;
    
    return true;
}

void LongQuery::query( LQ& lq, LQState lqe, CharId rank, CharId count, uint8_t c )
{
    if ( !advance( lq, lqe, c ) ) return;
    
    CharCount ranks, counts;
    ir_->countRange( lq.q_[lqe.len_++], rank, count, ranks, counts );
    
    if ( counts.endCounts && lqe.retrieve_ ) add( lq, lqe, ranks.endCounts, counts.endCounts );
    
    for ( uint8_t i = 0; i < 4; i++ ) if ( counts[i] ) query( lq, lqe, ranks[i], counts[i], i );
}

LongResult* LongQuery::yield( QueryBinaries* qb )
{
    result_.clear();
    while ( !flush( qb ) && p_ + blocks_[pBlock_] <= seq_.size() ) query();
    if ( result_.empty() && ( !matches_.empty() || !endRanks_[0].empty() || !endRanks_[1].empty() ) ) flush( qb, true );
    
    LongResult* result = result_.empty() ? NULL : new LongResult( result_ );
    return result;
}