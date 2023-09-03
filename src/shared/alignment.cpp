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

#include "alignment.h"
#include <iostream>
#include <algorithm>
#include <fstream>
#include <cassert>

Alignment::Alignment( string& a, string& b, int aStart, int aLen, int bStart, int bLen )
: a_( a.substr( aStart, aLen ) ), b_( b.substr( bStart, bLen ) ), m_( aLen+1, vector<int>( bLen+1 ) ), p_( aLen+1, vector<int>( bLen+1 ) ), scored_( false )
{
    freeEnds_[0] = freeEnds_[1] = false;
    gapExt_ = 1;
    gapOpen_ = 2;
    hit_ = 1;
    miss_ = 1;
    iMax_ = jMax_ = 0;
}

void Alignment::debug()
{
    ofstream ofs("/home/glen/Linda/glin_align_debug.csv");
    string line = ",-";
    for ( char c : a_ ) line += ","+string( 1, c );
    ofs << line << "\n";
    for ( int j = 0; j <= b_.size(); j++ )
    {
        line = j ? string( 1, b_[j-1] ) : "-";
        for ( int i = 0; i <= a_.size(); i++ ) line += "," + to_string( p_[i][j] );
        ofs << line << "\n";
    }
    line = "";
    ofs << line << "\n";
    line = ",-";
    for ( char c : a_ ) line += ","+string( 1, c );
    ofs << line << "\n";
    for ( int j = 0; j <= b_.size(); j++ )
    {
        line = j ? string( 1, b_[j-1] ) : "-";
        for ( int i = 0; i <= a_.size(); i++ ) line += "," + to_string( m_[i][j] );
        ofs << line << "\n";
    }
    ofs.close();
    assert( false );
}

bool Alignment::isFreeStart( int i, int j )
{
    if ( p_[i][j] != 0 || m_[i][j] != 0 ) return false;
    if ( i && j && m_[i][j] == miss_ ) return false;
    return true;
}

void Alignment::score()
{
    for ( int i = 0; i < a_.size()+1; i++ ) p_[i][0] = -i;
    for ( int j = 0; j < b_.size()+1; j++ ) p_[0][j] = j;
    for ( int i = 0; i < a_.size()+1; i++ ) m_[i][0] = freeEnds_[0] || !i ? 0 : -gapOpen_ - ( i * gapExt_ );
    for ( int j = 0; j < b_.size()+1; j++ ) m_[0][j] = freeEnds_[0] || !j ? 0 : -gapOpen_ - ( j * gapExt_ );
    vector<int> iMax( b_.size(), 0 ), jMax( a_.size(), 0 );
    for ( int i = 0; i < a_.size(); i++ )
    {
        for ( int j = 0; j < b_.size(); j++ )
        {
            m_[i+1][j+1] = m_[i][j] + ( a_[i] == b_[j] ? hit_ : -miss_ );
            p_[i+1][j+1] = 0;
            
            // Find the previous i value for which the current j was highest
            int aGapLen = j-jMax[i];
            int bGapLen = i-iMax[j];
            int aGapScore = m_[i+1][jMax[i]+1] - gapOpen_ - ( aGapLen * gapExt_ );
            int bGapScore = m_[iMax[j]+1][j+1] - gapOpen_ - ( bGapLen * gapExt_ );
            if ( aGapScore > m_[i+1][j+1] )
            {
                m_[i+1][j+1] = aGapScore;
                p_[i+1][j+1] = aGapLen;
            }
            if ( bGapScore > m_[i+1][j+1] )
            {
                m_[i+1][j+1] = bGapScore;
                p_[i+1][j+1] = -bGapLen;
            }
            if ( freeEnds_[0] && m_[i+1][j+1] < 0 )
            {
                m_[i+1][j+1] = 0;
                p_[i+1][j+1] = 0;
            }
            
            if ( p_[i+1][j+1] == 0 )
            {
                if ( m_[i+1][j+1] >= m_[i+1][jMax[i]+1] - ( aGapLen * gapExt_ ) ) jMax[i] = j;
                if ( m_[i+1][j+1] >= m_[iMax[j]+1][j+1] - ( bGapLen * gapExt_ ) ) iMax[j] = i;
                if ( m_[i+1][j+1] > m_[iMax_][jMax_])
                {
                    iMax_ = i+1;
                    jMax_ = j+1;
                }
            }
        }
    }
    scored_ = true;
}

AlignResult Alignment::align( int hit, int miss, int gapOpen, int gapExt, bool lAnchored, bool rAnchored )
{
    hit_ = hit;
    miss_ = miss;
    gapOpen_ = gapOpen;
    gapExt_ = gapExt;
    return align( lAnchored, rAnchored );
}

AlignResult Alignment::align( bool lAnchored, bool rAnchored )
{
    if ( lAnchored == freeEnds_[0] || rAnchored == freeEnds_[1] ) scored_ = false;
    freeEnds_[0] = !lAnchored;
    freeEnds_[1] = !rAnchored;
    if ( !scored_ ) score();
//    debug();
    
    AlignResult result;
    result.lIgnore[0] = result.lIgnore[1] = result.rIgnore[0] = result.rIgnore[1] = 0;
    
    string s[2]{ result.s_[0], result.s_[1] };
    int i = a_.size(), j = b_.size();
    if ( freeEnds_[1] )
    {
        result.rIgnore[0] = a_.size()-iMax_;
        result.rIgnore[1] = b_.size()-jMax_;
        while ( i-iMax_ > j-jMax_ ) s[0] += a_[--i];
        while ( j-jMax_ > i-iMax_ ) s[1] += b_[--j];
        if ( s[0].size() > s[1].size() ) s[1] += string( s[0].size(), '-');
        if ( s[1].size() > s[0].size() ) s[0] += string( s[1].size(), '-');
        while ( i > iMax_ ) s[0] += a_[--i];
        while ( j > jMax_ ) s[1] += b_[--j];
    }
    
    int start = 0, len = s[0].size();
    while ( i > 0 || j > 0 )
    {
        assert( i >= 0 && j >= 0 );
        if ( !i || !j || ( freeEnds_[0] && p_[i][j] == 0 && m_[i][j] == 0 && !( i && j && m_[i-1][j-1] == miss_ ) ) )
        {
            if ( freeEnds_[0] ) result.lIgnore[0] = i;
            if ( freeEnds_[0] ) result.lIgnore[1] = j;
            if ( freeEnds_[0] ) start = s[0].size();
            while ( i > 0 ) s[0] += a_[--i];
            while ( j > 0 ) s[1] += b_[--j];
            if ( s[0].size() > s[1].size() ) s[1] += string( s[0].size()-s[1].size(), '-');
            if ( s[1].size() > s[0].size() ) s[0] += string( s[1].size()-s[0].size(), '-');
            if ( freeEnds_[0] ) start = s[0].size() - start;
        }
        else if ( p_[i][j] == 0 )
        {
            s[0] += a_[--i];
            s[1] += b_[--j];
        }
        else if ( p_[i][j] < 0 )
        {
            for ( int k = p_[i][j]; k++ < 0; )
            {
                s[0] += a_[--i];
                s[1] += '-';
            }
        }
        else
        {
            for ( int k = p_[i][j]; k-- > 0; )
            {
                s[0] += '-';
                s[1] += b_[--j];
            }
        }
    }
//    len = s[0].size() - len - start;
//    reverse( s[0].begin(), s[0].end() );
//    reverse( s[1].begin(), s[1].end() );
//    cout << s[0] << endl << s[1] << endl ;
    
    result.s_[0] = string( s[0].rbegin(), s[0].rend() );
    result.s_[1] = string( s[1].rbegin(), s[1].rend() );
    result.start_ = start;
    result.len_ = s[0].size() - len - start;
    result.score_ = freeEnds_[1] ? m_[iMax_][jMax_] : m_[a_.size()][b_.size()];
    
//    return make_pair( start, len );
    return result;
}

AlignResult Alignment::alignBySeed( string& a, string& b, int aStart, int bStart, int len )
{
    Alignment l( a, b, 0, aStart, 0, bStart );
    Alignment r( a, b, aStart+len, a.size()-aStart-len, bStart+len, b.size()-bStart-len );
    AlignResult lAns = l.align( false, true );
    AlignResult rAns = r.align( true, false );
    AlignResult result;
    result.s_[0] = lAns.s_[0] + a.substr( aStart, len ) + rAns.s_[0];
    result.s_[1] = lAns.s_[1] + b.substr( bStart, len ) + rAns.s_[1];
    result.lIgnore[0] = lAns.lIgnore[0];
    result.lIgnore[1] = lAns.lIgnore[1];
    result.rIgnore[0] = rAns.rIgnore[0];
    result.rIgnore[1] = rAns.rIgnore[1];
    result.start_ = lAns.start_;
    result.len_ = lAns.len_ + len + rAns.len_;
    result.score_ = lAns.score_ + rAns.score_ + ( len * l.hit_ );
    
    return result;
}

//SnpAlignment::SnpAlignment( string& a, string& b, int aStart, int aLen, int bStart, int bLen )
//: a_( a.substr( aStart, aLen ) ), b_( b.substr( bStart, bLen ) ), m_( aLen+1, vector<int>( bLen+1 ) ), pp_( aLen+1, vector<AlignPointer>( bLen+1 ) ), snp_( aLen+1, vector< pair<int, int> >( bLen+1, make_pair( -1, -1 ) ) ), scored_( false )
//{
//    
//}

SnpAlignment::BubbleAlign::BubbleAlign( SnpAlignment* parent, Bubble* b, int coord, int parentLen  )
: parent_( parent ), bubble_( b ), snps_( NULL ), snp_( NULL ), start_( b->start_-coord ), end_( b->start_+b->len_-coord )
{
    if ( b->type_ != 2 )
    {
        assert( !b->len_ );
        if ( b->type_ ) end_ = parentLen;
        if ( !b->type_ ) start_ = 0;
    }
    start_ = max( 0, start_ );
    end_ = min( end_, parentLen );
    align_ = new SnpAlignment( this, b->template_, start_ );
    align_->wrapper_ = this;
    for ( Bubble* bub : b->bubs_ )
    {
        align_->bubbles_.push_back( new BubbleAlign( align_, bub, 0, b->template_.size() ) );
    }
}

SnpAlignment::BubbleAlign::BubbleAlign( SnpAlignment* parent, SNPs* snps, int i, int coord )
: align_( NULL ), parent_( parent ), bubble_( NULL ), snps_( snps ), snp_( &snps->snps_[i] ), start_( snps->start_-coord ), end_( snps->start_+snps->len_-coord )
{
    if ( !snps->snps_[i].seq_.empty() ) align_ = new SnpAlignment( this, snps->snps_[i].seq_, start_ );
    if ( align_ ) align_->wrapper_ = this;
}

bool SnpAlignment::BubbleAlign::isThis( SnpAlignResult::BubbleAlignCoords& bac )
{
    if ( bubble_ && bac.bubble_ != bubble_ ) return false;
    if ( snps_ && bac.snps_ != snps_ ) return false;
    if ( snp_ && bac.snp_ != snp_ ) return false;
    return true;
}


SnpAlignment::SnpAlignment( string& a, string& b, int aStart, int aLen, int bStart, int bLen, vector<SNPs*> snps, vector<Bubble*>& bubbles )
: Alignment( a, b, aStart, aLen, bStart, bLen ), wrapper_( NULL ), s_( aLen+1, vector<AlignPointer>( bLen+1 ) ), snps_( snps ), parent_( NULL )
{
    coord_[0] = aStart;
    coord_[1] = bStart;
    sort( snps_.begin(), snps_.end(), []( SNPs* a, SNPs* b ){ return a->start_+a->len_ < b->start_+b->len_; } );
    hit_ = 1;
    miss_ = 2;
    gapOpen_ = 3;
    gapExt_ = 3;
    for ( SNPs* s : snps ) for ( int i = 0; i < s->snps_.size(); i++ ) bubbles_.push_back( new BubbleAlign( this, s, i, coord_[0] ) );
    for ( Bubble* b : bubbles ) bubbles_.push_back( new BubbleAlign( this, b, coord_[0], a_.size() ) );
    for ( int i = 0; i < bubbles_.size(); i++ ) if ( bubbles_[i]->start_ < 0 || bubbles_[i]->end_ > a_.size() )
    {
        
    }
}

SnpAlignment::SnpAlignment( BubbleAlign* wrapper, string& a, int coord )
: Alignment( a, wrapper->parent_->b_, 0, a.size(), 0, wrapper->parent_->b_.size() ), s_( a.size()+1, vector<AlignPointer>( wrapper->parent_->b_.size()+1 ) ), parent_( wrapper->parent_ )
{
    coord_[0] = coord;
    coord_[1] = parent_->coord_[1];
    hit_ = parent_->hit_;
    miss_ = parent_->miss_;
    gapOpen_ = parent_->gapOpen_;
    gapExt_ = parent_->gapExt_;
    freeEnds_[0] = parent_->freeEnds_[0];
    freeEnds_[1] = parent_->freeEnds_[1];
}

//SnpAlignment::SnpAlignment( SnpAlignment* base, Bubble* bubble )
//: Alignment( bubble->template_, base->b_, 0, bubble->template_.size(), 0, base->b_.size() ), parent_( base )
//{
//    coord_[0] = base->parent_ ? bubble->start_ : bubble->start_-base->coord_[0];
//    coord_[1] = base->coord_[1];
//    hit_ = base->hit_;
//    miss_ = base->miss_;
//    gapOpen_ = base->gapOpen_;
//    gapExt_ = base->gapExt_;
//    freeEnds_[0] = base->freeEnds_[0];
//    freeEnds_[1] = base->freeEnds_[1];
//    assert( false );
//}
//
//SnpAlignment::SnpAlignment( SnpAlignment* base, SNPs* snps, int i )
//: Alignment( snps->snps_[i].seq_, base->b_, 0, snps->snps_[i].seq_.size(), 0, base->b_.size() ), parent_( base )
//{
//    coord_[0] = base->parent_ ? snps->start_ : snps->start_-base->coord_[0];
//    coord_[1] = base->coord_[1];
//    hit_ = base->hit_;
//    miss_ = base->miss_;
//    gapOpen_ = base->gapOpen_;
//    gapExt_ = base->gapExt_;
//    freeEnds_[0] = base->freeEnds_[0];
//    freeEnds_[1] = base->freeEnds_[1];
//    assert( false );
//}

void SnpAlignment::score()
{
    vector<GapPointer> iMax( b_.size(), GapPointer() );
    for ( int i = 0; i < a_.size()+1; i++ ) s_[i][0].set( this, 0, 0, freeEnds_[0] || !i ? 0 : -gapOpen_ - ( i * gapExt_ ) );
    for ( int j = 0; j < b_.size()+1; j++ ) s_[0][j].set( this, 0, 0, freeEnds_[0] || !j ? 0 : -gapOpen_ - ( j * gapExt_ ) );
    for ( int j = 0; j < iMax.size(); j++ ) iMax[j].set( this, 0, s_[0][j].score_ );
    score( iMax );
}

void SnpAlignment::score( vector<GapPointer>& iMax )
{
    vector<vector<BubbleAlign*>> bubbles[2]{ vector<vector<BubbleAlign*>>( a_.size()+1 ), vector<vector<BubbleAlign*>>( a_.size()+1 ) };
    for ( BubbleAlign* ba : bubbles_ ) if ( ( freeEnds_[0] || ba->start_ >= 0 ) && ( freeEnds_[1] || ba->end_ <= a_.size() ) )
    {
        bubbles[0][max( 0, ba->start_ )].push_back( ba );
        bubbles[1][min( (int)a_.size(), ba->end_ )].push_back( ba );
    }
    sort( bubbles_.begin(), bubbles_.end(), []( BubbleAlign* a, BubbleAlign* b ){ return a->end_ < b->end_; } );
    
    vector<int> jMax( a_.size(), 0 );
    for ( int i = 0; i < a_.size(); i++ )
    {
        for ( int j = 0; j < iMax.size(); j++ ) iMax[j].len_++;
        for ( BubbleAlign* ba : bubbles[0][i] ) setBubble( ba, iMax, i );
        
//        for ( int j = 0; j < b_.size(); j++ )
//        {
//            // Simple match
//            bool matched = s_[i+1][j+1].set( this, i, j, s_[i][j].score_ + ( a_[i] == b_[j] ? hit_ : -miss_ ) );
//            // Insert a_, gap b_
//            if ( s_[i+1][j+1].update( iMax[j].align_, iMax[j].i_+1, j+1, iMax[j].score( gapOpen_, gapExt_ ) ) ) matched = false;
//            // Insert b_, gap a_
//            if ( s_[i+1][j+1].update( this, i+1, jMax[i]+1, s_[i+1][jMax[i]+1].score_ - gapOpen_ - ( ( j-jMax[i] ) * gapExt_ ) ) ) matched = false;
//            // Try match bubbles or insert from bubble
//            for ( BubbleAlign* ba : bubbles[1][i] ) getBubble( ba, iMax, i, j, matched );
//            // Reset if free ends
//            if ( freeEnds_[0] && s_[i+1][j+1].score_ < 0 ) matched = s_[i+1][j+1].set( this, i+1, j+1, 0 );
//            
//            if ( matched )
//            {
//                if ( s_[i+1][j+1].score_ >= s_[i+1][jMax[i]+1].score_ - ( ( j-jMax[i] ) * gapExt_ ) ) jMax[i] = j;
//                if ( s_[i+1][j+1].score_ >= iMax[j].score_ - ( ( iMax[j].len_-iMax[j].del_ ) * gapExt_ ) ) iMax[j].set( s_[i+1][j+1].align_, s_[i+1][j+1].i_,  s_[i+1][j+1].score_ );
//                if ( s_[i+1][j+1].score_ >= s_[ijMax_.first][ijMax_.second].score_ ) ijMax_ = make_pair( i+1, j+1 );
//            }
////            for ( BubbleAlign* ba : dels ) if ( i == ba->end_ && i >= ba->snps_->len_ ) s_[i+1][j+1].update( this, i-ba->snps_->len_, j, s_[i-ba->snps_->len_][j].score_ + ( a_[i] == b_[j] ? hit_ : -miss_ ) );
//        }
        for ( int j = 0; j < b_.size(); j++ )
        {
            // Simple match
            bool matched = s_[i+1][j+1].set( this, i, j, s_[i][j].score_ + ( a_[i] == b_[j] ? hit_ : -miss_ ) );
            // Insert a_, gap b_
            if ( s_[i+1][j+1].update( iMax[j].align_, iMax[j].i_, j+1, iMax[j].score( gapOpen_, gapExt_ ) ) ) matched = false;
            // Insert b_, gap a_
            if ( s_[i+1][j+1].update( this, i+1, jMax[i]+1, s_[i+1][jMax[i]+1].score_ - gapOpen_ - ( ( j-jMax[i] ) * gapExt_ ) ) ) matched = false;
            // Try match bubbles or insert from bubble
            for ( BubbleAlign* ba : bubbles[1][i] ) getBubble( ba, iMax, i, j, matched );
            // Reset if free ends
            if ( freeEnds_[0] && s_[i+1][j+1].score_ < 0 ) matched = s_[i+1][j+1].set( this, i+1, j+1, 0 );
            
            if ( matched )
            {
                if ( s_[i+1][j+1].score_ >= s_[i+1][jMax[i]+1].score_ - ( ( j-jMax[i] ) * gapExt_ ) ) jMax[i] = j;
                if ( s_[i+1][j+1].score_ >= iMax[j].score_ - ( ( iMax[j].len_-iMax[j].del_ ) * gapExt_ ) ) iMax[j].set( this, i+1, s_[i+1][j+1].score_ );
                if ( s_[i+1][j+1].score_ >= s_[ijMax_.first][ijMax_.second].score_ ) ijMax_ = make_pair( i+1, j+1 );
            }
//            for ( BubbleAlign* ba : dels ) if ( i == ba->end_ && i >= ba->snps_->len_ ) s_[i+1][j+1].update( this, i-ba->snps_->len_, j, s_[i-ba->snps_->len_][j].score_ + ( a_[i] == b_[j] ? hit_ : -miss_ ) );
        }
    }
    for ( int j = 0; j < iMax.size(); j++ ) iMax[j].len_++;
    scored_ = true;
}

void SnpAlignment::setBubble( BubbleAlign* ba, vector<GapPointer> iMax, int ii )
{
    if ( ba->align_ )
    {
        for ( int i = 1; i < ba->align_->a_.size()+1; i++ ) ba->align_->s_[i][0].set( s_[ii][0].align_, 0, 0, freeEnds_[0] ? 0 : min( s_[ii][0].score_, -gapOpen_ ) - ( i * gapExt_ ) );
        for ( int j = 0; j < ba->align_->b_.size()+1; j++ ) ba->align_->s_[0][j] = s_[ii][j];
        ba->align_->score( iMax );
    }
    ba->iMax_ = iMax;
}

void SnpAlignment::getBubble( BubbleAlign* ba, vector<GapPointer>& iMax, int i, int j, bool& matched )
{
    if ( ba->align_ )
    {
        int ii = ba->align_->a_.size();
        if ( s_[i+1][j+1].update( ba->align_, ii, j, ba->align_->s_[ii][j].score_ + ( a_[i] == b_[j] ? hit_ : -miss_ ) ) ) matched = true;
        if ( s_[i+1][j+1].update( ba->iMax_[j].align_, ba->iMax_[j].i_, j+1, ba->iMax_[j].score( gapOpen_, gapExt_ ) ) ) matched = false;
        if ( iMax[j].score( gapOpen_, gapExt_ ) < ba->iMax_[j].score( gapOpen_, gapExt_ ) ) iMax[j] = ba->iMax_[j];
        for ( int k = ba->align_->bubbles_.size(); k-- > 0 && ba->align_->bubbles_[k]->end_ >= ii; )
        {
            getBubble( ba->align_->bubbles_[k], iMax, i, j, matched );
        }
    }
    else
    {
        if ( s_[i+1][j+1].update( this, ba->start_, j, s_[ba->start_][j].score_ + ( a_[i] == b_[j] ? hit_ : -miss_ ) ) ) s_[i+1][j+1].deletion_ = matched = true;
    }
}

//bool SnpAlignment::getStack( vector<BubbleAlign*>& stack, SnpAlignment* base )
//{
//    SnpAlignment* cur = this;
//    while ( cur && cur != base && cur->wrapper_ )
//    {
//        stack.insert( stack.end(), cur->wrapper_ );
//        cur = cur->wrapper_->parent_;
//    }
//    return cur == base;
//}

SnpAlignment::AlignPointer SnpAlignment::getEnd( SnpAlignResult& result )
{
    AlignPointer best;
    int i = freeEnds_[1] ? ijMax_.first : a_.size(), j = freeEnds_[1] ? ijMax_.second : b_.size();
    best.set( this, i, j, s_[i][j].score_ );
    for ( BubbleAlign* ba : bubbles_ )
    {
        if ( ba->align_ && ( freeEnds_[1] || ba->end_ == a_.size() ) )
        {
            AlignPointer alt = ba->align_->getEnd( result );
            if ( best.update( alt.align_, alt.i_, alt.j_, alt.score_ ) ) result.bubbles_.clear();
        }
        else if ( ba->snp_ && ba->snp_->seq_.empty() && ba->end_ == a_.size() && best.update( this, ba->start_, j, s_[ba->start_][j].score_ ) )
        {
            result.bubbles_.clear();
            SnpAlignResult::BubbleAlignCoords bac;
            bac.start_ = bac.end_ = 0;
            bac.coord_[0] = ba->start_;
            bac.coord_[1] = ba->end_;
            bac.bubble_ = ba->bubble_;
            bac.snps_ = ba->snps_;
            bac.snp_ = ba->snp_;
            result.bubbles_.push_back( bac );
        }
    }
    return best;
}

vector<SnpAlignment::BubbleAlign*> SnpAlignment::getStack( SnpAlignment* base )
{
    vector<BubbleAlign*> stack;
    SnpAlignment* cur = this;
    while ( cur && cur != base && cur->wrapper_ )
    {
//        stack.insert( stack.end(), cur->wrapper_ );
        stack.push_back( cur->wrapper_ );
        cur = cur->wrapper_->parent_;
    }
    return stack;
}

//void SnpAlignment::scoreOld()
//{
//    for ( int i = 0; i < a_.size()+1; i++ ) p_[i][0] = -i;
//    for ( int j = 0; j < b_.size()+1; j++ ) p_[0][j] = j;
//    for ( int i = 0; i < a_.size()+1; i++ ) m_[i][0] = freeEnds_[0] || !i ? 0 : -gapOpen_ - ( i * gapExt_ );
//    for ( int j = 0; j < b_.size()+1; j++ ) m_[0][j] = freeEnds_[0] || !j ? 0 : -gapOpen_ - ( j * gapExt_ );
//    vector<int> iMax( b_.size(), 0 ), jMax( a_.size(), 0 );
//    int iSnp = 0;
//    for ( int i = 0; i < a_.size(); i++ )
//    {
//        int coord = i+coord_[0];
//        while ( iSnp < snps_.size() && snps_[iSnp]->start_+snps_[iSnp]->len_ < coord+1 ) iSnp++;
//        for ( int j = 0; j < b_.size(); j++ )
//        {
//            m_[i+1][j+1] = m_[i][j] + ( a_[i] == b_[j] ? hit_ : -miss_ );
//            p_[i+1][j+1] = 0;
//            
//            // Find the previous i value for which the current j was highest
//            int aGapLen = j-jMax[i];
//            int bGapLen = i-iMax[j];
//            int aGapScore = m_[i+1][jMax[i]+1] - gapOpen_ - ( aGapLen * gapExt_ );
//            int bGapScore = m_[iMax[j]+1][j+1] - gapOpen_ - ( bGapLen * gapExt_ );
//            if ( aGapScore > m_[i+1][j+1] )
//            {
//                m_[i+1][j+1] = aGapScore;
//                p_[i+1][j+1] = aGapLen;
//            }
//            if ( bGapScore > m_[i+1][j+1] )
//            {
//                m_[i+1][j+1] = bGapScore;
//                p_[i+1][j+1] = -bGapLen;
//            }
//            
//            for ( int k = iSnp; k < snps_.size() && snps_[k]->start_+snps_[k]->len_ == coord+1; k++ )
//            {
//                for ( int kk = 0; kk < snps_[k]->snps_.size(); kk++ )
//                {
//                    int m, p = 0;
//                    if ( !snps_[k]->len_ )
//                    {
//                        assert( false );
//                        assert( snps_[k]->snps_[kk].seq_.size() == 1 );
//                        // gap in a
//                        m = m_[i+1][j] + ( snps_[k]->snps_[kk].seq_[0] == b_[i] ? hit_ : -miss_ );
//                        p = 1;
//                    }
//                    else if ( snps_[k]->snps_[kk].seq_.size() == 1 )
//                    {
//                        if ( snps_[k]->snps_[kk].seq_[0] != b_[j] ) continue;
//                        m = m_[i][j] + hit_;
//                    }
//                    else if ( snps_[k]->snps_[kk].seq_.empty() )
//                    {
//                        // gap in b
//                        m = m_[i][j+1];
//                        p = -1;
//                    }
//                    else
//                    {
//                        assert( false );
//                        // impossible!
//                    }
//                    
//                    if ( m > m_[i+1][j+1] )
//                    {
//                        assert( !p || m < 1 );
//                        m_[i+1][j+1] = m;
//                        p_[i+1][j+1] = p;
//                        snp_[i+1][j+1] = make_pair( k, kk );
//                    }
//                }
//            }
//            
//            if ( freeEnds_[0] && m_[i+1][j+1] < 0 )
//            {
//                m_[i+1][j+1] = 0;
//                p_[i+1][j+1] = 0;
//                snp_[i+1][j+1] = make_pair( -1, -1 );
//            }
//            
//            if ( p_[i+1][j+1] == 0 )
//            {
//                if ( m_[i+1][j+1] >= m_[i+1][jMax[i]+1] - ( aGapLen * gapExt_ ) ) jMax[i] = j;
//                if ( m_[i+1][j+1] >= m_[iMax[j]+1][j+1] - ( bGapLen * gapExt_ ) ) iMax[j] = i;
//                if ( m_[i+1][j+1] > m_[iMax_][jMax_])
//                {
//                    iMax_ = i+1;
//                    jMax_ = j+1;
//                }
//            }
//        }
//    }
//    scored_ = true;
//}

SnpAlignResult SnpAlignment::align( bool lAnchored, bool rAnchored )
{
    if ( lAnchored == freeEnds_[0] || rAnchored == freeEnds_[1] ) scored_ = false;
    freeEnds_[0] = !lAnchored;
    freeEnds_[1] = !rAnchored;
//    SnpAlignResult oldResult = alignOld( lAnchored, rAnchored );
    scored_ = false;
    if ( !scored_ ) score();
    
    SnpAlignResult result;
    AlignPointer start = getEnd( result );
    SnpAlignment* cur = start.align_;
    int i = start.i_, j = start.j_;
    
//    SnpAlignment* cur = this;
//    if ( freeEnds_[1] ) for ( pair<SnpAlignment*, int> b : findEnd() )
//    {
//        result.s_[0].insert( result.s_[0].end(), cur->a_.rbegin(), cur->a_.rbegin() + b.second );
//        cur = b.first;
//    }
//    
//    int i = freeEnds_[1] ? cur->ijMax_.first : cur->a_.size(), j = freeEnds_[1] ? cur->ijMax_.second : cur->b_.size();
//    result.s_[0].insert( result.s_[0].end(), cur->a_.rbegin(), cur->a_.rbegin() + ( cur->a_.size()-i ) );
//    result.s_[1].insert( result.s_[1].end(), b_.rbegin(), b_.rbegin() + ( b_.size()-j ) );
//    if ( result.s_[1].size() > result.s_[0].size() ) result.s_[0].insert( 0, string( result.s_[1].size()-result.s_[0].size(), '-' ) );
//    if ( result.s_[0].size() > result.s_[1].size() ) result.s_[1].insert( 0, string( result.s_[0].size()-result.s_[1].size(), '-' ) );
//    assert( result.s_[0].size() == result.s_[1].size() );
    
//    vector<BubbleAlign*> stack;
//    cur->getStack( stack, this );
    
    AlignPointer* ap = &cur->s_[i][j];
    while ( ap && ap != &s_[0][0] )
    {
        bool terminated = ap->align_ == cur && ap->i_ == i && ap->j_ == j;
        if ( terminated ) ap->set( this, 0, 0, 0 );
        if ( terminated || !i || !j ) result.len_ = result.s_[0].size();
        SnpAlignment* base = ap->align_ == cur ? cur : this;
        
        // Set bubble stacks
        vector<BubbleAlign*> stacks[2];
        if ( ap->align_ != base ) stacks[0] = ap->align_->getStack( base );
        if ( cur != base ) stacks[1] = cur->getStack( this );
        while ( !stacks[0].empty() && !stacks[1].empty() && stacks[0].back()->align_ == stacks[1].back()->align_ ) for ( int k : { 0, 1 } ) stacks[k].pop_back();
        for ( int k : { 0, 1 } ) if ( !stacks[k].empty() ) base = stacks[k].back()->parent_;
//        while ( !stacks[0].empty() && !stacks[1].empty() && stacks[0][0]->align_ == stacks[1][0]->align_ ) for ( int k : { 0, 1 } ) stacks[k].erase( stacks[k].begin() );
//        for ( int k : { 0, 1 } ) if ( !stacks[k].empty() ) base = stacks[k][0]->parent_;
        for ( int k : { 0, 1 } ) if ( !stacks[k].empty() ) assert( base == stacks[k].back()->parent_ );
        
        int len[2]{ (int)result.s_[0].size(), (int)result.s_[1].size() };
        assert( len[0] == len[1] );
        assert( result.s_[0].size() == result.s_[1].size() );
        
        // Descend bubbles
        for ( int k = 0; k < stacks[1].size(); k++ ) stacks[1][k]->align_->fill( result, i, 0, ap->j_ == j || ap->deletion_, true, false );
        // Fill base
        base->fill( result, i, stacks[0].empty() ? ap->i_ : stacks[0].back()->end_, ap->j_ == j || ap->deletion_, false, false );
        // Ascend bubbles
        for ( int k = stacks[0].size(); k-- > 0; ) stacks[0][k]->align_->fill( result, i, k ? stacks[0][k-1]->end_ : ap->i_, ap->j_ == j || ap->deletion_, false, true );
        
        // Fill b
        while ( ap->j_ < j ) result.s_[1] += b_[--j];
        if ( result.s_[1].size() > result.s_[0].size() && !terminated ) assert( result.s_[0].size() == len[0] );
        // Add gaps
        for ( int k : { 0, 1 } ) if ( result.s_[k].size() < result.s_[!k].size() ) result.s_[k] += string( result.s_[!k].size()-result.s_[k].size(), '-' );
        
        cur = ap->align_;
        ap = &cur->s_[ap->i_][ ap->j_];
    }
    
    finish( result, start );
//    assert( oldResult.s_[0] == result.s_[0] );
//    assert( oldResult.s_[1] == result.s_[1] );
    
    return result;
}

void SnpAlignment::add( SnpAlignResult& result, int i )
{
    if ( wrapper_ )
    {
        vector<BubbleAlign*> stack = getStack( NULL );
        vector<SnpAlignResult::BubbleAlignCoords>* cur = &result.bubbles_;
        for ( int k = stack.size(); k-- > 0; )
        {
            if ( cur->empty() || !stack[k]->isThis( cur->back() ) )
            {
                SnpAlignResult::BubbleAlignCoords bac;
                bac.start_ = bac.end_ = result.s_[0].size();
                bac.coord_[0] = stack[k]->start_;
                bac.coord_[1] = stack[k]->end_;
                bac.bubble_ = stack[k]->bubble_;
                bac.snps_ = stack[k]->snps_;
                bac.snp_ = stack[k]->snp_;
                cur->push_back( bac );
            }
            assert( cur->back().end_ == result.s_[0].size() );
            cur->back().end_++;
            cur = &cur->back().bubbles_;
        }
    }
    result.s_[0] += a_[i];
}

void SnpAlignment::fill( SnpAlignResult& result, int& i, int iEnd, bool deletion, bool descending, bool ascending )
{
    if ( ascending ) i = a_.size();
    if ( deletion )
    {
        for ( BubbleAlign* ba : bubbles_ ) if ( ba->snp_ && ba->snp_->seq_.empty() && iEnd <= ba->start_ && ba->end_ <= i )
        {
            while ( ba->end_ < i ) add( result, --i );
            SnpAlignResult::BubbleAlignCoords bac;
            bac.start_ = bac.end_ = result.s_[0].size();
            bac.coord_[0] = ba->start_;
            bac.coord_[1] = ba->end_;
            bac.bubble_ = ba->bubble_;
            bac.snps_ = ba->snps_;
            bac.snp_ = ba->snp_;
            result.bubbles_.push_back( bac );
            i = ba->start_;
            break;
        }
    }
    while ( iEnd < i ) add( result, --i );
//    while ( iEnd < i ) result.s_[0] += a_[--i];
    if ( descending ) i = wrapper_->start_;
}

void SnpAlignment::finish( SnpAlignResult& result, AlignPointer& start )
{
    result.reverse();
    
    int i = start.i_, len = result.start_ + result.len_;
    for ( BubbleAlign* ba : start.align_->getStack( this ) )
    {
        while ( i < ba->align_->a_.size() ) ba->align_->add( result, i++ );
        i = ba->end_;
    }
    if ( !result.bubbles_.empty() && result.bubbles_.back().coord_[0] == i && result.bubbles_.back().snp_ && result.bubbles_.back().snp_->seq_.empty())
    {
        i = result.bubbles_.back().coord_[1];
    }
    result.s_[0] += a_.substr( i );
    result.s_[1] += b_.substr( start.j_ );
    for ( int s : { 0, 1 } ) result.rIgnore[s] = result.s_[s].size() - len;
    for ( int s : { 0, 1 } ) if ( result.s_[s].size() < result.s_[!s].size() ) result.s_[s] += string( result.s_[!s].size()-result.s_[s].size(), '-' );
    result.score_ = start.score_;
}

//bool SnpAlignment::build( SnpAlignment& result, AlignPointer*& ap, int& i, int& j )
//{
//    while ( ap->align_ == this && ( i > 0 || j > 0 ) )
//    {
//        if ( s_[i][j].i_ == i && s_[i][j].j_ == j ) assert( freeEnds_[0] );
//        if ( s_[i][j].i_ == i && s_[i][j].j_ == j ) return false;
//        
//        assert( ( s_[i][j].i_ == i-1 && s_[i][j].j_ == j-1 ) || s_[i][j].i_ == i || s_[i][j].j_ == j );
//        if ( s_[i][j].i_ == i ) for ( int k = 0; k < j-s_[i][j].j_; k++ ) result.s_[0] += '-';
//        if ( s_[i][j].j_ == j ) for ( int k = 0; k < i-s_[i][j].i_; k++ ) result.s_[1] += '-';
//        while ( i > s_[i][j].i_ ) result.s_[0] += a_[--i];
//        while ( j > s_[i][j].j_ ) result.s_[1] += b_[--j];
//        assert( result.s_[0].size() == result.s_[1].size() );
//        ap = s_[i][j];
//        src = this;
//    }
//    return true;
//}

//SnpAlignment::AlignPointer* SnpAlignment::build( SnpAlignment& result, vector<BubbleAlign*> stack, int i, int& j )
//{
//    AlignPointer* ap;
//    if ( !stack.empty() )
//    {
//        ap = stack[0]->align_->build( result, vector<BubbleAlign*>( stack.begin()+1, stack.end() ), i, j );
//        i = stack[0]->start_;
//    }
//    
//    while ( ap )
//    {
//        if ( ap->align_ != this )
//        {
//            vector<BubbleAlign*> stack;
//            if ( !ap->align_->getStack( stack, this ) )
//            {
//                assert( !i );
//                while ( i > 0 ) result.s_[0].push_back( a_[--i] );
//                return ap;
//            }
//            
//            // Ascend bubble stack
//            while ( stack[0]->end_ < i ) result.s_[0].push_back( a_[--i] );
//            ap = stack[0]->align_->build( result, vector<BubbleAlign*>( stack.begin()+1, stack.end() ), i, j );
//            i = stack[0]->start_;
//        }
//        else
//        {
//            ap = s_[ap->i_][ap->j_];
//        }
//    }
//    
//    return &s_[i][j];
//}

vector<pair<SnpAlignment*, int>> SnpAlignment::findEnd()
{
    SnpAlignment* best = this;
    vector<pair<SnpAlignment*, int>> bestList;
    
    for ( BubbleAlign* ba : bubbles_ ) if ( ba->align_ )
    {
        vector<pair<SnpAlignment*, int>> bubbleList = ba->align_->findEnd();
        bubbleList.insert( bubbleList.begin(), make_pair( ba->align_, a_.size()-ba->end_ ) );
        SnpAlignment* cur = bubbleList.back().first;
        if ( cur->s_[cur->ijMax_.first][cur->ijMax_.second].score_ <= best->s_[best->ijMax_.first][best->ijMax_.second].score_ ) continue;
        assert( false );
        best = cur;
        bubbleList = bestList;
    }
    
    return bestList;
}

//SnpAlignResult SnpAlignment::alignOld( bool lAnchored, bool rAnchored )
//{
//    if ( lAnchored == freeEnds_[0] || rAnchored == freeEnds_[1] ) scored_ = false;
//    freeEnds_[0] = !lAnchored;
//    freeEnds_[1] = !rAnchored;
//    if ( !scored_ ) scoreOld();
//    
//    SnpAlignResult result;
//    result.lIgnore[0] = result.lIgnore[1] = result.rIgnore[0] = result.rIgnore[1] = 0;
//    
//    string s[2]{ result.s_[0], result.s_[1] };
//    int i = a_.size(), j = b_.size();
//    if ( freeEnds_[1] )
//    {
//        result.rIgnore[0] = a_.size()-iMax_;
//        result.rIgnore[1] = b_.size()-jMax_;
//        while ( i-iMax_ > j-jMax_ ) s[0] += a_[--i];
//        while ( j-jMax_ > i-iMax_ ) s[1] += b_[--j];
//        if ( s[0].size() > s[1].size() ) s[1] += string( s[0].size(), '-');
//        if ( s[1].size() > s[0].size() ) s[0] += string( s[1].size(), '-');
//        while ( i > iMax_ ) s[0] += a_[--i];
//        while ( j > jMax_ ) s[1] += b_[--j];
//    }
//    
//    int start = 0, len = s[0].size();
//    while ( i > 0 || j > 0 )
//    {
//        assert( i >= 0 && j >= 0 );
//        if ( !i || !j || ( freeEnds_[0] && p_[i][j] == 0 && m_[i][j] == 0 && !( i && j && m_[i-1][j-1] == miss_ ) ) )
//        {
//            if ( freeEnds_[0] ) result.lIgnore[0] = i;
//            if ( freeEnds_[0] ) result.lIgnore[1] = j;
//            if ( freeEnds_[0] ) start = s[0].size();
//            while ( i > 0 ) s[0] += a_[--i];
//            while ( j > 0 ) s[1] += b_[--j];
//            if ( s[0].size() > s[1].size() ) s[1] += string( s[0].size()-s[1].size(), '-');
//            if ( s[1].size() > s[0].size() ) s[0] += string( s[1].size()-s[0].size(), '-');
//            if ( freeEnds_[0] ) start = s[0].size() - start;
//            break;
//        }
//        
//        if ( snp_[i][j].first >= 0 )
//        {
//            result.snps_.push_back( make_pair( snps_[snp_[i][j].first], snp_[i][j].second ) );
//        }
//        
//        if ( p_[i][j] == 0 )
//        {
//            s[0] += a_[--i];
//            s[1] += b_[--j];
//        }
//        else if ( p_[i][j] < 0 )
//        {
//            for ( int k = p_[i][j]; k++ < 0; )
//            {
//                s[0] += a_[--i];
//                s[1] += '-';
//            }
//        }
//        else
//        {
//            for ( int k = p_[i][j]; k-- > 0; )
//            {
//                s[0] += '-';
//                s[1] += b_[--j];
//            }
//        }
//    }
//    
//    result.s_[0] = string( s[0].rbegin(), s[0].rend() );
//    result.s_[1] = string( s[1].rbegin(), s[1].rend() );
//    result.start_ = start;
//    result.len_ = s[0].size() - len - start;
//    result.score_ = freeEnds_[1] ? m_[iMax_][jMax_] : m_[a_.size()][b_.size()];
//    
//    return result;
//}
