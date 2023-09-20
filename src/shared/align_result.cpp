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

#include "align_result.h"
#include "bubble.h"
#include <cassert>
#include <algorithm>

AlignResult::AlignResult()
{
    lIgnore[0] = lIgnore[1] = rIgnore[0] = rIgnore[1] = score_ = 0;
    start_ = len_ = -1;
}

int AlignResult::getSeqLen( int i )
{
    int len = 0;
    for ( int j = 0; j < s_[i].size(); j++ ) if ( s_[i][j] != '-' ) len++;
    return len;
}

void AlignResult::trimFromEnd( int sIndex, int trimLen, bool drxn )
{
    int len = 0;
    for ( int i = drxn ? s_[0].size()-1 : 0; i < s_[0].size() && i >= 0 && trimLen > 0; drxn ? i-- : i++ )
    {
        if ( s_[sIndex][i] != '-' ) trimLen--;
        for ( int s : { 0, 1 } ) if ( s_[sIndex][i] != '-' && ( drxn ? rIgnore[s] : lIgnore[s] ) > 0 ) ( drxn ? rIgnore[s] : lIgnore[s] )--;
        len++;
    }
    for ( int s : { 0, 1 } ) s_[s].erase( drxn ? s_[s].end()-len : s_[s].begin(), drxn ? s_[s].end() : s_[s].begin()+len );
    if ( start_ > s_[0].size() ) start_ = s_[0].size();
    if ( drxn ) len_ = min( len_, (int)s_[0].size()-start_ );
    if ( !drxn && len > start_ ) len_ -= min( len_, len-start_ );
    if ( !drxn ) start_ = max( start_-len, 0 );
}

void SnpAlignResult::BubbleAlignCoords::reverse( vector<BubbleAlignCoords>& bubbles, int base )
{
    std::reverse( bubbles.begin(), bubbles.end() );
    for ( BubbleAlignCoords& bc : bubbles )
    {
        int len = bc.end_ - bc.start_;
        bc.start_ = base - bc.end_;
        bc.end_ = bc.start_ + len;
        BubbleAlignCoords::reverse( bc.bubbles_, base );
    }
}

bool SnpAlignResult::BubbleAlignCoords::isDeletion()
{
    return snp_ ? snp_->seq_.empty() : bubble_->template_.empty();
}

void SnpAlignResult::reverse()
{
    for ( int k : { 0, 1 } ) std::reverse( s_[k].begin(), s_[k].end() );
    
    if ( len_ < 0 ) len_ = s_[0].size();
    start_ = s_[0].size() - len_;
    for ( int k = 0; k < start_; k++ ) for ( int s : { 0, 1 } ) if ( s_[s][k] != '-' ) lIgnore[s]++;
    
    BubbleAlignCoords::reverse( bubbles_, s_[0].size() );
}