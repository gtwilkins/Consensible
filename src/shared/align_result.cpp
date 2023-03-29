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