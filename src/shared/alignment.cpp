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
                if ( m_[i+1][j+1] >= m_[i+1][jMax[i]+1] - ( bGapLen * gapExt_ ) ) jMax[i] = j;
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

pair<int,int> Alignment::align( bool lAnchored, bool rAnchored )
{
    if ( lAnchored == freeEnds_[0] || rAnchored == freeEnds_[1] ) scored_ = false;
    freeEnds_[0] = !lAnchored;
    freeEnds_[1] = !rAnchored;
    if ( !scored_ ) score();
//    debug();
    
    string s[2];
    int i = a_.size(), j = b_.size();
    if ( freeEnds_[1] )
    {
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
        if ( freeEnds_[0] && p_[i][j] == 0 && m_[i][j] == 0 && !( i && j && m_[i][j] == miss_ ) )
        {
            start = s[0].size();
            while ( i > 0 ) s[0] += a_[--i];
            while ( j > 0 ) s[1] += b_[--j];
            if ( s[0].size() > s[1].size() ) s[1] += string( s[0].size()-s[1].size(), '-');
            if ( s[1].size() > s[0].size() ) s[0] += string( s[1].size()-s[0].size(), '-');
            start = s[0].size() - start;
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
    len = s[0].size() - len - start;
    reverse( s[0].begin(), s[0].end() );
    reverse( s[1].begin(), s[1].end() );
    cout << s[0] << endl << s[1] << endl ;
    return make_pair( start, len );
}