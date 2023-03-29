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

SnpAlignment::SnpAlignment( string& a, string& b, int aStart, int aLen, int bStart, int bLen, vector<SNPs*>& snps )
: Alignment( a, b, aStart, aLen, bStart, bLen ), snp_( aLen+1, vector< pair<int, int> >( bLen+1, make_pair( -1, -1 ) ) ), snps_( snps )
{
    coord_[0] = aStart;
    coord_[1] = bStart;
    sort( snps_.begin(), snps_.end(), []( SNPs* a, SNPs* b ){ return a->start_+a->len_ < b->start_+b->len_; } );
    hit_ = 1;
    miss_ = 3;
    gapOpen_ = 3;
    gapExt_ = 3;
}

void SnpAlignment::score()
{
    for ( int i = 0; i < a_.size()+1; i++ ) p_[i][0] = -i;
    for ( int j = 0; j < b_.size()+1; j++ ) p_[0][j] = j;
    for ( int i = 0; i < a_.size()+1; i++ ) m_[i][0] = freeEnds_[0] || !i ? 0 : -gapOpen_ - ( i * gapExt_ );
    for ( int j = 0; j < b_.size()+1; j++ ) m_[0][j] = freeEnds_[0] || !j ? 0 : -gapOpen_ - ( j * gapExt_ );
    vector<int> iMax( b_.size(), 0 ), jMax( a_.size(), 0 );
    int iSnp = 0;
    for ( int i = 0; i < a_.size(); i++ )
    {
        int coord = i+coord_[0];
        while ( iSnp < snps_.size() && snps_[iSnp]->start_+snps_[iSnp]->len_ < coord+1 ) iSnp++;
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
            
            for ( int k = iSnp; k < snps_.size() && snps_[k]->start_+snps_[k]->len_ == coord+1; k++ )
            {
                for ( int kk = 0; kk < snps_[k]->snps_.size(); kk++ )
                {
                    int m, p = 0;
                    if ( !snps_[k]->len_ )
                    {
                        assert( false );
                        assert( snps_[k]->snps_[kk].seq_.size() == 1 );
                        // gap in a
                        m = m_[i+1][j] + ( snps_[k]->snps_[kk].seq_[0] == b_[i] ? hit_ : -miss_ );
                        p = 1;
                    }
                    else if ( snps_[k]->snps_[kk].seq_.size() == 1 )
                    {
                        if ( snps_[k]->snps_[kk].seq_[0] != b_[j] ) continue;
                        m = m_[i][j] + hit_;
                    }
                    else if ( snps_[k]->snps_[kk].seq_.empty() )
                    {
                        // gap in b
                        m = m_[i][j+1];
                        p = -1;
                    }
                    else
                    {
                        assert( false );
                        // impossible!
                    }
                    
                    if ( m > m_[i+1][j+1] )
                    {
                        assert( !p || m < 1 );
                        m_[i+1][j+1] = m;
                        p_[i+1][j+1] = p;
                        snp_[i+1][j+1] = make_pair( k, kk );
                    }
                }
            }
            
            if ( freeEnds_[0] && m_[i+1][j+1] < 0 )
            {
                m_[i+1][j+1] = 0;
                p_[i+1][j+1] = 0;
                snp_[i+1][j+1] = make_pair( -1, -1 );
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

SnpAlignResult SnpAlignment::align( bool lAnchored, bool rAnchored )
{
    if ( lAnchored == freeEnds_[0] || rAnchored == freeEnds_[1] ) scored_ = false;
    freeEnds_[0] = !lAnchored;
    freeEnds_[1] = !rAnchored;
    if ( !scored_ ) score();
    
    SnpAlignResult result;
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
            break;
        }
        
        if ( snp_[i][j].first >= 0 )
        {
            result.snps_.push_back( make_pair( snps_[snp_[i][j].first], snp_[i][j].second ) );
        }
        
        if ( p_[i][j] == 0 )
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
    
    result.s_[0] = string( s[0].rbegin(), s[0].rend() );
    result.s_[1] = string( s[1].rbegin(), s[1].rend() );
    result.start_ = start;
    result.len_ = s[0].size() - len - start;
    result.score_ = freeEnds_[1] ? m_[iMax_][jMax_] : m_[a_.size()][b_.size()];
    
    return result;
}
