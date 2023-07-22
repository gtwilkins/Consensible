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

#include "match.h"
#include "read.h"
#include "target.h"
#include "alignment.h"
#include <cassert>
#include <iostream>

Match::Match( Target* tar, MappedRead* read, AlignResult& ar, int tarCoord )
: tar_( tar ), read_( read ), score_( ar.score_ )
{
    anchors_[0] = ar.lIgnore[1];
    anchors_[1] = read->seq_.size() - ar.rIgnore[1] - 1;
    for ( int i = 0; i < ar.s_[0].size(); i++ )
    {
        if ( ar.s_[1][i] != '-' ) coords_.push_back( tarCoord );
        if ( ar.s_[0][i] != '-' ) tarCoord++;
    }
    read_->matches_.push_back( this );
}

//Match::Match( Target* tar, MappedRead* read, string t, string r, vector<pair<int,int>> &anchors, int tarCoord )
//: tar_( tar ), read_( read ), score_( 0 )
//{
//    assert( t.size() == r.size() );
//    assert( !anchors.empty() );
//    anchors_[0] = anchors[0].first;
//    anchors_[1] = anchors.back().second;
//    for ( int i = 0; i < t.size(); i++ )
//    {
//        if ( i == anchors[0].first ) anchors_[0] = coords_.size();
//        if ( i == anchors.back().second ) anchors_[1] = coords_.size();
//        if ( i >= anchors[0].first && i <= anchors.back().second ) score_ += t[i] == r[i] ? 1 : -2;
//        if ( r[i] != '-' ) coords_.push_back( tarCoord );
//        if ( t[i] != '-' ) tarCoord++;
//    }
//    assert( coords_.size() == read->seq_.size() );
//    read_->matches_.push_back( this );
//    
//    int curI = 0, curLen = 0, curScore = 0;
//    int bestI = 0, bestLen = 0, bestScore = 0;
//    for ( int i = 0; i < coords_.size(); i++ )
//    {
//        int gap = i ? coords_[i]-1-coords_[i-1] : 0, match = read->seq_[i] == tar->seq_[coords_[i]] ? 1 : -4;
//        for ( ; i+1 < coords_.size() && coords_[i] == coords_[i+1]; i++ ) gap++;
//        int gapScore = ( gap ? -4 : 0 ) - ( gap * 4 ), matchScore = i && read->seq_[i] == read->seq_[i-1] ? min( 0, match ) : match;
//        
//        if ( curScore + gapScore < 0 )
//        {
//            curI = i;
//            curScore = curLen = 0;
//            gapScore = 0;
//        }
//        curScore += gapScore + matchScore;
//        curLen++;
//        
//        if ( curScore < 0 )
//        {
//            curI = i+1;
//            curScore = curLen = 0;
//        }
//        else if ( curScore == bestScore ? curLen > bestLen : curScore > bestScore )
//        {
//            bestI = curI;
//            bestScore = curScore;
//            bestLen = curLen;
//        }
//    }
//    assert( bestLen > 15 );
//    anchors_[0] = bestLen > 10 ? bestI : coords_.size();
//    anchors_[1] = bestLen > 10 ? bestI + bestLen - 1 : 0;
//}

void Match::anchorCoords( Match* match, vector< pair<int,int>>& coords )
{
    if ( coords.empty() ) return;
    int limits[2]{ coords[0].first, coords.back().first }, tlimits[2]{ (int)tar_->seq_.size(), 0 };
    int drxn = anchors_[1] < limits[0];
    assert( drxn || limits[1] < anchors_[0] );
    
    bool changed = false;
    for ( int i = 0; i < coords.size(); i++ )
    {
        int ins = coords[i].second < 0;
        if ( ins ) while ( ++i < coords.size() && coords[i].second < 0 ) ins++;
        assert( !ins );
        if ( i >= coords.size() ) assert( false );
        if ( coords_[coords[i].first] != match->coords_[coords[i].second] ) changed = true;
        coords_[coords[i].first] = match->coords_[coords[i].second];
        if ( coords_[coords[i].first] < tlimits[0] ) tlimits[0] = coords_[coords[i].first];
        if ( coords_[coords[i].first] > tlimits[1] ) tlimits[1] = coords_[coords[i].first];
        if ( match->anchors_[0] <= coords[i].second && coords[i].second <= match->anchors_[1] )
        {
            if ( coords[i].first < anchors_[0] ) anchors_[0] = coords[i].first;
            if ( anchors_[1] < coords[i].first ) anchors_[1] = coords[i].first;
        }
    }
    assert( tlimits[1] );
    
    int len = drxn ? coords_.size() - limits[1] - 1 : limits[0];
    if ( !changed || !len ) return;
    if ( drxn ? coords_[limits[1]] < coords_[limits[1]+1] : coords_[limits[0]-1] < coords_[limits[0]] ) return;
//    if ( drxn && !tlimits[1] ) tlimits[1] = coords_[limits[1]+1]-1;
//    if ( !drxn && tlimits[0] == tar_->seq_.size() ) tlimits[0] = coords_[limits[0]-1]+1;
    
    assert( false );
    int tlen = min( drxn ? (int)tar_->seq_.size() - tlimits[1] - 1 : tlimits[0], 30 + len );
    AlignResult result = Alignment( read_->seq_, tar_->seq_, drxn ? limits[1]+1 : 0, len, drxn ? tlimits[1]+1 : tlimits[0]-tlen, tlen ).align( drxn, !drxn );
    assert( false );
    
}

int Match::getAnchor( int drxn )
{
    return coords_[anchors_[drxn]];
}

vector<pair<int, int>> Match::getGaps()
{
    vector<pair<int, int>> gaps;
    for ( int i = 1; i < coords_.size(); i++ ) if ( coords_[i] == coords_[i-1] )
    {
        int gap = 1;
        for ( ; i+1 < coords_.size() && coords_[i] == coords_[i+1]; i++ ) gap++;
        gaps.push_back( make_pair( coords_[i], gap ) );
    }
    
    return gaps;
}

bool Match::isBridge( Match* l, Match* r, AlignResult& result, int competingCons )
{
    int curI = 0, curLen = 0, curScore = 0;
    int bestI = 0, bestLen = 0, bestScore = -100;
    int gapD = 0, gapLen = 0;
    int counts[2]{0};
    int minMismatch = abs( l->getAnchor( 1 ) - r->getAnchor( 0 ) ) + l->coords_.size() + r->coords_.size();
    for ( int i = 0; i < result.s_[0].size(); i++ )
    {
        int gapScore = gapLen && result.s_[gapD][i] != '-' ? -( gapLen+1 ) * 4 : 0;
        if ( gapScore ) gapLen = 0;
        if ( !gapLen ) for ( int d : { 0, 1 } ) if ( result.s_[d][i] == '-' && ( gapLen = 1 ) ) gapD = d;
        
        int matchScore = gapLen ? 0 : ( result.s_[0][i] == result.s_[1][i] ? 1 : -4 );
        if ( counts[0] && counts[1] && matchScore > 0 && ( 
                l->read_->seq_[counts[0]] == l->read_->seq_[counts[0]-1] || 
                r->read_->seq_[counts[1]] == r->read_->seq_[counts[1]-1] ) ) matchScore = 0;
        
        if ( curScore + gapScore < 0 && counts[1] <= r->anchors_[0] )
        {
            curI = i;
            curScore = curLen = 0;
            gapScore = 0;
        }
        curScore += gapScore + matchScore;
        curLen++;
        if ( curScore < 0 && counts[1] < r->anchors_[0] )
        {
            curI = i+1;
            curScore = curLen = 0;
        }
        
        if ( counts[0] >= l->anchors_[1] && ( curScore == bestScore ? curLen > bestLen : curScore > bestScore ) && !gapLen )
        {
            bestI = curI;
            bestScore = curScore;
            bestLen = curLen;
        }
        
        if ( minMismatch >= 20 && result.s_[0][i] == result.s_[1][i] )
        {
            int lCoord = counts[0] <= l->anchors_[1] ? l->coords_[counts[0]] : l->getAnchor( 1 ) + ( counts[0] - l->anchors_[1] );
            int rCoord = r->anchors_[0] <= counts[1] ? r->coords_[counts[1]] : r->getAnchor( 0 ) - ( r->anchors_[0] - counts[1] );
            minMismatch = min( minMismatch, abs( lCoord-rCoord ) );
        }
        
        for ( int d : { 0, 1 } ) if ( result.s_[d][i] != '-' ) counts[d]++;
    }
    
    int score = bestScore+bestLen;
    int competingPenalty = min( 10, max( 0, competingCons-5 ) );
    int mismatchPenalty = min( 10, max( 0, ( minMismatch-10 ) / 10 ) );
    if ( bestScore <= 0 || bestLen < 10 ) return false;
    return ( bestScore+bestLen ) >= ( 20 + competingPenalty + mismatchPenalty );
}

bool Match::isInsert( int coord )
{
    return coord+1 < coords_.size() && coords_[coord] == coords_[coord+1];
}

int Match::size()
{
    return read_->seq_.size();
}

//void Match::updateCoords( Match* match, vector<int> coords, int start )
//{
//    vector<pair<int, int>> diffs[2];
//    for ( int i = 0; i < coords.size(); i++ )
//    {
//        int j = i+start, k = coords[i];
//        int gap[2]{ j+1 >= coords.size() ? 1 : coords_[j+1]-coords_[j], k+1 >= match->coords_.size() ? 1 : match->coords_[k+1]-match->coords_[k] };
//        bool matched = coords_[j] == match->coords_[k] && ( gap[0] == gap[1] || ( gap[0] && gap[1] ) );
//        if ( !diffs[0].empty() && matched )
//        {
//            updateCoords( match, diffs[0], diffs[1] );
//            for ( int d : { 0, 1 } ) diffs[d].clear();
//        }
//        else if ( !matched )
//        {
//            diffs[0].push_back( make_pair( j, k ) );
//            if ( !diffs[1].empty() ) for ( int kk = diffs[1].back().first+1; kk < k; kk++ ) assert( false );
//            if ( !diffs[1].empty() ) for ( int kk = diffs[1].back().first+1; kk < k; kk++ ) diffs[1].push_back( make_pair( kk, j ) );
//            diffs[1].push_back( make_pair( k, j ) );
//        }
//    }
//    
////    int cur = 0, len = 0;
////    vector<pair<int, int>> diffs;
////    for ( int i = 0; i < coords.size(); i++ )
////    {
////        int j = i+start, k = coords[i];
////        int gap[2]{ j+1 >= coords.size() ? 1 : coords_[j+1]-coords_[j], k+1 >= match->coords_.size() ? 1 : match->coords_[k+1]-match->coords_[k] };
////        bool matched = coords_[j] == match->coords_[k] && ( gap[0] == gap[1] || ( gap[0] && gap[1] ) );
////        if ( len && matched )
////        {
////            diffs.push_back( make_pair( cur, len ) );
////            len = 0;
////        }
////        else if ( !matched )
////        {
////            if ( !len ) cur = i;
////            len++;
////        }
//////        cout << read_->seq_.substr( 0 ) << endl;
//////        cout << tar_->seq_.substr( coords_[0], 100 ) << endl;
//////        cout << match->read_->seq_.substr( 0 ) << endl;
//////        cout << tar_->seq_.substr( match->coords_[0], 100 ) << endl;
////    }
////    for ( pair<int,int> diff : diffs )
////    {
////        int prv[2]{ diff.first+start ? coords_[diff.first+start-1]: coords_[0]-1, coords[diff.first] ? match->coords_[coords[diff.first]-1] : match->coords_[0]-1 };
////        int prvK = coords[diff.first]-1;
////        for ( int i = diff.first; i < diff.first+diff.second; i++ )
////        {
////            int j = i+start, k = coords[i];
////            int cur[2]{ coords_[j], match->coords_[k] };
////            int nxt[2]{ j+1 < coords_.size() ? coords_[j+1] : coords_[j]+1, k+1 < match->coords_.size() ? match->coords_[k+1] : match->coords_[k]+1 };
////            int curGap[2]{ prv[0] < cur[0] ? cur[0]-1-prv[0] : nxt[0]-1-cur[0], prv[1] < cur[1] ? cur[1]-1-prv[1] : nxt[1]-1-cur[1]  };
////            int altGap[2]{0};
////            assert( !curGap[0] && !curGap[1] );
////            for ( ; prvK+1 < k; k++ )
////            {
////                assert( false );
////            }
////            for ( int d : { 0, 1 } ) prv[d] = cur[d];
////            prvK = k;
////        }
////        assert( false );
////    }
//}
//
//void Match::updateCoords( Match* match, vector< pair<int,int>>& lcoords, vector< pair<int,int>>& rcoords )
//{
//    Match* m[2]{ this, match };
//    vector< pair<int,int>> coords[2]{ lcoords, rcoords };
//    int scores[2][2]{0};
//    for ( int d : { 0, 1 } )
//    {
////        int gap[2]{ coords[d][0].first ? m[d]->coords_[coords[d][0].first]-m[d]->coords_[coords[d][0].first-1]-1: 0
////                  , coords[d][0].first ? m[!d]->coords_[coords[d][0].second]-m[d]->coords_[coords[d][0].first-1]-1: 0 };
////        for ( int i = 0; i < coords[d].size(); i++ )
////        {
////            int j = coords[d][i].first, k = coords[d][i].second;
////            int coord[2]{ m[d]->coords_[j], m[!d]->coords_[k] };
////            int x = 0;
////        }
//        int gap = coords[d][0].first ? m[d]->coords_[coords[d][0].first]-m[d]->coords_[coords[d][0].first-1]-1 : 0;
//        int score = gap > 0 ? -4 - ( gap * 4 ) : 0;
//        for ( int i = 0; i < coords[d].size(); i++ )
//        {
//            gap = 0;
//            for ( ; i < coords[d].size() && m[d]->isInsert( coords[d][i].first ); i++ ) gap--;
//            if ( !gap && coords[d][i].first+1 < m[d]->coords_.size() ) gap = m[d]->coords_[coords[d][i].first+1]-m[d]->coords_[coords[d][i].first]-1;
//            if ( gap ) score -= 4 + abs( gap * 4 );
//            if ( i >= coords[d].size() ) break;
//            
//            int j = coords[d][i].first, coord = m[d]->coords_[coords[d][i].first];
//            score += m[d]->read_->seq_[j] == m[d]->tar_->seq_[coord] ? 1 : -1;
//        }
//        scores[d][0] = score;
//        
//        gap = coords[d][0].first ? m[!d]->coords_[coords[d][0].second]-m[d]->coords_[coords[d][0].first-1]-1: 0;
//        score = gap ? -4 - abs( gap * 4 ) : 0;
//        assert( gap >= 0 );
//        for ( int i = 0; i < coords[d].size(); i++ )
//        {
//            gap = i+1 < coords[d].size() ? coords[d][i+1].second - coords[d][i].second - 1 : 0;
//            if ( gap < 0 ) for ( ; ++i < coords[d].size() && coords[d][i+1].second == coords[d][i].second; i++ ) gap--;
//            if ( gap ) score -= 4 + abs( gap * 4 );
//            if ( i >= coords[d].size() ) break;
//            
//            int j = coords[d][i].second, coord = m[!d]->coords_[coords[d][i].second];
//            score += m[d]->read_->seq_[j] == m[d]->tar_->seq_[coord] ? 1 : -1;
//        }
//        scores[d][1] = score;
//        int x = 0;
//    }
//    assert( false );
//}
