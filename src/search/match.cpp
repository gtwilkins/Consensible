/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "match.h"
#include "read.h"
#include <cassert>

Match::Match( Target* tar, MappedRead* read, string t, string r, vector<pair<int,int>> &anchors, int tarCoord )
: tar_( tar ), read_( read ), score_( 0 )
{
    assert( t.size() == r.size() );
    assert( !anchors.empty() );
    anchors_[0] = anchors[0].first;
    anchors_[1] = anchors.back().second;
    for ( int i = 0; i < t.size(); i++ )
    {
        if ( i == anchors[0].first ) anchors_[0] = coords_.size();
        if ( i == anchors.back().second ) anchors_[1] = coords_.size();
        if ( i >= anchors[0].first && i <= anchors.back().second ) score_ += t[i] == r[i] ? 1 : -2;
        if ( r[i] != '-' ) coords_.push_back( tarCoord );
        if ( t[i] != '-' ) tarCoord++;
    }
    assert( coords_.size() == read->seq_.size() );
    read_->matches_.push_back( this );
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