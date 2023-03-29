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

#include "consensus_map.h"
#include "match.h"
#include "read.h"
#include "alignment.h"
#include <cassert>
#include <algorithm>
#include <iostream>

int ConMap::excess( int d )
{
    return d ? node_->coords_.size() - mapped_[1] - 1 : mapped_[0];
}

bool ConMap::match( ConMap* cm, vector<pair<ConMap*, AlignResult>>& hits, int drxn )
{
    assert( drxn ? coord_[1] <= cm->coord_[1] : cm->coord_[0] <= cm->coord_[0] );
    int i = cm->range_[drxn];
    if ( drxn ) while ( i && coord_[1] < cm->node_->coords_[i] ) i--;
    else while ( i+1 < cm->node_->coords_.size() && cm->node_->coords_[i] < coord_[0] ) i++;
    int len[2]{ drxn ? node_->size() - range_[1] - 1 : range_[0], drxn ? cm->node_->size() - i - 1 : i };
    AlignResult result = Alignment( node_->read_->seq_, cm->node_->read_->seq_, drxn ? range_[1]+1 : 0, len[0], drxn ? i+1 : 0, len[1] ).align( 1, 4, 4, 4, drxn, !drxn );
    if ( result.len_ > 1 )
    {
        hits.push_back( make_pair( cm, result ) );
        return true;
    }
    return false;
}

void ConMap::updateCoords( ConMap* donor, AlignResult& result, int dIndex, bool drxn )
{
    assert( result.start_ >= 0 && result.start_+result.len_ <= result.s_[0].size() );
    assert( drxn ? result.start_ == 0 : result.start_+result.len_ == result.s_[0].size() );
    
    int len[2]{ result.getSeqLen( !dIndex ), result.getSeqLen( dIndex ) };
    int index[2]{ drxn ? node_->size()-len[0]-1 : len[0], drxn ? donor->node_->size()-len[1]-1 : len[1] };
    int initRange = index[0];
    vector<int> unmatched;
    for ( int i = ( drxn ? 0 : result.s_[0].size()-1 ); i < result.start_+result.len_ && i >= result.start_; drxn ? i++ : i-- )
    {
        for ( int s : { 0, 1 } ) if ( result.s_[s][i] != '-' ) index[s==dIndex] += ( drxn ? 1 : -1 );
        if ( result.s_[dIndex][i] == '-' ) unmatched.push_back( index[0] );
        else if( result.s_[!dIndex][i] != '-' && ( drxn ? range_[1] < index[0] : index[0] < range_[0] ) )
        {
            node_->coords_[index[0]] = donor->node_->coords_[index[1]];
            if ( donor->range_[0] <= index[1] && index[1] <= donor->range_[1] ) range_[drxn] = drxn ? max( range_[1], index[0] ) : min( range_[0], index[0] );
        }
    }
    mapped_[drxn] = drxn ? max( range_[1], mapped_[1] ) : min( range_[0], mapped_[0] );
    coord_[drxn] = node_->coords_[ range_[drxn] ];
    node_->anchors_[drxn] = range_[drxn];
    for ( int i = 0; i < unmatched.size(); )
    {
        int coord = unmatched[i], len = 1;
        while ( ++i < unmatched.size() && unmatched[i] == coord+len ) len++;
        int limits[2]{ coord ? node_->coords_[coord-1] : node_->coords_[coord+len]-len-1, coord+len < node_->size() ? node_->coords_[coord+len] : node_->coords_[coord-1]+len+1 };
        for ( int j = coord; j < coord+len; j++ )
        {
            if ( node_->coords_[j] <= limits[0] ) node_->coords_[j] = limits[0]+1;
            if ( node_->coords_[j] > limits[1] ) node_->coords_[j] = limits[1];
            
        }
    }
    if ( range_[drxn] != initRange ) result.trimFromEnd( !dIndex, abs( initRange - range_[drxn] ), !drxn );
    assert( result.len_ >= 0 );
}

void ConMap::updateCoords( SnpAlignResult& result, bool drxn )
{
    sort( result.snps_.begin(), result.snps_.end(), []( pair<SNPs*, int>& a, pair<SNPs*, int>& b ){ return a.first->start_ == b.first->start_ ? a.first->len_ < b.first->len_ : a.first->start_ < b.first->start_;} );
    int iSnp = 0, coord[2]{ drxn ? coord_[1]+1 : coord_[0]-result.getSeqLen( 0 )+result.lIgnore[0], drxn ? range_[1]+1 : result.lIgnore[1] };
    for ( int i = result.start_; i < result.start_ + result.len_; i++ )
    {
        assert( node_->coords_[coord[1]] == coord[0] );
        node_->coords_[coord[1]] = coord[0];
        mapped_[drxn] = drxn ? max( mapped_[1], coord[1] ) : min( mapped_[0], coord[1] );
        if ( result.s_[0][i] == result.s_[1][i] )
        {
            coord_[drxn] = drxn ? max( coord_[1], coord[0] ) : min( coord_[0], coord[0] );
            range_[drxn] = drxn ? max( range_[1], coord[1] ) : min( range_[0], coord[1] );
        }
        for ( ; iSnp < result.snps_.size() && result.snps_[iSnp].first->start_+result.snps_[iSnp].first->len_ <= coord[0]+(result.s_[0][i] != '-'?1:0); iSnp++ )
        {
            SNPs* snps = result.snps_[iSnp].first;
            SNP* snp = &snps->snps_[result.snps_[iSnp].second];
            assert( coord[0] == snps->start_ );
            if ( !snps->len_ ) assert( result.s_[0][i] == '-' && result.s_[1][i] != '-' );
            if ( snp->seq_.empty() ) assert( result.s_[0][i] != '-' && result.s_[1][i] == '-' );
            if ( result.s_[1][i] != '-' ) assert( string( 1, result.s_[1][i] ) == snp->seq_ );
            snps->addSnp( result.s_[1][i] == '-' ? "" : string( 1, result.s_[1][i] ), node_, coord[1] );
        }
        for ( int s : { 0, 1 } ) if ( result.s_[s][i] != '-' ) coord[s]++;
    }
    node_->anchors_[drxn] = drxn ? max( range_[1], node_->anchors_[1] ) : min( range_[0], node_->anchors_[0] );
}
