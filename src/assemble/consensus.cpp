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

#include "consensus.h"
#include "target.h"
#include "alignment.h"
#include "bubble_resolution.h"
#include <cassert>
#include <algorithm>
#include <iostream>

Consensus::Consensus( vector<Match*>& matches, Target* tar )
: tar_( tar ), kmers_( NULL ), template_( tar->seq_ )
{
    coord_[0] = matches[0]->coords_[matches[0]->anchors_[0]];
    coord_[1] = matches[0]->coords_[matches[0]->anchors_[1]];
    for ( Match* m : matches ) addMatch( m );
    sort( maps_.begin(), maps_.end(), []( ConMap* a, ConMap* b ){ return a->coord_[0] == b->coord_[0] ? a->coord_[1] > b->coord_[1] : a->coord_[0] < b->coord_[0]; } );
}

Consensus::~Consensus()
{
    if ( kmers_ ) delete kmers_;
}

void Consensus::addBranch( ConMap* cm, vector<pair<ConMap*, AlignResult>>& hits, bool drxn )
{
    sort( hits.begin(), hits.end(), []( pair<ConMap*, AlignResult>& a, pair<ConMap*, AlignResult>& b){ return a.second.len_ == b.second.len_ ? a.second.score_ > b.second.score_ : a.second.len_ > b.second.len_; });
    
    int bestI = 0, bestLen = -1, bestScore = 0;
    for ( int i = 0; i < hits.size(); i++ )
    {
        int len = hits[i].second.getSeqLen( 1 );
        int anchor = drxn ? hits[i].first->node_->size() - len - 1 : len;
        int anchorable = drxn ? hits[i].first->range_[1]-anchor : anchor-hits[i].first->range_[0];
        int anchorLen = max( anchorable, len-( drxn ? hits[i].second.rIgnore[1] : hits[i].second.lIgnore[1] ) );
        if ( anchorLen == bestLen ? hits[i].second.score_ > bestScore : anchorLen > bestLen )
        {
            bestI = i;
            bestLen = anchorLen;
            bestScore = hits[i].second.score_;
        }
    }
    if ( bestLen > 1 )
    {
        cm->updateCoords( hits[bestI].first, hits[bestI].second, 1, drxn );
        for ( int i = 0; i < hits.size(); i++ ) if ( i != bestI )
        {
            hits[i].first->updateCoords( cm, hits[i].second, 0, drxn );
        }
    }
    
    for ( int i = 0; i < hits.size(); i++ ) if ( !hits[i].second.len_ ) hits.erase( hits.begin() + i-- );
    if ( hits.empty() ) return;
    
    sort( hits.begin(), hits.end(), []( pair<ConMap*, AlignResult>& a, pair<ConMap*, AlignResult>& b){ return a.second.len_ == b.second.len_ ? a.second.score_ > b.second.score_ : a.second.len_ > b.second.len_; });
    
    Bubble* bub = new Bubble( hits[0].second, cm, hits[0].first, drxn );
    for ( int i = 1; i < hits.size(); i++ ) bub->addMatch( hits[i].second, cm, hits[i].first, drxn );
    branch_[drxn].push_back( bub );
}

bool Consensus::addBubble( Match* l, Match* r, AlignResult& result )
{
    ConMap* maps[2]{ NULL, NULL };
    for ( ConMap* cm : maps_ ) if ( cm->node_ == l ) maps[0] = cm;
    for ( ConMap* cm : maps_ ) if ( cm->node_ == r ) maps[1] = cm;
    maps[0]->mapped_[1] = max( maps[0]->mapped_[1], (int)l->coords_.size()-1-result.rIgnore[0] );
    maps[1]->mapped_[0] = max( maps[1]->mapped_[0], result.lIgnore[1] );
    
    int counts[2]{0}, pSeq = 0;
    int excess[2]{ (int)l->coords_.size()-1-l->anchors_[1]-result.rIgnore[0], -result.lIgnore[1] };
//    int bridgeEnds[2]{ (int)result.s_[0].size(), (int)result.s_[1].size() };
    int bridgedCoord[2]{ (int)l->coords_.size(), (int)r->coords_.size() };
//    int remapAnchors[2]{ max( l->anchors_[1], result.lIgnore[0]-1 ), min( r->anchors_[0], (int)r->coords_.size()-result.rIgnore[1] ) };
    int anchors[2]{ l->anchors_[1], r->anchors_[0] };
    int tarCoords[2]{ l->getAnchor( 1 ), r->getAnchor( 0 ) };
    int gap = r->getAnchor( 0 ) - l->getAnchor( 1 ) - 1;
    string bridge = "", tar = gap > 0 ? template_.substr( l->getAnchor( 1 )+1, gap ) : "";
//    vector<int> remap[2];
//    vector<int> bridgeCoords[2];
    vector< pair<int,int> > remaps[2][2], bridgeCoords[2];
    for ( int i = 0; i < result.s_[0].size(); i++ )
    {
        if ( !pSeq && ( ( result.s_[0] == result.s_[1] && excess[0] <= excess[1] ) || excess[0] <= 0 ) ) pSeq = 1;
        if ( anchors[0] < counts[0] && counts[1] < anchors[1] )
        {
            for ( int d : { 0, 1 } ) if ( result.s_[d][i] != '-' && result.s_[pSeq][i] != '-' ) bridgeCoords[d].push_back( make_pair( counts[d], bridge.size() ) );
            if ( result.s_[pSeq][i] != '-' ) bridge += result.s_[pSeq][i];
//            for ( int d : { 0, 1 } ) bridgeCoords[d].push_back( result.s_[d][i] != '-' ? counts[d] : -1 );
        }
        
        bool anchored[2]{ l->anchors_[0] <= counts[0] && counts[0] <= l->anchors_[1], r->anchors_[0] <= counts[1] && counts[1] <= r->anchors_[1] };
        if ( result.start_ <= i && i < result.start_+result.len_ && l->anchors_[0] <= counts[0] && counts[1] <= r->anchors_[1] ) for ( int d : { 0, 1 } )
        {
            if ( !anchored[d] && result.s_[d][i] != '-' ) remaps[d][anchored[!d]].push_back( make_pair( counts[d], result.s_[!d][i] != '-' ? counts[!d] : -1 ) );
        }
//        if ( remapAnchors[0] < counts[0] && result.s_[0][i] != '-' && i < result.start_+result.len_ ) remap[0].push_back( counts[1] );
//        if ( counts[1] < remapAnchors[1] && result.s_[1][i] != '-' && i >= result.start_ ) remap[1].push_back( counts[0] );
        
        if ( anchors[0] < counts[0] && result.s_[0][i] != '-' ) excess[0]--;
        if ( counts[1] < anchors[1] && result.s_[1][i] != '-' ) excess[1]++;
//        if ( counts[0] == anchors[0] && result.s_[0][i] != '-' ) bridgeEnds[0] = i;
//        if ( counts[1] == anchors[1] && result.s_[1][i] != '-' ) bridgeEnds[1] = i;
        if ( anchors[1] <= counts[1] && counts[0] < bridgedCoord[0] && result.s_[0][i] != '-' ) bridgedCoord[0] = counts[0];
        if ( counts[0] <= anchors[0] && result.s_[1][i] != '-' ) bridgedCoord[1] = counts[1];
        for ( int d : { 0, 1 } ) if ( result.s_[d][i] != '-' ) counts[d]++;
    }
    
    assert( gap >= 0 || gap == bridgedCoord[1] - bridgedCoord[0] - 1 );
    
//    l->anchorCoords( r, remaps[0][1] );
//    r->anchorCoords( l, remaps[1][1] );
    
//    l->updateCoords( r, remap[0], remapAnchors[0]+1 );
//    r->updateCoords( l, remap[1], result.lIgnore[1] );
    
    Bubble* bub = new Bubble();
    bub->template_ = bridge;
    bub->start_ = tarCoords[0]+1;
    bub->len_ = gap;
    vector<int> bubbleCoords;
    if ( !bridge.empty() )
    {
        AlignResult bridgeResult = Alignment( bridge, tar, 0, bridge.size(), 0, tar.size() ).align( true, true );
        int tarCoord = bub->start_;
        for ( int i = 0; i < bridgeResult.s_[0].size(); i++ )
        {
            if ( bridgeResult.s_[0][i] != '-' ) bubbleCoords.push_back( tarCoord );
            if ( bridgeResult.s_[1][i] != '-' ) tarCoord++;
        }
    }
    
    Match* m[2]{ l, r };
    for ( int d : { 0, 1 } ) for ( pair<int, int> coord : bridgeCoords[d] )
    {
        m[d]->coords_[coord.first] = coord.second < bubbleCoords.size() ? bubbleCoords[coord.second] : tarCoords[1];
    }
    for ( int d : { 0, 1 } ) m[d]->anchorCoords( m[!d], remaps[d][1] );
    
    bool bridged[2]{ bridgedCoord[0] < l->coords_.size(), bridgedCoord[1] < r->coords_.size() };
    if ( bridged[0] ) addMatch( l, bridgedCoord[0], l->coords_.size()-1-result.rIgnore[0], true );
    if ( bridged[1] ) addMatch( r, result.lIgnore[1], bridgedCoord[1], true );
    if ( bridged[0] && bridged[1] && gap < 2 && gap >= 0 && bridgedCoord[0]-anchors[0] < 2 && anchors[1]-bridgedCoord[1] < 2  )
    {
        addMismatch( l, anchors[0], bridgedCoord[0] );
        addMismatch( r, bridgedCoord[1], anchors[1] );
        return true;
    }
    
    bub->maps_.push_back( Bubble::BubbleMap( l, bridgeCoords[0], anchors[0], bridged[0] ? bridgedCoord[0] : -1 ) );
    bub->maps_.push_back( Bubble::BubbleMap( r, bridgeCoords[1], bridged[1] ? bridgedCoord[1] : -1, anchors[1] ) );
    bubble_.push_back( bub );
    
//    if ( !bridge.empty() )
//    {
//        Bubble::BubbleMap maps[2]{ Bubble::BubbleMap( l, anchors[0], bridged[0] ? bridgedCoord[0] : -1 ), Bubble::BubbleMap( r, anchors[1], bridged[1] ? bridgedCoord[1] : -1 )  };
//        AlignResult bridgeResult = Alignment( bridge, tar, 0, bridge.size(), 0, tar.size() ).align( true, true );
//        counts[0] = 0;
//        counts[1] = tarCoords[0]+1;
//        Match* m[2]{ l, r };
//        for ( int i = 0; i < bridgeResult.s_[0].size(); i++ )
//        {
//            if ( bridgeResult.s_[0][i] != '-' )
//            {
//    //            for ( int d : { 0, 1 } ) if ( bridgeCoords[d][counts[0]] >= 0 ) assert( m[d]->coords_[bridgeCoords[d][counts[0]]] == counts[1] );
//                for ( int d : { 0, 1 } ) if ( bridgeCoords[d][counts[0]] >= 0 )
//                {
//                    m[d]->coords_[bridgeCoords[d][counts[0]]] = counts[1];
//                    
//                }
//                bub.template_ += bridgeResult.s_[0][i];
//            }
//            for ( int d : { 0, 1 } ) if ( bridgeResult.s_[d][i] != '-' ) counts[d]++;
//        }
//    }
    int x = 0;
    return true;
    
    
//    int gap = r->getAnchor( 0 ) - l->getAnchor( 1 ) - 1,  coords[2]{0}, counts[2]{ 0 };
//    int good[2][2]{ { l->anchors_[1], l->coords_.size() }, { -1, r->anchors_[0] } };
//    for ( int i = 0; i < result.s_[0].size(); i++ )
//    {
//        if ( i >= result.start_ && i < result.start_+result.len_ )
//        {
//            if ( counts[0] <= l->anchors_[1] && counts[1] < r->anchors_[0] )
//            {
//                assert( l->coords_[counts[0]] == r->coords_[counts[1]]);
//            }
//            if ( l->anchors_[1] < counts[0] && r->anchors_[0] <= counts[1] )
//            {
//                assert( l->coords_[counts[0]] == r->coords_[counts[1]]);
//            }
//        }
////        if ( counts[0] == stops[0] && result.s_[0][i] != '-' ) coords[0] = i;
////        if ( counts[1] == stops[1] && result.s_[1][i] != '-' ) coords[1] = i;
//        if ( counts[0] == l->anchors_[1] && result.s_[0][i] != '-' ) good[1][0] = result.s_[1][i] != '-' ? counts[1] : counts[1]-1;
//        if ( counts[1] == r->anchors_[0] && result.s_[1][i] != '-' ) good[0][1] = result.s_[0][i] != '-' ? counts[0] : counts[0]+1;
//        for ( int d : { 0, 1 } ) if ( result.s_[d][i] != '-' ) counts[d]++;
//    }
////    addMismatch( l, good[0][0], good[0][1] );
////    addMismatch( r, good[1][0], good[1][1] );
//    
//    string seq[2];
//    for ( int i = coords[0]+1; i < coords[1]; i++ )
//    {
//        if ( result.s_[0][i] != '-' ) seq[0] += result.s_[0][i];
//        if ( result.s_[1][i] != '-' ) seq[1] += result.s_[1][i];
//    }
//    if ( max( seq[0].size(), seq[1].size() ) > 1 || gap > 1 )
//    {
//        assert( false );
//        return true;
//    }
//    
//    SNPs snps;
//    snps.start_ = l->getAnchor( 1 )+1;
//    snps.len_ = gap;
//    snps.snps_.push_back( SNPs::SNP() );
//    snps.snps_.back().seq_ = seq[0];
//    snps.snps_.back().matches_.push_back( make_pair( l, stops[0]+1 ) );
//    if ( seq[0] != seq[1] ) snps.snps_.push_back( SNPs::SNP() );
//    snps.snps_.back().seq_ = seq[1];
//    snps.snps_.back().matches_.push_back( make_pair( r, stops[1]-seq[1].size() ) );
//    snps_.push_back( snps );
//    assert( false );
    return false;
}

void Consensus::addMatch( Match* m )
{
    ConMap* cm = new ConMap;
    cm->node_ = m;
    cm->coord_[0] = m->coords_[m->anchors_[0]];
    cm->coord_[1] = m->coords_[m->anchors_[1]];
    cm->range_[0] = cm->mapped_[0] = m->anchors_[0];
    cm->range_[1] = cm->mapped_[1] = m->anchors_[1];
    coord_[0] = min( coord_[0], cm->coord_[0] );
    coord_[1] = max( coord_[1], cm->coord_[1] );
    maps_.push_back( cm );
    
    addMatch( m, m->anchors_[0], m->anchors_[1], false );
    
//    int prv = m->anchors_[1], cur = m->anchors_[1];
//    int minLen = 1;
//    for ( int i = m->anchors_[0]; i <= m->anchors_[1]; i++ )
//    {
//        if ( ( i == m->anchors_[1] || i == cur + minLen ) && prv < cur ) addMismatch( m, prv, cur );
//        if ( i >= m->anchors_[1] ) break;
//        bool ins = i+1 < m->coords_.size() && m->coords_[i] == m->coords_[i+1];
//        bool match = !ins && m->read_->seq_[i] == m->tar_->seq_[m->coords_[i]];
//        bool del = i && ( m->coords_[i-1]+1 < m->coords_[i] );
//        
//        if ( ( !match || del ) && cur + minLen <= i ) prv = i-1;
//        if ( !match ) cur = m->anchors_[1];
//        if ( match && ( i < cur || del ) ) cur = i;
//    }
}

void Consensus::addMatch( Match* m, int blockStart, int blockLast, bool preexisting )
{
    if ( preexisting )
    {
        if ( blockStart < m->anchors_[0] ) m->anchors_[0] = blockStart;
        if ( blockLast > m->anchors_[1] ) m->anchors_[1] = blockLast;
        for ( ConMap* cm : maps_ ) if ( cm->node_ == m )
        {
            cm->coord_[0] = min( cm->coord_[0], m->coords_[m->anchors_[0]] );
            cm->coord_[1] = max( cm->coord_[1], m->coords_[m->anchors_[1]] );
            cm->range_[0] = min( cm->range_[0], m->anchors_[0] );
            cm->range_[1] = max( cm->range_[1], m->anchors_[1] );
            cm->mapped_[0] = min( cm->mapped_[0], cm->range_[0] );
            cm->mapped_[1] = max( cm->mapped_[1], cm->range_[1] );
            coord_[0] = min( coord_[0], cm->coord_[0] );
            coord_[1] = max( coord_[1], cm->coord_[1] );
            break;
        }
    }
    int prv = blockLast, cur = blockLast;
    int minLen = 1;
    for ( int i = blockStart; i <= blockLast; i++ )
    {
        if ( i == min( blockLast, cur + minLen ) && prv < cur )
        {
            addMismatch( m, prv, cur );
        }
        if ( i >= blockLast ) break;
        bool ins = i+1 < m->coords_.size() && m->coords_[i] == m->coords_[i+1];
        bool match = !ins && m->read_->seq_[i] == m->tar_->seq_[m->coords_[i]];
        bool del = i && ( m->coords_[i-1]+1 < m->coords_[i] );
        
        if ( ( !match || del ) && cur + minLen <= i ) prv = i-1;
        if ( !match ) cur = blockLast;
        if ( match && ( i < cur || del ) ) cur = i;
    }
}

void Consensus::addMismatch( Match* match, int lGood, int rGood )
{
    int matchLen = rGood - lGood - 1;
    int tarLen = lGood < 0 || rGood == match->anchors_[1] ? matchLen : match->coords_[rGood] - match->coords_[lGood] - 1;
    if ( matchLen < 2 && tarLen < 2 )
    {
        int start = match->coords_[lGood+1];
        int i = 0;
        while ( i < snps_.size() && ( snps_[i]->start_ == start ? snps_[i]->len_ < tarLen : snps_[i]->start_ < start ) ) i++;
        if ( i == snps_.size() || snps_[i]->start_ != start || snps_[i]->len_ != tarLen  )
        {
            snps_.insert( snps_.begin() + i, new SNPs() );
            snps_[i]->start_ = start;
            snps_[i]->len_ = tarLen;
        }
        SNP snp;
        snp.seq_ = matchLen ? match->read_->seq_.substr( lGood+1, 1 ) : "";
        int j = 0;
        while ( j < snps_[i]->snps_.size() && snps_[i]->snps_[j].seq_ != snp.seq_ ) j++;
        if ( j == snps_[i]->snps_.size() ) snps_[i]->snps_.insert( snps_[i]->snps_.begin() + j, snp );
        snps_[i]->snps_[j].matches_.push_back( make_pair( match, lGood+1 ) );
    }
    else
    {
        Bubble* bubble = new Bubble();
        bubble->start_ = match->coords_[lGood+1];
        bubble->len_ = tarLen;
        bubble->template_ = match->read_->seq_.substr( lGood+1, matchLen );
        bubble_.push_back( bubble );
    }
}

void Consensus::addSNP( Match*, int start, int tarLen, string seq, int coord )
{
    int i = 0;
    while ( i < snps_.size() && ( snps_[i]->start_ == start ? snps_[i]->len_ < tarLen : snps_[i]->start_ < start ) ) i++;
    if ( i == snps_.size() || snps_[i]->start_ != start || snps_[i]->len_ != tarLen  )
    {
        snps_.insert( snps_.begin() + i, new SNPs() );
        snps_[i]->start_ = start;
        snps_[i]->len_ = tarLen;
    }
}

bool Consensus::bridge( Consensus* lhs, Consensus* rhs )
{
    if ( !lhs->kmers_ ) lhs->mapKmers();
    if ( !rhs->kmers_ ) rhs->mapKmers();
    
    ConsensusKmers* kmers[2] = { lhs->kmers_, rhs->kmers_ };
    unordered_map<Match*, unordered_map<Match*,vector< pair<int, int> > > > hits;
    for ( int d : { 0, 1 } )
    {
        Kmer* hit = NULL;
        for ( const pair<uint16_t, Kmer>& k : kmers[!d]->ends_[d] ) if ( hit = kmers[d]->getKmer( k.first ) )
        {
            for ( pair<Match*, int> a : k.second.coords_ ) for ( pair<Match*, int> b : hit->coords_ )
            {
                if ( d && b.second < b.first->anchors_[0] ) continue;
                Match* matches[2]{ d ? a.first : b.first, d ? b.first : a.first };
                pair<int,int> coords{ d ? a.second : b.second, d ? b.second : a.second };
                auto ins1 = hits.insert( make_pair( matches[0], unordered_map<Match*, vector< pair<int, int> > >{} ) );
                auto ins2 = ins1.first->second.insert( make_pair( matches[1], vector< pair<int, int> >{} ) );
                ins2.first->second.push_back( coords );
            }
        }
    }
    
    typedef pair<pair<Match*, Match*>, vector< pair<int, int> > > MatchHits;
    vector< MatchHits > hitLists;
    for ( const pair<Match*, unordered_map<Match*,vector< pair<int, int> > > >& l : hits ) for ( const pair< Match*,vector< pair<int, int> > >& r : l.second )
    {
        vector< pair<int, int> > coords;
        for ( pair<int, int> coord : r.second )
        {
            int lSpare = l.first->anchors_[1] < coord.first ? l.first->anchors_[1]+1-coord.first : max( 0, l.first->anchors_[1]-coord.first-7 );
            int rSpare = r.first->anchors_[0] < coord.second ? coord.second-r.first->anchors_[0] : min( 0, coord.second+8-r.first->anchors_[0] );
            int overlap = max( 0, l.first->getAnchor( 1 ) - r.first->getAnchor( 0 ) );
            if ( lSpare + rSpare - overlap < 10 ) coords.push_back( coord );
        }
        if ( !coords.empty() )
        {
            sort( coords.begin(), coords.end(), []( pair<int, int>& a, pair<int, int>& b ){ return a.first == b.first ? a.second < b.second : a.first < b.first; } );
            hitLists.push_back( make_pair( make_pair( l.first, r.first ), coords ) );
        }
    }
    
    sort( hitLists.begin(), hitLists.end(), []( MatchHits& a, MatchHits& b ){ return a.second.size() > b.second.size(); } );
    for ( MatchHits& hit : hitLists )
    {
        vector< pair< pair<int, int>, int> > blocks;
        for ( int i = 0; i < hit.second.size(); i++ )
        {
            pair<int, int> coords = hit.second[i];
            int len = 8;
            for ( ; i+1 < hit.second.size(); i++ )
            {
                if ( hit.second[i].first+1 == hit.second[i+1].first && hit.second[i].second+1 == hit.second[i+1].second ) len++;
                else if ( hit.second[i].first+9 == hit.second[i+1].first && hit.second[i].second+9 == hit.second[i+1].second ) len += 9;
                else break;
            }
            blocks.push_back( make_pair( coords, len ) );
        }
        sort( blocks.begin(), blocks.end(), []( pair< pair<int, int>, int>& a, pair< pair<int, int>, int>& b ){ return a.second > b.second; } );
        for ( int i = 0; i < blocks.size(); i++ )
        {
            AlignResult result = Alignment::alignBySeed( hit.first.first->read_->seq_, hit.first.second->read_->seq_, blocks[i].first.first, blocks[i].first.second, blocks[i].second );
            int unalignedCons[2]{ (int)hit.first.first->coords_.size() - hit.first.first->anchors_[1] - 1, hit.first.second->anchors_[0] };
            int excessCons[2]{ lhs->coord_[1] - hit.first.first->getAnchor( 1 ), hit.first.second->getAnchor( 0 )-rhs->coord_[0] };
            int overlap = max( 0, hit.first.first->getAnchor( 1 ) + unalignedCons[0] - hit.first.second->getAnchor( 0 ) + unalignedCons[1] );
            int sparePenalty = max( 0, result.rIgnore[0] - unalignedCons[0] ) + max( 0, result.lIgnore[1] - unalignedCons[1] );
            if ( Match::isBridge( hit.first.first, hit.first.second, result, excessCons[0]+excessCons[1] ) )
            {
                return lhs->merge( rhs, result, hit.first.first, hit.first.second );
            }
        }
        int x = 0;
    }
    int x = 0;
    return false;
}

bool Consensus::foldEnd( ConMap* cm, vector<Bubble*>* branches, bool force, bool d )
{
//    int tLen = min( cm->unmapped( d )+10, d ? (int)template_.size()-cm->coord_[1]-1 : cm->coord_[0] );
//    int tStart = d ? cm->coord_[1]+1 : cm->coord_[0]-tLen;
//    
//    vector<SNPs*> snps;
//    vector<Bubble*> bubbles;
    
    
    
//    if ( branches )
//    {
//        for ( vector<Bubble*>* bs : { branches, &branch_[d] } ) for ( Bubble* b : *bs ) if ( tStart <= b->start_ && b->start_ <= tStart+tLen ) bubbles.push_back( b );
//        if ( bubbles.empty() ) return false;
//    }
//    for ( SNPs* s : snps_ ) if ( tStart <= s->start_ && s->start_ + s->len_ <= tStart + tLen ) snps.push_back( s );
//    for ( Bubble* b : bubble_ ) if ( ( tStart <= b->start_ && b->start_ <= tStart + tLen ) || ( tStart <= b->start_ + b->len_ && b->start_ + b->len_ <= tStart+tLen ) ) bubbles.push_back( b );
//    if ( !force && snps.empty() && bubbles.empty() ) return false;
//    
//    SnpAlignment align( template_, cm->node_->read_->seq_, tStart, tLen, d ? cm->range_[1]+1 : 0, cm->unmapped( d ), snps, bubbles );
//    SnpAlignResult result = align.align( d, !d );
    
    int tLen = min( cm->unmapped( d )+10, d ? (int)template_.size()-cm->coord_[1]-1 : cm->coord_[0] );
    int limits[2]{ d ? cm->coord_[1]+1 : cm->coord_[0]-tLen, d ? cm->coord_[1]+1+tLen : cm->coord_[0] };
    int read[2]{ d ? cm->range_[1]+1 : 0, d ? cm->node_->size() : cm->range_[0] };
    
    vector<SNPs*> snps;
    vector<Bubble*> bubbles;
    for ( pair<bool,bool> firstSec = make_pair( true, false ); firstSec.first || firstSec.second; firstSec.first = false )
    {
        bubbles.clear();
        if ( firstSec.second ) cm->setFinalFoldCoords( limits, read, d );
        if ( branches ) for ( vector<Bubble*>* bs : { branches, &branch_[d] } ) for ( Bubble* b : *bs ) if ( b->isFoldTarget( cm, limits, firstSec.second, d ) ) bubbles.push_back( b );
        if ( branches && bubbles.empty() ) return false;
        for ( Bubble* b : bubble_ ) if ( b->isFoldTarget( cm, limits, firstSec.second, d ) ) bubbles.push_back( b );
        firstSec.second = ( firstSec.first && ( d ? limits[0] <= cm->coord_[1] : cm->coord_[0] < limits[1] ) );
    }
    for ( SNPs* s : snps_ ) if ( s->isFoldTarget( limits ) ) snps.push_back( s );
    if ( !force && snps.empty() && bubbles.empty() ) return false;
    
    SnpAlignment align( template_, cm->node_->read_->seq_, limits[0], limits[1]-limits[0], read[0], read[1]-read[0], snps, bubbles );
    SnpAlignResult result = align.align( d, !d );
    
    if ( !result.len_ ) return false;
    cm->unsetMappedEnd( read[!d], &bubbles, &snps, d );
    updateMapped( cm, result, d );
    return true;
}

void Consensus::foldEnds()
{
    for ( ConMap* cm : maps_ ) assert( cm->range_[1] < cm->node_->size() );
    for ( ConMap* cm : maps_ ) assert( cm->mapped_[1] < cm->node_->size() );
    sort( maps_.begin(), maps_.end(), []( ConMap* a, ConMap* b ){ return a->coord_[0] == b->coord_[0] ? a->coord_[1] > b->coord_[1] : a->coord_[0] < b->coord_[0]; } );
    sort( snps_.begin(), snps_.end(), []( SNPs* a, SNPs* b ){ return a->start_ == b->start_ ? a->start_+a->len_ > b->start_+b->len_ : a->start_ < b->start_; } );
    
    vector<ConMap*> ends[2];
    for ( ConMap* cm : maps_ ) for ( int d : { 0, 1 } ) if ( cm->unmapped( d ) ) ends[d].push_back( cm );
    bool force = true;
    while ( !ends[0].empty() || !ends[1].empty() )
    {
        unordered_map<SNPs*,int> oldSnps;
        unordered_map<Bubble*,int> oldBubs;
        for ( SNPs* s : snps_ ) oldSnps.insert( make_pair( s, s->snps_.size() ) );
        for ( Bubble* b : bubble_ ) oldBubs.insert( make_pair( b, b->maps_.size() ) );
        
        for ( int d : { 0, 1 } ) for ( ConMap* cm : ends[d] ) foldEnd( cm, NULL, force, d );
        
        vector<pair<int,int>> remapped;
        SNPs::setRemapped( remapped, oldSnps, snps_ );
        Bubble::setRemapped( remapped, oldBubs, bubble_ );
        for ( int d : { 0, 1 } ) ends[d] = ConMap::getFoldableEnds( maps_, remapped, d );
        force = false;
    }
    
    for ( int d : { 0, 1 } ) for ( ConMap* cm : maps_ ) if ( cm->unmapped( d ) && cm->range_[d] == cm->mapped_[d] ) ends[d].push_back( cm );
    sort( ends[0].begin(), ends[0].end(), [] ( ConMap* a, ConMap* b ){ return a->coord_[0] == b->coord_[0] ? a->uncoorded( 0 ) > b->uncoorded( 0 ) : a->coord_[0] < b->coord_[0]; } );
    sort( ends[1].begin(), ends[1].end(), [] ( ConMap* a, ConMap* b ){ return a->coord_[1] == b->coord_[1] ? a->uncoorded( 1 ) > b->uncoorded( 1 ) : a->coord_[1] > b->coord_[1]; } );
    
    for ( int d : { 0, 1 } )
    {
        vector<Bubble*> branches;
        for ( ConMap* cm : ends[d] )
        {
            foldEnd( cm, &branches, false, d );
            if ( cm->unmapped( d ) && cm->range_[d] == cm->mapped_[d] ) branches.push_back( new Bubble( cm, d ) );
        }
        for ( Bubble* b : branches )
        {
            if ( b->maps_.size() > 1 ) branch_[d].push_back( b );
            else b->abort();
        }
        
    }
    int x = 0;
    
//    unordered_map<ConMap*, AlignResult> aligns[2];
//    if ( !snps_.empty() ) for ( int d : { 0, 1 } )
//    {
//        sort( snps_.begin(), snps_.end(), [&]( SNPs* a, SNPs* b ){ return d ? a->start_ < b->start_ : a->start_+a->len_ < b->start_+b->len_; } );
//        for ( ConMap* cm : maps_ ) if ( cm->unmapped( d ) )
//        {
//            int tLen = min( cm->unmapped( d )+10, d ? (int)template_.size()-cm->coord_[1]-1 : cm->coord_[0] );
//            int tStart = d ? cm->coord_[1]+1 : cm->coord_[0]-tLen;
//            vector<SNPs*> snps;
//            vector<Bubble*> bubbles;
//            for ( SNPs* s : snps_ ) if ( tStart <= s->start_ && s->start_ + s->len_ <= tStart + tLen ) snps.push_back( s );
//            for ( Bubble* b : bubble_ ) if ( ( tStart <= b->start_ && b->start_ <= tStart + tLen ) || ( tStart <= b->start_ + b->len_ && b->start_ + b->len_ <= tStart+tLen ) ) bubbles.push_back( b );
////            for ( Bubble* b : branch_[0] ) if ( tStart <= b->start_ || b->start_ + b->len_ <= tStart+tLen ) bubbles.push_back( b );
//            if ( snps.empty() && bubbles.empty() ) continue;
//            SnpAlignment align( template_, cm->node_->read_->seq_, tStart, tLen, d ? cm->range_[1]+1 : 0, cm->unmapped( d ), snps, bubbles );
//            SnpAlignResult result = align.align( d, !d );
//            if ( !result.len_ ) continue;
//            updateMapped( cm, result, d );
////            cm->updateCoords( result, d );
//            
//            
////            int limits[2]{ cm->coord_[d] + ( d ? 1 : -cm->excess( 0 )-5 ), cm->coord_[d] + ( d ? cm->excess( 1 )+5 : 0 ) };
////            int snps[2]{ (int)snps_.size(), 0 };
////            for ( int i = 0; i < snps_.size() && snps_[i].start_ < limits[1]; i++ ) if ( limits[0] <= snps_[i].start_ && snps_[i].start_+snps_[i].len_ <= limits[1] )
////            {
////                snps[0] = min( snps[0], i );
////                snps[1] = max( snps[1], i );
////            }
////            if ( snps[1] < snps[0] ) continue;
////            AlignResult* result = NULL;
////
////            int tLen = min( cm->excess( d )+10, d ? (int)template_.size()-cm->coord_[1]-1 : cm->coord_[0] );
////            int tStart = d ? cm->coord_[1]+1 : cm->coord_[0]-tLen;
////
////            auto it = aligns[d].find( cm );
////            if ( it != aligns[d].end() ) result = &it->second;
////            else
////            {
////                AlignResult ar = Alignment( template_, cm->node_->read_->seq_, tStart, tLen, d ? cm->range_[1]+1 : 0, cm->excess( d ) ).align( d, !d );
////                auto ins = aligns[d].insert( make_pair( cm, ar ) );
////                result = &ins.first->second;
////            }
////            
////            int coord = tStart + ( d ? 0 : tLen-1 ), excess = cm->excess( d ), score = 0;
////            for ( int i = d ? 0 : result->s_[0].size()-1; i >= 0 && i < result->s_[0].size(); d ? i++ : i-- )
////            {
////                int gap[2]{ 0 };
////                for ( int s : { 0, 1 } ) for ( ; result->s_[s][i] == '-' && i >= 0 && i < result->s_[0].size(); d ? i++ : i-- ) gap[s]++;
////                if ( i < 0 || i >= result->s_[0].size() ) break;
////                assert( !gap[0] && !gap[1] );
////
////                bool match = result->s_[0][i] == result->s_[1][i];
////                score += match ? 1 : -4;
////            }
//
////            if ( !d ) for ( ; i < result->s_[0].size() && result->s_[1] == '-' && result->s_[0] != '-'; i++ ) coord++;
////            for ( ; i < result->s_[0].size() && ( !d || excess > 0 ); i++ ) if ( result->s_[0][i] != result->s_[1][i] )
////            {
////
////            }
////            assert( false );
//        }
//    }
}

void Consensus::mapKmers()
{
    if ( kmers_ ) delete kmers_;
    kmers_ = new ConsensusKmers();
    
    for ( ConMap* cm : maps_ )
    {
        kmers_->addMatch( cm->node_, cm->range_ );
    }
}

bool Consensus::merge( Consensus* rhs, AlignResult& result, Match* l, Match* r )
{
    assert( tar_ && tar_ == rhs->tar_ );
    coord_[1] = max( coord_[1], rhs->coord_[1] );
//    template_ += tar_->seq_.substr( coord_[1], rhs->coord_[1]- coord_[1] );
    maps_.insert( maps_.end(), rhs->maps_.begin(), rhs->maps_.end() );
    snps_.insert( snps_.end(), rhs->snps_.begin(), rhs->snps_.end() );
    bubble_.insert( bubble_.end(), rhs->bubble_.begin(), rhs->bubble_.end() );
    for ( int d : { 0, 1 } ) branch_[d].insert( branch_[d].end(), rhs->branch_[d].begin(), rhs->branch_[d].end() );
    if ( kmers_ ) delete kmers_;
    kmers_ = NULL;
    rhs->maps_.clear();
    rhs->snps_.clear();
    rhs->bubble_.clear();
    for ( int d : { 0, 1 } ) rhs->branch_[d].clear();
    addBubble( l, r, result );
    return true;
}

void Consensus::resolveBubbles()
{
    for ( Bubble* b : bubble_ )
    {
        BubbleResolution* br = new BubbleResolution( b, template_ );
        assert( false );
    }
}

void Consensus::resolveBranches()
{
//    vector<ConMap> map( maps_.begin(), maps_.end() );
//    for ( int d : { 0, 1 } )
//    {
//        sort( map.begin(), map.end(), [&]( ConMap& a, ConMap& b ){ return d ? a.coord_[0] < b.coord_[0] : a.coord_[1] < b.coord_[1]; } );
//    }
//    for ( int d : { 0, 1 } ) for ( Branch& b : branch_[d] )
//    {
//        
//    }
}

string Consensus::resolve()
{
    foldEnds();
//    setBubbles();
    setBranches();
//    resolveBranches();
    resolveBubbles();
    consensus_ = "";
    int i = 0, j = 0, coord = coord_[0];
    while ( i < bubble_.size() || j < snps_.size() )
    {
        for ( ; i < bubble_.size() && ( j >= snps_.size() || bubble_[i]->start_ < snps_[j]->start_+snps_[j]->len_ ) ; i++ )
        {
            
        }
        
        for ( ; j < snps_.size() && ( i >= bubble_.size() || snps_[j]->start_+snps_[j]->len_ < bubble_[i]->start_ ) ; j++ )
        {
            consensus_ += template_.substr( coord, snps_[j]->start_-coord ) + snps_[j]->resolve( template_.substr( snps_[j]->start_, snps_[j]->len_ ) );
            coord = snps_[j]->start_ + snps_[j]->len_;
        }
    }
    consensus_ += template_.substr( coord, coord_[1]-coord );
//    cout << ">Template_" << coord_[0] << endl;
//    cout << template_.substr( coord_[0], coord_[1]-coord_[0] ) << endl;
//    cout << ">Consensus_" << coord_[0] << endl;
//    cout << consensus_ << endl;
//    int x = 0;
    return consensus_;
}

void Consensus::setBranches()
{
    vector<ConMap*> tangents[2];
    sort( maps_.begin(), maps_.end(), [&]( ConMap* a, ConMap* b ){ return a->coord_[0] == b->coord_[0] ? a->coord_[1] > b->coord_[1] : a->coord_[0] < b->coord_[0]; } );
    for ( int d : { 1, 0 } )
    {
        vector<ConMap*> maps = d ? vector<ConMap*>( maps_.rbegin(), maps_.rend() ) : maps_;
        if ( d ) sort( maps.begin(), maps.end(), []( ConMap* a, ConMap* b ){ return a->coord_[1] == b->coord_[1] ? a->coord_[0] < b->coord_[0] : a->coord_[1] > b->coord_[1]; } );
        sort( snps_.begin(), snps_.end(), [&]( SNPs* a, SNPs* b ){ return d ? a->start_ < b->start_ : a->start_+a->len_ < b->start_+b->len_; } );
//        sort( bubble_.begin(), bubble_.end(), [&]( Bubble& a, Bubble& b ){ return d ? a.start_ < b.start_ : a.start_+a.len_ < b.start_+b.len_; } );
        
        int excess;
        for ( int i = 1; i < maps.size(); i++ ) if ( excess = maps[i]->unmapped( d ) )
        {
            int limits[2]{ d ? maps[i]->coord_[1] : maps[i]->coord_[0]-excess, d ? maps[i]->coord_[1]+excess : maps[i]->coord_[0] };
            
            // Check nearby bubbles
            
            // Check nearby SNPs
            
            // Check nearby unmapped ends
            vector<pair<ConMap*, AlignResult>> hits;
            for ( int j = i; j-- > 0 && abs( maps[j]->coord_[d] - maps[i]->coord_[d] ) < excess; )
            {
                maps[i]->match( maps[j], hits, d );
            }
            if ( !hits.empty() ) addBranch( maps[i], hits, d );
        }
    }
    int x = 0;
}

void Consensus::setBubbles()
{
    if ( bubble_.empty() ) return;
    sort( maps_.begin(), maps_.end(), [&]( ConMap* a, ConMap* b ){ return a->coord_[0] == b->coord_[0] ? a->coord_[1] > b->coord_[1] : a->coord_[0] < b->coord_[0]; } );
    sort( bubble_.begin(), bubble_.end(), []( Bubble* a, Bubble* b ){ return a->start_ == b->start_ ? a->len_ > b->len_ : a->start_ < b->start_; } );
    int i = 0;
    for ( Bubble* b : bubble_ )
    {
        while ( i < maps_.size() && maps_[i]->coord_[1]+maps_[i]->unmapped( 1 ) >= b->start_ ) i++;
        for ( int j = i; j < maps_.size() && maps_[j]->coord_[1] <= b->start_; j++ ) if ( maps_[j]->coord_[1] + maps_[j]->unmapped( 1 ) > b->start_ )
        {
            assert( false );
        }
    }
}

void Consensus::updateMapped( ConMap* cm, SnpAlignResult& result, bool drxn )
{
    int base = drxn ? cm->coord_[1]+1 : cm->coord_[0]-result.getSeqLen( 0 );
    int coord[2]{ base, drxn ? cm->range_[1]+1 : 0 };
    int badCoord[2]{ coord[0], coord[1] }, badLen = 0;
    string badStr = "";
    
    vector<SnpAlignResult::BubbleAlignCoords>::iterator it = result.bubbles_.begin();
    for ( int i = 0; i < result.s_[0].size(); i++ )
    {
        for ( ; it != result.bubbles_.end() && i == it->start_; it++ )
        {
            if ( it->bubble_ ) it->bubble_->addMatch( result, *it, cm, i, coord[1] );
            if ( it->snps_ ) it->snps_->addMatch( result, *it, cm, i, coord[1] );
            coord[0] = base + it->coord_[1];
        }
//        cout << tar_->seq_.substr( coord[0], 10 ) << endl;
//        cout << result.s_[0].substr( i ) << endl;
        if ( i >= result.start_ && i < result.start_+result.len_ )
        {
            cm->updateCoord( coord[0], coord[1] );
            if ( result.s_[0][i] != result.s_[1][i] )
            {
                if ( !badLen && badStr.empty() ) for ( int s : { 0, 1 } ) badCoord[s] = coord[s];
                if ( result.s_[0][i] != '-' ) badLen++;
                if ( result.s_[1][i] != '-' ) badStr += result.s_[1][i];
            }
            bool badRun = badLen || !badStr.empty();
            bool ending = i+1 == result.s_[0].size();
            bool bubbleNext = it != result.bubbles_.end() && i+1 == it->start_;
            if ( badRun && ( ending || bubbleNext || result.s_[0][i+1] == result.s_[1][i+1] ) )
            {
                if ( badLen <= 1 && badStr.size() <= 1 )
                {
                    SNPs::addSnp( snps_, badStr, cm, badCoord[0], badCoord[1], badLen );
                }
                else
                {
                    Bubble::addBubble( bubble_, badStr, cm, badCoord[0], badCoord[1], badLen );
                }
                badStr.clear();
                badLen = 0;
            }
        }
        
        for ( int s : { 0, 1 } ) if ( result.s_[s][i] != '-' ) coord[s]++;
    }
//    vector<int> starts;
//    
//    bool bubbled = false;
//    for ( int i = 0; i < result.s_[0].size(); i++ )
//    {
//        if ( bubbled && i == result.bubbles_[starts.size()-1].start_ ) bubbled = false;
//        if ( !bubbled && ( bubbled == starts.size() < result.bubbles_.size() && i == result.bubbles_[starts.size()].start_ ) ) starts.push_back( coord[1] );
//        
//        if ( !bubbled ) assert( cm->node_->coords_[coord[1]] == coord[0] );
//        if ( !bubbled ) cm->node_->coords_[coord[1]] = coord[0];
//        
//        if ( i >= result.start_ && i < result.start_+result.len_ )
//        {
//            cm->mapped_[drxn] = drxn ? max( cm->mapped_[1], coord[1] ) : min( cm->mapped_[0], coord[1] );
//            if ( !bubbled && result.s_[0][i] == result.s_[1][i] )
//            {
//                cm->coord_[drxn] = drxn ? max( cm->coord_[1], coord[0] ) : min( cm->coord_[0], coord[0] );
//                cm->range_[drxn] = drxn ? max( cm->range_[1], coord[1] ) : min( cm->range_[0], coord[1] );
//            }
//        }
//        
//        while ( starts.size() < result.bubbles_.size() && i == result.bubbles_[starts.size()].start_ ) starts.push_back( make_pair( coord[0], coord[1] ) );
//        for ( int s : { 0, 1 } ) if ( result.s_[s] != '-' ) coord[s]++;
//    }
//    
//    for ( int b = 0; b < starts.size(); b++ )
//    {
//        coord[0] = starts[b].first;
//        coord[1] = starts[b].second;
//        int j = 0;
//        vector<BubbleAlignCoords*> stack = { &result.bubbles_ };
//        for ( int i = result.bubbles_[b].start_; i < result.bubbles_[b].end_; i++ )
//        {
//            
//        }
//    }
//    
//    vector< pair< vector<SnpAlignResult::BubbleAlignCoords>*, int> > bubbles{ make_pair( &result.bubbles_, 0 ) };
//    vector<SnpAlignResult::BubbleAlignCoords*> stack{ NULL };
//    for ( int i = 0; i < result.s_[0].size(); i++ )
//    {
//        for ( int& j = bubbles.back().second; j < bubbles.back().first->size() && (*bubbles.back().first)[j].start_ >= i; j++ );
////        if ( bubbles.back().second < bubbles.back().first->size() && );
//    }
//    assert( false );
}
