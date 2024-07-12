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

#include "bundle.h"
#include "alignment.h"
#include <algorithm>
#include <cassert>
#include <iostream>

//Bundle::BundleEdge::BundleEdge( BundleNode* l, BundleNode* r )
//: score_( -1 ), errs_( max( l->node_->size(), r->node_->size() ) )
//{
//    start_[0] = start_[1] = end_[0] = end_[1] = 0;
//    
//    ReadMatch* rm[2]{ l->node_, r->node_ };
//    assert( !rm[0]->align_.empty() && !rm[1]->align_.empty() );
//    int p[2]{0}, s[2]{0};
//    while ( p[0] < rm[0]->align_.size() && p[1] < rm[1]->align_.size() )
//    {
//        rm[0]->align_[p[0]];
//        rm[1]->align_[p[1]];
//        size_t tar = max( rm[0]->align_[p[0]]->tar_, rm[1]->align_[p[1]]->tar_ );
//        bool failed[2]{ false, false };
//        for ( int i : { 0, 1 } ) if ( failed[i] = rm[i]->align_[p[i]]->tar_+rm[i]->align_[p[i]]->len_ <= tar ) p[i]++;
//        assert( !failed[0] || !failed[1] );
//        if ( failed[0] || failed[1] ) continue;
//        Align al;
//        al.len_ = min( rm[0]->align_[p[0]]->tar_ + rm[0]->align_[p[0]]->len_ - tar
//                     , rm[1]->align_[p[1]]->tar_ + rm[1]->align_[p[1]]->len_ - tar );
//        for ( int i : { 0, 1 } ) al.s_[i] = rm[i]->align_[p[i]]->read_ + tar - rm[i]->align_[p[i]]->tar_;
//        assert( al.len_ );
//        if ( al.len_ ) align_.push_back( al );
//    }
//    assert( !align_.empty() );
//    for ( int i = 0; i <= align_.size(); i++ )
//    {
//        bool anchored[2]{ i, i < align_.size() };
//        int starts[2]{ s[0], s[1] };
//        int ends[2]{ i < align_.size() ? align_[i].s_[0] : rm[0]->size(), i < align_.size() ? align_[i].s_[1] : rm[1]->size() };
//        assert( starts[0] <= ends[0] && starts[1] <= ends[1] );
//        for ( int j : { 0, 1 } ) s[j] = align_[i].s_[j] + align_[i].len_;
//        if ( starts[0] == ends[0] || starts[1] == ends[1] ) continue;
//        AlignResult ar = Alignment( rm[0]->seq(), rm[1]->seq_, starts[0], ends[0]-starts[0], starts[1], ends[1]-starts[1] ).align( 1, 2, 3, 3, anchored[0], anchored[1]);
//        vector<Align> realigns;
//        for ( int j = 0; j < ar.size(); j++ )
//        {
//            if ( ar.isValid( j ) && ar.isMatch( j ) )
//            {
//                if ( !realigns.empty() && realigns.back().s_[0]+realigns.back().len_ == s[0] && realigns.back().s_[1]+realigns.back().len_ == s[1] ) realigns.back().len_++;
//                else realigns.push_back( Align( s[0], s[1], 1 ) );
//            }
//        }
//        
//        assert( realigns.empty() );
//        if ( !realigns.empty() ) align_.insert( align_.begin() + i, realigns.begin(), realigns.end() );
//        i += realigns.size();
//    }
//    
//    for ( int i = 1; i < align_.size(); i++ )
//    {
//        if ( align_[i-1].s_[0] + align_[i-1].len_ == align_[i].s_[0] && align_[i-1].s_[1] + align_[i-1].len_ == align_[i].s_[1] )
//        {
//            align_[i-1].len_ += align_[i].len_;
//            align_.erase( align_.begin() + i-- );
//        }
//    }
//    
//    if ( align_.empty() ) return;
//    
//    int score = score_ = align_[0].len_, errs = errs_ = 0, start[2];
//    for ( int s : { 0, 1 } ) start_[s] = start[s] = align_[0].s_[s];
//    for ( int s : { 0, 1 } ) end_[s] = align_[0].s_[s] + align_[0].len_;
//    
//    
//    for ( int i = 1; i < align_.size(); i++ )
//    {
//        int gap = max( align_[i].s_[0]-end[0], align_[0].s_[1]-end[1] );
//        errs += gap;
//        assert( gap >= 1 );
//        score -= gap * 4;
//        if ( score < 0 ) for ( int s : { 0, 1 } ) start[s] = align_[i].s_[s];
//        score = max( score, 0 ) + align_[i].len_;
//        if ( score > score_ )
//        {
//            for ( int s : { 0, 1 } ) start_[s] = start[s];
//            for ( int s : { 0, 1 } ) end_[s] = align_[0].s_[s] + align_[0].len_;
//            errs_ = errs;
//            score_ = score;
//        }
//        if ( score > score_ ) score_ = score;
//    }
//    
//}

vector<Bundle::BBubble*> Bundle::BBubble::getBubbles( vector<BundleNode*>& nodes )
{
    vector<BBubble*> bubs;
    for ( BundleNode* node : nodes )
    {
        ReadMatch* rm  = node->node_;
        for ( int i = 1; i < rm->align_.size(); i++ )
        {
            string seq = rm->seq_.substr( rm->align_[i-1].readEnd(), rm->align_[i].read_-rm->align_[i-1].readEnd() );
            size_t coord[2]{ rm->align_[i-1].tarEnd(), rm->align_[i].tar_ };
            bool added = false;
            BBubble* bub = NULL;
            for ( int j = 0; !added && j <= bubs.size(); j++ )
            {
                if ( j >= bubs.size() || ( coord[0] == bubs[j]->start_ && coord[1] > bubs[j]->end_ ) || coord[0] < bubs[j]->start_ )
                {
                    bubs.insert( bubs.begin()+j, bub = new BBubble( seq, rm, coord[0], coord[1], 2 ) );
                }
                else if ( bubs[j]->setMatch( seq, rm, coord[0], coord[1], 2 ) ) bub = bubs[j];
            }
            assert( bub );
            node->bubs_.push_back( bub );
        }
    }
    
    for ( BundleNode* node : nodes )
    {
        ReadMatch* rm  = node->node_;
        for ( BBubble* b : bubs  )
        {
            if ( rm->tar_[1] <= b->start_ ) break;
            if ( rm->tar_[0] < b->end_ ) b->setNonmatch( rm );
        }
    }
    
    return bubs;
}

bool Bundle::BBubble::setMatch( string& seq, ReadMatch* rm, size_t start, size_t end, size_t drxn )
{
    if ( start != start_ || end != end_ || drxn != drxn_ || seq != seq_ ) return false;
    matches_.push_back( rm );
    return true;
}

bool Bundle::BBubble::setNonmatch( ReadMatch* rm )
{
    if ( find( matches_.begin(), matches_.end(), rm ) != matches_.end() ) return false;
    matches_.push_back( rm );
    return true;
}

//void Bundle::add( ReadMatch* rm, vector<BBubble*>& bubs )
//{
//    for ( int i = 1; i < rm->align_.size(); i++ )
//    {
//        string seq = rm->seq_.substr( rm->align_[i-1].readEnd(), rm->align_[i].read_-rm->align_[i-1].readEnd() );
//        size_t coord[2]{ rm->align_[i-1].tarEnd(), rm->align_[i].tar_ };
//        bool added = false;
//        for ( int j = 0; !added && j <= bubs.size(); j++ )
//        {
//            if ( added = ( j >= bubs.size() || ( coord[0] == bubs[j]->start_ && coord[1] > bubs[j]->end_ ) || coord[0] < bubs[j]->start_ ) )
//            {
//                bubs.insert( bubs.begin()+j, new BBubble( seq, rm, coord[0], coord[1], 2 ) );
//            }
//            else added = bubs[j]->setMatch( seq, rm, coord[0], coord[1], 2 );
//        }
//    }
//}

Bundle::BundleNode::BundleNode( ReadMatch* rm )
: node_( rm ), score_( 0 )
{
    tar_[0] = rm->align_[0].tar_;
    tar_[1] = rm->align_.back().tarEnd();
    for ( int i = 0; i < rm->align_.size(); i++ )
    {
        score_ += rm->align_[i].len_;
        if ( i ) score_ -= rm->align_[i-1].gapLen( rm->align_[i] ) * 4;
    }
    score_ -= ( rm->align_[0].read_ + ( rm->size() - rm->align_.back().readEnd() ) ) * 2;
}

Bundle::PPath::PPath( BundleNode* base, vector<BundleNode*>& nodes, vector<BBubble*>& bubs )
{
    vector<BundleNode*> extNodes[2];
    vector<BBubble*> extBubs[2];
    for ( BundleNode* node : nodes ) if ( node->tar_[0] < base->tar_[0] ) extNodes[0].push_back( node );
    for ( BundleNode* node : nodes ) if ( node->tar_[1] > base->tar_[1] ) extNodes[1].push_back( node );
    for ( BBubble* b : bubs ) if ( b->end_ <= base->tar_[0] ) extBubs[0].push_back( b );
    for ( BBubble* b : bubs ) if ( b->start_ >= base->tar_[1] ) extBubs[1].push_back( b );
    reverse( extNodes[0].begin(), extNodes[0].end() );
    reverse( extBubs[0].begin(), extBubs[0].end() );
    sort( extNodes[0].begin(), extNodes[0].end(), []( BundleNode* a, BundleNode* b ){ return a->tar_[1] == b->tar_[1] ? a->tar_[0] < b->tar_[0] : a->tar_[1] > b->tar_[1]; } );
    sort( extNodes[1].begin(), extNodes[1].end(), []( BundleNode* a, BundleNode* b ){ return a->tar_[0] == b->tar_[0] ? a->tar_[1] > b->tar_[1] : a->tar_[0] < b->tar_[0]; } );
    sort( extBubs[0].begin(), extBubs[0].end(), []( BBubble* a, BBubble* b ){ return a->end_ == b->end_ ? a->start_ < b->start_ : a->end_ > b->end_; } );
    sort( extBubs[1].begin(), extBubs[1].end(), []( BBubble* a, BBubble* b ){ return a->start_ == b->start_ ? a->end_ > b->end_ : a->start_ < b->start_; } );
    
    
//    while ( extend( nodes, bubs ) );
    // Add specific read when adding bubble = most supported, least controversial
}

bool Bundle::PPath::extend( vector<BundleNode*>& nodes, vector<BBubble*>& bubs, bool drxn )
{
//    vector<Branch> branches;
//    for ( BBubble* bub : bubs ) if ( !bub->matches_.empty() && !bub->nonmatches_.empty() )
//    {
//        vector<BundleNode*> overlaps = bub->matches_;
//        overlaps.insert( overlaps.end(), bub->nonmatches_.begin(), bub->nonmatches_.end() );
//        sort( overlaps.begin(), overlaps.end(), [&]( BundleNode* a, BundleNode* b ){ return drxn ? a->tar_[0] < b->tar_[0] : a->tar_[1] > b->tar_[1]; } );
//        int i = 0;
//        for ( int j = i; j < overlaps.size(); )
//        {
//            while ( j < overlaps.size() && ) j++;
//            sort( overlaps.begin()+i, overlaps.begin()+j, [&]( BundleNode* a, BundleNode* b ){ return drxn ? a->tar_[0] < b->tar_[0] : a->tar_[1] > b->tar_[1]; } );
//            for ( int k = i; k < j; k++ )
//            {
//                for ( Branch branch : branches )
//                {
//                    
//                }
//            }
//        }
//    }
//    return false;
}

//Bundle::BundleRead::BundleRead( ReadMatch* rm, vector<Node*>& bubs )
//: read_( rm )
//{
//    for ( int i = 0; i <= rm->align_.size(); i++ ) if ( ( i || rm->align_[0].read_ ) && ( i < rm->align_.size() || rm->align_.back().readEnd() < rm->size() ) )
//    {
//        int type = !i ? 0 : ( i < rm->align_.size() ? 2 : 1 );
//        size_t read[2]{ i ? rm->align_[i-1].readEnd() : 0, i < rm->align_.size() ? rm->align_[i].read_ : rm->size() };
//        size_t tar[2]{ i ? rm->align_[i-1].tarEnd() : rm->align_[0].tar_, i < rm->align_.size() ? rm->align_[i].tar_ : rm->align_[i-1].tarEnd() };
//        string seq = rm->seq_.substr( read[0] , read[1] - read[0] );
//        
//        int j = bubs.size();
//        Node* node = NULL;
//        for ( ; !node && j > 0 && bubs[j-1]->coord_[0] == tar[0] ? bubs[j-1]->coord_[1] < tar[1] : bubs[j-1]->coord_[0] > tar[0]; j-- )
//        {
//            if ( type != 0 && bubs[j-1]->coord_[0] != tar[0] ) continue;
//            if ( type != 1 && bubs[j-1]->coord_[1] != tar[1] ) continue;
//            if ( seq != bubs[j-1]->seq_ ) continue;
//            node = bubs[j-1];
//            node->reads_.push_back( rm );
//        }
//        if ( !node ) bubs.insert( bubs.begin() + j, node = new Node( seq, rm, tar, type ) );
//        nodes_.push_back( node );
//    }
//}

Bundle::BundleRead::BundleRead( ReadMatch* rm, vector<Node*>& base, vector<Node*>& bubs, int& p )
: read_( rm ), score_( 0 )
{
    assert( rm->align_.size() );
    tar_[0] = rm->align_[0].tar_;
    tar_[1] = rm->align_.back().tarEnd();
    while( p < base.size() && base[p]->coord_[1] <= rm->align_[0].tar_ ) p++;
    int i = p;
    
    Node* cur = NULL;
    size_t read = 0, tar = rm->align_[0].tar_;
    for ( ReadMatch::Align& al : rm->align_ )
    {
        if ( al.read_ ) nodes_.push_back( cur = Node::setBubble( bubs, cur, this, rm->seq_.substr( read, al.read_-read ), tar, al.tar_, read ? 2 : 0 ) );
        
        for ( ; i < base.size() && base[i]->coord_[0] < al.tarEnd(); i++ ) if ( al.tar_ < base[i]->coord_[1] )
        {
            base[i]->reads_.push_back( this );
            Edge::thread( cur, base[i], this );
            cur = base[i];
        }
        read = al.readEnd();
        tar = al.tarEnd();
    }
    if ( read < rm->size() ) nodes_.push_back( Node::setBubble( bubs, cur, this, rm->seq_.substr( read, rm->size()-read ), tar, tar, 1 ) );
}

int Bundle::BundleRead::getMismatches()
{
    int mismatches = 0;
    for ( Node* node : nodes_ ) mismatches += node->mismatches_;
    return mismatches;
}

bool Bundle::BundleRead::isGood( size_t errorRate )
{
    int acceptable = ( read_->size() * errorRate ) / 100;
    return getMismatches() <= acceptable;
}

vector<Bundle::Node*> Bundle::BundleRead::getPath( vector<Node*>& base )
{
    vector<Node*> path = !nodes_.empty() && nodes_[0]->coord_[1] == tar_[0] ? vector<Node*>( 1, nodes_[0] ) : vector<Node*>();
    int i = 0, j = path.size(), tar = tar_[0];
    while ( i < base.size() && base[i]->coord_[1] <= tar_[0] ) i++;
    while ( ( i < base.size() && base[i]->coord_[0] < tar_[1] ) || ( j < nodes_.size() ) )
    {
        if ( j < nodes_.size() && tar == nodes_[j]->coord_[0] )
        {
            tar = nodes_[j]->coord_[1];
            path.push_back( nodes_[j++] );
            while ( i < base.size() && base[i]->coord_[0] < tar ) i++;
        }
        else
        {
            assert( i < base.size() && ( path.empty() ? base[i]->coord_[0] <= tar : base[i]->coord_[0] == tar ) );
            tar = base[i]->coord_[1];
            path.push_back( base[i++] );
        }
    }
    assert( !path.empty() );
    assert( path[0]->coord_[0] <= tar_[0] && tar_[0] <= path[0]->coord_[1] );
    assert( path.back()->coord_[0] <= tar_[1] && tar_[1] <= path.back()->coord_[1] );
    for ( Node* node : nodes_ ) assert( find( path.begin(), path.end(), node ) != path.end() );
    for ( int k = 1; k < path.size(); k++ ) assert( path[k-1]->coord_[1] == path[k]->coord_[0] );
    return path;
}

int Bundle::BundleRead::getRemainder( Edge* edge, bool drxn )
{
    for ( int d : { 0, 1 } ) if ( edge->edge_[d]->type_ == -1 )
    {
        for ( ReadMatch::Align& al : read_->align_ ) if ( d ? al.tar_ == edge->edge_[1]->coord_[0] : al.tarEnd() == edge->edge_[0]->coord_[1] )
        {
            int coord = d ? al.read_ : al.readEnd();
            return drxn ? read_->size() - coord : coord;
        }
        assert( false );
    }
    assert( false );
    return 0;
}

bool Bundle::BundleRead::isPerfect()
{
    return !getMismatches();
}

void Bundle::BundleRead::sort( vector<BundleRead*>& reads )
{
    std::sort( reads.begin(), reads.end(), []( BundleRead* a, BundleRead* b ){ return a->tar_[0] == b->tar_[0] ? a->tar_[1] > b->tar_[1] : a->tar_[0] < b->tar_[0];} );
}

Bundle::Node::Node( string& seq, Node* lEdge, BundleRead* read, size_t lCoord, size_t rCoord, int type )
: seq_( type == -1 ? seq.substr( lCoord, rCoord-lCoord ) : seq ), mismatches_( 0 ), type_( type )
{
    if ( lEdge ) Edge::thread( lEdge, this, read );
    if ( read ) reads_.push_back( read );
    coord_[0] = lCoord;
    coord_[1] = rCoord;
}

Bundle::Node::~Node()
{
    for ( int d : { 0, 1 } ) for ( Edge* edge : edges_[d] )
    {
        edge->edge_[!d] = NULL;
        if ( !edge->edge_[d] ) delete edge;
    }
}

bool Bundle::Node::addBubble( string& seq, BundleRead* read, size_t lCoord, size_t rCoord, int type )
{
    if ( type != type_ || lCoord != coord_[0] || rCoord != coord_[1] || seq != seq_ ) return false;
    reads_.push_back( read );
    return true;
}

Bundle::Node* Bundle::Node::setBubble( vector<Node*>& bubs, Node* lEdge, BundleRead* read, string seq, size_t lCoord, size_t rCoord, int type )
{
    int i = bubs.size();
    for ( ; i > 0 && ( bubs[i-1]->coord_[0] == lCoord ? bubs[i-1]->coord_[1] < rCoord : bubs[i-1]->coord_[0] > lCoord ); i-- )
    {
        if ( bubs[i-1]->addBubble( seq, read, lCoord, rCoord, type ) ) 
        {
            Edge::thread( lEdge, bubs[i-1], read );
            return bubs[i-1];
        }
    }
    
    bubs.insert( bubs.begin() + i, new Node( seq, lEdge, read, lCoord, rCoord, type ) );
    return bubs[i];
}

void Bundle::Node::setMismatches( string& base )
{
    mismatches_ = 0;
    if ( type_ == -1 || coord_[1] == 0 || coord_[0] == base.size() ) return;
    
    int coord[2]{ type_ != 0 ? (int)coord_[0] : max( 0, (int)coord_[1]-(int)seq_.size()-5 ), type_ != 1 ? (int)coord_[1] : (int)min( base.size(), coord_[0]+seq_.size()+5 ) };
    AlignResult ar = Alignment( base, seq_, coord[0], coord[1]-coord[0], 0, seq_.size() ).align( 1, 2, 3, 3, type_ != 0, type_ != 1 );
    
    mismatches_ = ar.getScore( 1, 2, 3, 3, type_ ).miss_;
}

//void Bundle::Node::setNodes( string& seq, vector<Node*>& base, vector<Node*>& bubbles, vector<ReadMatch*>& matches, vector<BundleRead*>& reads )
//{
//    size_t limits[2]{ matches[0]->align_[0].tar_, matches.back()->align_.back().tarEnd() };
//    unordered_set<size_t> joins;
//    for ( ReadMatch* rm : matches )
//    {
//        limits[0] = min( limits[0], rm->align_[0].tar_ );
//        limits[1] = max( limits[1], rm->align_.back().tarEnd() );
//        for ( Align& al : rm->align_ ) if ( al.read_ ) joins.insert( al.tar_ );
//        for ( Align& al : rm->align_ ) if ( al.readEnd() < rm->size() ) joins.insert( al.readEnd() );
//    }
//    joins.insert( limits[0] );
//    joins.insert( limits[1] );
//    
//    vector<size_t> breaks( joins.begin(), joins.end() );
//    sort( breaks.begin(), breaks.end() );
//    for ( int i = 1; i < breaks.size(); i++ ) base.push_back( new Node( seq, base, breaks[i-1], breaks[i] ) );
//    
//    int p = 0;
//    for ( ReadMatch* rm : matches ) reads.push_back( new BundleRead( rm, base, bubbles, p ) );
//}

//void Bundle::Node::setBubbles( vector<Node*>& bubbles )
//{
//    
//}

Bundle::Edge::Edge( Node* lNode, Node* rNode )
{
    edge_[0] = lNode;
    edge_[1] = rNode;
    for ( int d : { 0, 1 } ) edge_[d]->edges_[!d].push_back( this );
}

bool Bundle::Edge::thread( Node* l, Node* r, BundleRead* read )
{
    Edge* edge = NULL;
    if ( l ) for ( int i = 0; !edge && i < l->edges_[1].size(); i++ ) if ( l->edges_[1][i]->edge_[1] == r ) edge = l->edges_[1][i];
    if ( !edge && l && r ) edge = new Edge( l, r );
    if ( edge && read )
    {
        assert( find( edge->reads_.begin(), edge->reads_.end(), read ) == edge->reads_.end() );
        edge->reads_.push_back( read );
    }
}

Bundle::Path::Extend::Extend( BundleRead* read, Edge* edge, bool drxn )
: read_( read )
{
    len_[0] = len_[1] = 0;
    len_[2] = read->read_->size();
    if ( edge )
    {
        ext_[drxn].push_back( edge->edge_[drxn] );
        len_[drxn] = read->getRemainder( edge, drxn );
        len_[2] -= len_[drxn];
    }
}

Bundle::Path::Extend::Extend( Extend* ext, unordered_set<BundleRead*> reads )
: reads_( reads )
{
    if ( ext ) for ( int d : { 0, 1 } ) ext_[d] = ext->ext_[d];
    if ( ext ) for ( BundleRead* br : reads ) ext->reads_.erase( br );
}

//void Bundle::Path::Extend::init( vector<Node*>& nodes, vector<Extend*>& ext, vector<Extend*>& splinters, bool drxn )
//{
//    Node* prv = NULL,* cur = NULL,* nxt = drxn ? nodes.back() : nodes[0];
//    unordered_map<BundleRead*, Extend*> extends;
//    unordered_set<BundleRead*> ignore;
//    for ( Edge* e : nxt->edges_[drxn] ) for ( BundleRead* br : e->reads_ ) extends.insert( make_pair( br, new Extend( br, e, drxn ) ) );
//    
//    int dist = 0;
//    for ( int i = 0; i < nodes.size(); i++ )
//    {
//        cur = nxt;
//        nxt = i+1 < nodes.size() ? ( d ? nodes.end()[-i-2] : nodes[i+1] ) : NULL;
//        
//        if ( prv ) for ( Edge* e : cur->edges_[drxn] ) if ( e->edge_[drxn] != prv ) ignore.insert( e->reads_.begin(), e->reads_.end() );
//        
//        for ( BundleRead* br : cur->reads_ ) if ( ignore.find( br ) == ignore.end() && extends.find( br ) == extends.end() ) extends.insert( make_pair( br, new Extend( br, NULL, drxn ) ) );
//        dist += cur->seq_.size();
//        
//        for ( Edge* e : cur->edges_[!drxn] ) if ( e->edge_[!drxn] != nxt ) for ( BundleRead* br : e->reads_ ) if ( ignore.find( br ) == ignore.end() )
//        {
//            auto it = extends.find( br );
//            if ( it == extends.end() ) continue;
//        }
//        if ( nxt );
//        
//        prv = cur;
//    }
//}

bool Bundle::Path::Extend::shrink( Node* node, bool drxn )
{
    if ( ext_[drxn].empty() ) return true;
    if ( ext_[drxn][0] == node )
    {
        ext_[drxn].erase( ext_[drxn].begin() );
        if ( !ext_[0].empty() || !ext_[1].empty() ) return true;
    }
    delete this;
    return false;
}

bool Bundle::Path::Extend::splinter( Edge* edge, int dist, bool drxn )
{
    assert( ext_[drxn].empty() );
    len_[drxn] = read_->getRemainder( edge, drxn );
    len_[2] -= len_[drxn];
    if ( ext_[!drxn].empty() ) len_[!drxn] = dist - len_[drxn];
    assert( ext_[!drxn].empty() ? len_[!drxn] <= 0 : len_[!drxn] > 0 );
    ext_[drxn].push_back( edge->edge_[drxn] );
    return true;
}

Bundle::Path::Path( vector<Node*> path, BundleRead* br )
: nodes_( path )
{
    coords_[0] = br->tar_[0];
    coords_[1] = br->tar_[1];
    
    addReads();
    
//    unordered_set<BundleRead*> miss = getTangents( true ), added;
//    for ( int i = 0; i < path.size(); i++ )
//    {
//        vector<BundleRead*> reads;
//        for ( BundleRead* read : path[i]->reads_ ) if ( miss.find( read ) == miss.end() ) reads.push_back( read );
//        if ( i && i+1 < path.size() ) for ( BundleRead* read : reads ) if ( added.insert( read ).second ) reads_.push_back( read );
//        if ( !i )
//        {
//            sort( reads.begin(), reads.end(), []( BundleRead* a, BundleRead* b ){ return a->tar_[1] == b->tar_[1] ? a->tar_[0] > b->tar_[0] : a->tar_[1] < b->tar_[1]; } );
//            for ( int j = reads.size(); j-- > 0 && coords_[0] < reads[j]->tar_[1]; )
//            {
//                coords_[0] = min( coords_[0], reads[j]->tar_[0] );
//                if ( added.insert( reads[j] ).second ) reads_.push_back( reads[j] );
//            }
//        }
//        if ( i+1 == path.size() )
//        {
//            sort( reads.begin(), reads.end(), []( BundleRead* a, BundleRead* b ){ return a->tar_[0] == b->tar_[0] ? a->tar_[1] > b->tar_[1] : a->tar_[0] < b->tar_[0]; } );
//            for ( int j = 0; j < reads.size() && reads[j]->tar_[0] < coords_[1]; j++ )
//            {
//                coords_[1] = max( coords_[1], reads[j]->tar_[1] );
//                if ( added.insert( reads[j] ).second ) reads_.push_back( reads[j] );
//            }
//        }
//    }
    assert( !reads_.empty() );
}

void Bundle::Path::addReads()
{
    unordered_set<BundleRead*> miss = getTangents( true ), added( reads_.begin(), reads_.end() );
    vector<BundleRead*> reads;
    for ( Node* node : nodes_ ) for ( BundleRead* br : node->reads_ ) if ( miss.find( br ) == miss.end() && added.insert( br ).second ) reads.push_back( br );
    
    sort( reads.begin(), reads.end(), []( BundleRead* a, BundleRead* b ){ return a->tar_[1] == b->tar_[1] ? a->tar_[0] > b->tar_[0] : a->tar_[1] < b->tar_[1]; } );
    for ( int i = reads.size(); i-- > 0 && coords_[0] < reads[i]->tar_[1]; ) coords_[0] = min( coords_[0], reads[i]->tar_[0] );
    sort( reads.begin(), reads.end(), []( BundleRead* a, BundleRead* b ){ return a->tar_[0] == b->tar_[0] ? a->tar_[1] > b->tar_[1] : a->tar_[0] < b->tar_[0]; } );
    for ( int i = 0; i < reads.size() && reads[i]->tar_[0] < coords_[1]; i++ )  coords_[1] = max( coords_[1], reads[i]->tar_[1] );
    
    for ( BundleRead* br : reads ) if ( br->tar_[0] < coords_[1] && coords_[0] < br->tar_[1] ) reads_.push_back( br );
}

void Bundle::Path::evaluate( Branch& branch, bool drxn )
{
    
}

void Bundle::Path::extend()
{
    for ( int d : { 0, 1 } )
    {
        Node* cur,* nxt = d ? nodes_.back() : nodes_[0];
        for ( int i = 0; i < nodes_.size(); i++ )
        {
            cur = nxt;
            nxt = i+1 < nodes_.size() ? ( d ? nodes_.end()[-i-2] : nodes_[i+1] ) : NULL;
        }
    }
//    unordered_set<BundleRead*> incl, miss = getTangents();
//    for ( int d : { 0, 1 } ) for ( BundleRead* br : getExtending( miss, d ) ) incl.insert( br );
//    if ( incl.empty() ) return;
//    
//    vector<Extend*> extends = { new Extend( NULL, incl ) };
//    setExtend( extends[0], NULL, extends, 0 );
//    
    unordered_set<BundleRead*> miss = getTangents();
    vector<Extend*> extends = getExtends( miss );
    vector<Branch> branches[2]{ setBranches( extends, 0 ), setBranches( extends, 1 ) };
    for ( bool extended = true; extended && !( extended=false ) && ( branches[0].size() + branches[1].size() ); )
    {
        while ( branches[0].size() == 1 || branches[1].size() == 1  )
        {
            for ( int d : { 0, 1 } ) while ( branches[d].size() == 1 ) extended = extend( branches[d][0], branches, extends, miss, d );
        }
        for ( int d : { 0, 1 } ) if ( extended = resolve( branches[d], branches[!d], extends, d ) ) break;
        assert( branches[0].empty() && branches[1].empty() );
    }
}

bool Bundle::Path::extend( Branch& branch, vector<Branch> (&branches)[2], vector<Extend*>& extends, unordered_set<BundleRead*>& miss, bool drxn )
{
    Node* fork = drxn ? nodes_.back() : nodes_[0];
    for ( Edge* e : fork->edges_[drxn] ) if ( e->edge_[drxn] != branch.node_ ) miss.insert( e->reads_.begin(), e->reads_.end() );
    for ( Edge* e : branch.node_->edges_[!drxn] ) if ( e->edge_[!drxn] != fork ) miss.insert( e->reads_.begin(), e->reads_.end() );
    
    nodes_.insert( drxn ? nodes_.end() : nodes_.begin(), branch.node_ );
    addReads();
    for ( int i = 0; i < extends.size(); i++ ) if ( !extends[i]->shrink( branch.node_, drxn ) ) extends.erase( extends.begin() + i-- );
    
    unordered_set<BundleRead*> used( miss.begin(), miss.end() ), added;
    for ( Extend* ext : extends ) used.insert( ext->reads_.begin(), ext->reads_.end() );
    for ( Edge* e : branch.node_->edges_[drxn] ) for ( BundleRead* br : e->reads_ ) if ( used.find( br ) == used.end() ) added.insert( br );
    
    if ( !added.empty() )
    {
        extends.push_back( new Extend( NULL, added ) );
        setExtend( extends.back(), NULL, extends, 0 );
        assert( false );
    }
    
    for ( int d : { 0, 1 } ) branches[d] = setBranches( extends, d );
    return true;
}

void Bundle::Path::initExtends( vector<Extend*>& exts, vector<Extend*>& splinters, bool drxn )
{
    Node* prv = NULL,* cur = NULL,* nxt = drxn ? nodes_.back() : nodes_[0];
    unordered_map<BundleRead*, Extend*> extends;
    unordered_set<BundleRead*> ignore( reads_.begin(), reads_.end() );
    for ( Edge* e : nxt->edges_[drxn] ) for ( BundleRead* br : e->reads_ ) extends.insert( make_pair( br, new Extend( br, e, drxn ) ) );
    
    int dist = 0;
    for ( int i = 0; i < nodes_.size(); i++ )
    {
        cur = nxt;
        nxt = i+1 < nodes_.size() ? ( drxn ? nodes_.end()[-i-2] : nodes_[i+1] ) : NULL;
        
        if ( prv ) for ( Edge* e : cur->edges_[drxn] ) if ( e->edge_[drxn] != prv ) ignore.insert( e->reads_.begin(), e->reads_.end() );
        
        for ( BundleRead* br : cur->reads_ ) if ( ignore.find( br ) == ignore.end() && extends.find( br ) == extends.end() ) extends.insert( make_pair( br, new Extend( br, NULL, drxn ) ) );
        dist += cur->seq_.size();
        
        for ( Edge* e : cur->edges_[!drxn] ) if ( e->edge_[!drxn] != nxt ) for ( BundleRead* br : e->reads_ ) if ( !ignore.insert( br ).second )
        {
            auto it = extends.find( br );
            if ( it == extends.end() ) continue;
            it->second->splinter( e, dist, !drxn );
        }
        
        prv = cur;
    }
}

vector<Bundle::Path::Extend*> Bundle::Path::getExtends( unordered_set<BundleRead*>& miss )
{
    vector<Extend*> extends;
    unordered_set<BundleRead*> incl;
    for ( int d : { 0, 1 } )
    {
        Node* node = d ? nodes_.back() : nodes_[0];
        for ( Edge* e : node->edges_[d] ) for ( BundleRead* read : e->reads_ ) if ( miss.find( read ) == miss.end() ) incl.insert( read );
    }
    
    if ( !incl.empty() ) extends.push_back( new Extend( NULL, incl ) );
    if ( !extends.empty() ) setExtend( extends[0], NULL, extends, 0 );
    
    return extends;
}

//vector<Bundle::BundleRead*> Bundle::Path::getExtending( unordered_set<BundleRead*>& ignore, bool drxn )
//{
//    vector<BundleRead*> reads;
//    Node* node = drxn ? nodes_.back() : nodes_[0];
//    for ( Edge* e : node->edges_[drxn] ) for ( BundleRead* read : e->reads_ ) if ( ignore.find( read ) == ignore.end() ) reads.push_back( read );
//    return reads;
//}

unordered_set<Bundle::BundleRead*> Bundle::Path::getTangents( bool inclEnds )
{
    unordered_set<BundleRead*> miss;
    for ( int i = 0; i < nodes_.size(); i++ )
    {
        Node* nxt[2]{ i ? nodes_[i-1] : NULL, i+1 < nodes_.size() ? nodes_[i+1] : NULL };
        for ( int d : { 0, 1 } ) if ( inclEnds || nxt[d] ) for ( Edge* edge : nodes_[i]->edges_[d] ) if ( edge->edge_[d] != nxt[d] ) miss.insert( edge->edge_[d]->reads_.begin(), edge->edge_[d]->reads_.end() );
    }
    return miss;
}

bool Bundle::Path::resolve( vector<Branch>& fwd, vector<Branch>& rev, vector<Extend*>& extends, bool drxn )
{
    assert( fwd.size() != 1 );
    if ( fwd.size() < 2 ) return false;
    
    unordered_set<BundleRead*> ignore( reads_.begin(), reads_.end() ), reads;
    unordered_map<BundleRead*, Edge*> alts;
    int ol = 0, maxOl = 0;
    for ( Branch b : fwd )
    {
        for ( Extend* ext : b.extends_ ) ignore.insert( ext->reads_.begin(), ext->reads_.end() );
        for ( BundleRead* br : b.edge_->reads_ ) if ( ignore.insert( br ).second ) assert( false );
        for ( BundleRead* br : b.edge_->reads_ ) if ( ignore.insert( br ).second && reads.insert( br ).second ) alts.insert( make_pair( br, b.edge_ ) );
        maxOl = max( maxOl, b.extends_[0]->ol_[drxn] );
    }
    Node* nxt = drxn ? nodes_.back() : nodes_[0],* cur = NULL,* prv = NULL;
    for ( int i = 0; i < nodes_.size(); i++ )
    {
        prv = cur;
        cur = nxt;
        nxt = i+1 < nodes_.size() ? ( drxn ? nodes_.end()[-i-2] : nodes_[i+1] ) : NULL;
        ol += cur->seq_.size();
        
        if ( i ) for ( Edge* e : cur->edges_[drxn] ) if ( e->edge_[drxn] != prv ) ignore.insert( e->reads_.begin(), e->reads_.end() );
        for ( Edge* e : cur->edges_[!drxn] ) if ( e->edge_[!drxn] != nxt )
        {
            unordered_set<BundleRead*> reads;
            for ( BundleRead* read : e->reads_ ) if ( ignore.find( read ) == ignore.end() ) reads.insert( read );
            assert( reads.empty() );
        }
    }
    
    assert( false );
    return true;
}

vector<Bundle::Path::Branch> Bundle::Path::setBranches( vector<Extend*>& extends, bool drxn )
{
    vector<Branch> branches;
    if ( extends.empty() ) return branches;
    setOverlaps( extends );
    for ( Edge* edge : ( drxn ? nodes_.back() : nodes_[0] )->edges_[drxn] )
    {
        Branch branch;
        branch.edge_ = edge;
        branch.node_ = edge->edge_[drxn];
        for ( Extend* ext : extends ) if ( !ext->ext_[drxn].empty() && ext->ext_[drxn][0] == branch.node_ ) branch.extends_.push_back( ext );
        sort( branch.extends_.begin(), branch.extends_.end(), [&]( Extend* a, Extend* b ){ return a->ol_[drxn] > b->ol_[drxn]; } );
        if ( !branch.extends_.empty() && branch.extends_[0]->ol_[drxn] > 0 ) branches.push_back( branch );
    }
    return branches;
}

void Bundle::Path::setExtend( Extend* ext, Edge* edge, vector<Extend*>& extends, bool drxn )
{
    if ( edge ) ext->ext_[drxn].push_back( edge->edge_[drxn] );
    Node* cur = edge ? edge->edge_[drxn] : ( drxn ? nodes_.back() : nodes_[0] );
    
    bool ended = true;
    for ( Edge* e : cur->edges_[drxn] )
    {
        unordered_set<BundleRead*> reads;
        for ( BundleRead* read : e->reads_ ) if ( ext->reads_.find( read ) != ext->reads_.end() ) reads.insert( read );
        if ( reads.size() == ext->reads_.size() )
        {
            assert( ended );
            ended = false;
            setExtend( ext, e, extends, drxn );
        }
        else if ( !reads.empty() )
        {
            extends.push_back( new Extend( ext, reads ) );
            setExtend ( extends.back(), e, extends, drxn );
        }
    }
    
    if ( ended ) assert( !ext->reads_.empty() );
    if ( ended && !drxn ) setExtend( ext, NULL, extends, 1 );
}

void Bundle::Path::setOverlaps( vector<Extend*>& extends )
{
    for ( Extend* ext : extends ) for ( int d : { 0, 1 } )
    {
        ext->ol_[d] = 0;
        bool set = false;
        if ( !ext->ext_[d].empty() ) for ( BundleRead* br : ext->reads_ )
        {
            int ol = d ? (int)coords_[1] - (int)br->tar_[0] : (int)br->tar_[1] - (int)coords_[0];
            if ( !set || ol > ext->ol_[d] ) ext->ol_[d] = ol;
            set = true;
        }
    }
}

Bundle::Bundle( BundleNode* bn, vector<BundleNode*>& nodes, vector<BBubble*>& bubs )
{
    PPath bp( bn, nodes, bubs );
}

vector<Bundle*> Bundle::bundle( string& seq, LongResult* lr )
{
    vector<Node*> base, bubbles;
    vector<BundleRead*> reads = Bundle::setGraph( seq, base, bubbles, lr->reads_ );
    vector<Bundle*> bundles;
    Bundle::getBase( bundles, base, bubbles, reads );
//    while ( !nodes.empty() )
//    {
//        BundleNode* node = nodes[0];
//        for ( int i = 1; i < nodes.size(); i++ ) if ( nodes[i]->score_ > node->score_ ) node = nodes[i];
//        
//        Bundle* bundle = new Bundle( node, nodes, bubs );
//    }
    for ( Node* node : base ) delete node;
    for ( Node* node : bubbles ) delete node;
    for ( BundleRead* read : reads ) delete read;
    return bundles;
//    vector<BundleNode*> nodes;
//    for ( ReadMatch* rm : lr.reads_ ) nodes.push_back( new BundleNode( rm ) );
    
}

void Bundle::build( BundleRead* read, vector<Node*>& base, unordered_set<BundleRead*>& used )
{
    if ( used.find( read ) != used.end() ) return;
    if ( read->read_->id_ == 4764655 )
    {
        int x = 0;
    }
    Path path( read->getPath( base ), read );
    path.extend();
}

void Bundle::getBase( vector<Bundle*>& bundles, vector<Node*>& base, vector<Node*>& bubbles, vector<BundleRead*>& reads )
{
    vector<BundleRead*> perfect, good, bad;
    for ( BundleRead* read : reads )
    {
        if ( read->isPerfect() ) perfect.push_back( read );
        else if ( read->isGood( 2 ) ) good.push_back( read );
        else bad.push_back( read );
    }
    unordered_set<BundleRead*> used;
    for ( BundleRead* read : perfect ) build( read, base, used );
    for ( BundleRead* read : good ) build( read, base, used );
    
    cout << "p: " << perfect.size() << "  g: " << good.size() << "  b: " << bad.size() << endl;
    int x = 0;
}

vector<Bundle::BundleRead*> Bundle::setGraph( string& seq, vector<Node*>& base, vector<Node*>& bubbles, vector<ReadMatch*>& matches )
{
    unordered_set<size_t> joins;
    vector<pair<size_t, size_t>> limits{ make_pair( matches[0]->align_[0].tar_, matches[0]->align_[0].tarEnd()) };
    sort( matches.begin(), matches.end(), []( ReadMatch* a, ReadMatch* b ){ return a->align_[0].tar_ == b->align_[0].tar_ ? a->align_.back().tarEnd() > b->align_.back().tarEnd() : a->align_[0].tar_ < b->align_[0].tar_; } );
    for ( ReadMatch* rm : matches )
    {
        if ( rm->align_[0].tar_ < limits.back().second )
        {
            limits.back().first = min( limits.back().first, rm->align_[0].tar_ );
            limits.back().second = max( limits.back().second, rm->align_.back().tarEnd() );
        }
        else limits.push_back( make_pair( rm->align_[0].tar_, rm->align_.back().tarEnd() ) );
        
        for ( ReadMatch::Align& al : rm->align_ ) if ( al.read_ ) joins.insert( al.tar_ );
        for ( ReadMatch::Align& al : rm->align_ ) if ( al.readEnd() < rm->size() ) joins.insert( al.tarEnd() );
    }
    for ( pair<size_t, size_t>& limit : limits ) joins.insert( limit.first );
    for ( pair<size_t, size_t>& limit : limits ) joins.insert( limit.second );
    
    vector<size_t> breaks( joins.begin(), joins.end() );
    sort( breaks.begin(), breaks.end() );
    for ( int i = 1; i < breaks.size(); i++ ) base.push_back( new Node( seq, base.empty() ? NULL : base.back(), NULL, breaks[i-1], breaks[i], -1 ) );
    
    int p = 0;
    vector<BundleRead*> reads;
    for ( ReadMatch* rm : matches ) reads.push_back( new BundleRead( rm, base, bubbles, p ) );
    for ( Node* node : bubbles ) node->setMismatches( seq );
    
    for ( Node* node : bubbles )
    {
        for ( int d : { 0, 1 } ) if ( node->type_ == d || node->type_ == 2 ) assert( !node->edges_[!d].empty() && node->edges_[!d][0]->edge_[!d]->type_ == -1 );
    }
    for ( int i = 1; i < base.size(); i++ ) assert( !base[i]->edges_[0].empty() && base[i]->edges_[0][0]->edge_[0]->type_ == -1 );
    for ( int i = 0; i+1 < base.size(); i++ ) assert( !base[i]->edges_[1].empty() && base[i]->edges_[1][0]->edge_[1]->type_ == -1 );
    
    return reads;
}

void Bundle::setPaths( vector<Node*>& base, vector<Node*>& bubbles, vector<BundleRead*>& reads )
{
    
}

