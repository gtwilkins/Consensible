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

#ifndef BUNDLE_H
#define BUNDLE_H

#include "types.h"
#include "long_query.h"


class Bundle
{
//    struct BundleEdge;
//    struct BundleNode
//    {
//        BundleNode( ReadMatch* rm ): node_( rm ){};
//        ReadMatch* node_;
//        vector<BundleEdge*> edges_;
//    };
//    
//    struct BundleEdge
//    {
//        struct Align
//        {
//            Align(){};
//            Align( size_t lStart, size_t rStart, size_t len ): len_( len ){ s_[0] = lStart; s_[1] = rStart; };
//            size_t s_[2], len_;
//        };
//        BundleEdge( BundleNode* l, BundleNode* r );
//        vector<Align> align_;
//        int score_, errs_, start_[2], end_[2];
//    };
    struct BundleNode;
    
    struct BBubble
    {
        BBubble( string& seq, ReadMatch* rm, size_t start, size_t end, size_t drxn ): seq_( seq ), matches_{ rm }, start_( start ), end_( end ), drxn_( drxn ){};
        static vector<BBubble*> getBubbles( vector<BundleNode*>& nodes );
        bool setMatch( string& seq, ReadMatch* rm, size_t start, size_t end, size_t drxn );
        bool setNonmatch( ReadMatch* rm );
        string seq_;
        vector<ReadMatch*> matches_, nonmatches_;
        size_t start_, end_, drxn_;
    };
    
    struct BundleNode
    {
        BundleNode( ReadMatch* rm );
        ReadMatch* node_;
        vector<BBubble*> bubs_;
        size_t tar_[2];
        int score_;
    };
    
    struct PPath
    {
        struct Branch
        {
            bool add( BundleNode* bn );
            vector<BundleNode*> nodes_;
            vector< pair<BBubble*, bool> > overlaps_, extends_;
            int disagrees_;
        };
        PPath( BundleNode* base, vector<BundleNode*>& nodes, vector<BBubble*>& bubs );
        bool extend( vector<BundleNode*>& nodes, vector<BBubble*>& bubs );
        bool extend( vector<BundleNode*>& nodes, vector<BBubble*>& bubs, bool drxn );
        ReadMatch* base_;
        vector<BBubble*> hit_, miss_;
        size_t tar_[2];
    };
    struct BundleRead;
    struct Edge;
    
    struct Node
    {
        Node( string& seq, Node* lEdge, BundleRead* read, size_t lCoord, size_t rCoord, int type );
        ~Node();
        bool addBubble( string& seq, BundleRead* read, size_t lCoord, size_t rCoord, int type );
//        Node( string seq, BundleRead* rm, size_t* coord, int type ): seq_( seq ), reads_( 1, rm ), type_( type ) { coord_[0] = coord[0], coord_[1] = coord[1]; };
        static Node* setBubble( vector<Node*>& bubbles, Node* lEdge, BundleRead* read, string seq, size_t lCoord, size_t rCoord, int type );
       void setMismatches( string& base );
//        static void setNodes( string& seq, vector<Node*>& base, vector<Node*>& bubbles, vector<ReadMatch*>& reads );
        string seq_;
        vector<BundleRead*> reads_;
        vector<Edge*> edges_[2];
        size_t coord_[2];
        int mismatches_, type_;
    };
    
    struct Edge
    {
        static bool thread( Node* l, Node* r, BundleRead* read );
        vector<BundleRead*> reads_;
        Node* edge_[2];
    private:
        Edge( Node* l, Node* r );
    };
    
    struct BundleRead
    {
//        BundleRead( ReadMatch* rm, vector<Node*>& bubs );
        BundleRead( ReadMatch* rm, vector<Node*>& base, vector<Node*>& bubs, int& p );
        int getMismatches();
        vector<Node*> getPath( vector<Node*>& base );
        int getRemainder( Edge* edge, bool drxn );
        bool isGood( size_t errorRate );
        bool isPerfect();
        static void sort( vector<BundleRead*>& reads );
        ReadMatch* read_;
        vector<Node*> nodes_;
        size_t tar_[2];
        int score_;
    };
    
    struct Path
    {
        struct Extend
        {
            Extend( BundleRead* read, Edge* edge, bool drxn );
            Extend( Extend* ext, unordered_set<BundleRead*> reads );
//            static void init( vector<Node*>& nodes, vector<Extend*>& ext, vector<Extend*>& splinters );
            bool shrink( Node* node, bool drxn );
            bool splinter( Edge* edge, int dist, bool drxn );
            BundleRead* read_;
            vector<Node*> ext_[2];
            unordered_set<BundleRead*> reads_;
            int ol_[2], len_[3]; // [0] = left ext, [1] = right ext, [2] = aligned to path
        };
        struct Branch
        {
            Edge* edge_;
            Node* node_;
            vector<Extend*> extends_;
        };
        Path( vector<Node*> path, BundleRead* br );
        void addReads();
        void extend();
        bool extend( Branch& branch, vector<Branch> (&branches)[2], vector<Extend*>& extends, unordered_set<BundleRead*>& miss, bool drxn );
        void evaluate( Branch& branch, bool drxn );
        void initExtends( vector<Extend*>& exts, vector<Extend*>& splinters, bool drxn );
        vector<Extend*> getExtends( unordered_set<BundleRead*>& miss );
        vector<BundleRead*> getExtending( unordered_set<BundleRead*>& ignore, bool drxn );
        unordered_set<BundleRead*> getTangents( bool inclEnds=false );
        bool resolve( vector<Branch>& fwd, vector<Branch>& rev, vector<Extend*>& extends, bool drxn );
        vector<Branch> setBranches( vector<Extend*>& extends, bool drxn );
        void setExtend( Extend* ext, Edge* edge, vector<Extend*>& extends, bool drxn );
        void setOverlaps( vector<Extend*>& extends );
        vector<Node*> nodes_;
        vector<BundleRead*> reads_;
        size_t coords_[2];
    };
    
    Bundle( BundleNode* bn, vector<BundleNode*>& nodes, vector<BBubble*>& bubs );
    static void add( ReadMatch* rm, vector<BBubble*>& bubs );
    static void build( BundleRead* read, vector<Node*>& base, unordered_set<BundleRead*>& used );
    static void getBase( vector<Bundle*>& bundles, vector<Node*>& base, vector<Node*>& bubbles, vector<BundleRead*>& reads );
    static vector<BundleRead*> setGraph( string& seq, vector<Node*>& base, vector<Node*>& bubbles, vector<ReadMatch*>& matches );
    static void setPaths( vector<Node*>& base, vector<Node*>& bubbles, vector<BundleRead*>& reads );
    static void setScores( vector<Node*>& bubbles );
public:
    static vector<Bundle*> bundle( string& seq, LongResult* lr );
    vector<ReadMatch*> matches_;
};

#endif /* BUNDLE_H */

