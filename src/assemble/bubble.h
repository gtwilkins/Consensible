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


#ifndef BUBBLE_H
#define BUBBLE_H

#include "types.h"
#include "align_result.h"
#include "consensus_map.h"

struct Match;

struct SNP
{
    vector< pair<Match*, int> > matches_;
    string seq_;
};

struct SNPs
{
    SNPs():score_( 0 ){};
    static void addSnp( vector<SNPs*>& snps, string seq, ConMap* cm, int baseCoord, int matchCoord, int len );
    void addSnp( string seq, Match* m, int coord );
    void addMatch( SnpAlignResult& result, SnpAlignResult::BubbleAlignCoords& bac, ConMap* cm, int& i, int& read );
    string getConsensus();
    bool isFoldTarget( int (&limits)[2] );
    string resolve( string base );
    int score( unordered_set<Match*>& matches );
    static void setRemapped( vector<pair<int,int>>& remapped, unordered_map<SNPs*,int>& oldSnps, vector<SNPs*>& curSnp );
    vector<SNP> snps_;
    int start_, len_, score_;
};

struct Bubble
{
//    struct IntraBubble
//    {
//        Bubble* bubble_;
//        int start_, end_;
//    };
    struct BubbleCoord
    {
        BubbleCoord(){};
        BubbleCoord( int bubble, int match, int len ): bubble_( bubble ), match_( match ), len_( len ){ };
        int bubble_, match_, len_;
    };
    struct BubbleMap
    {
        BubbleMap( ConMap* cm ): map_( cm ){};
        BubbleMap( Match* match, vector< pair<int,int> >& coords, int lanchor, int ranchor );
        void addCoord( int match, int bubble );
        ConMap* map_;
        Match* match_;
        vector<BubbleCoord> coords_;
//        vector<IntraBubble*> bubbles_;
        int anchors_[2];
    };
    struct BubbleMapDetails
    {
        BubbleMapDetails( Bubble* bub, BubbleMap& map );
        int setEnd( Bubble* bub, BubbleMap& map );
        ConMap* map_;
        string seq_;
        int type_, start_, end_, off_;
    };
    struct BubbleAltPath
    {
        BubbleAltPath(){};
        BubbleAltPath( BubbleMapDetails& map ): seq_( map.seq_ ), map_( 1, map ){};
        bool add( BubbleMapDetails& map );
        void assign( int* coords );
        static vector<BubbleAltPath> build( Bubble* b );
        int getOffset( BubbleMapDetails& map );
        vector<Bubble*> getStack( int d );
        void remap( SnpAlignResult::BubbleAlignCoords& bac, int& coord );
        bool set( SnpAlignResult result );
        vector<Bubble*> split( Bubble* base, SnpAlignResult::BubbleAlignCoords& bac );
        string seq_;
        SnpAlignResult result_;
        SnpAlignResult::BubbleAlignCoords bestAlign_;
        vector<BubbleMapDetails> map_;
        vector< pair<Bubble*,pair<int,int>> > path_;
        int cut_[2], coord_[2], score_;
    };
    Bubble(): type_( 2 ), score_( 0 ){};
    Bubble( Bubble* base, BubbleAltPath& path, SnpAlignResult::BubbleAlignCoords& bac, int cut, bool drxn );
    Bubble( ConMap* cm, bool drxn );
    Bubble( AlignResult& result, ConMap* a, ConMap* b, bool drxn );
    Bubble( AlignResult& result, ConMap* l, ConMap* r );
    ~Bubble();
    void abort();
    static void addBubble( vector<Bubble*>& bubbles, string seq, ConMap* cm, int baseCoord, int matchCoord, int len );
    void addBubble( pair<int,int> bubble, pair<int,int> read, ConMap* cm );
    static Bubble* addBubble( string seq, vector<Bubble*>& bubbles, int start, int len );
    void addMatch( ConMap* cm, int start );
    void addMatch( AlignResult& result, ConMap* a, ConMap* b, bool drxn );
    void addMatch( SnpAlignResult& result, SnpAlignResult::BubbleAlignCoords& bac, ConMap* cm, int& i, int& read );
    bool consolidate( vector<SNPs*>& snps, vector<Bubble*>& bubs, vector<Bubble*>& baseBubs, string& baseSeq );
    void cut( vector<Bubble*> cutStack, vector<Bubble*> keepStack, int cut, bool drxn );
//    vector<Bubble::BubbleAltPath> getAlignableSeqs();
    string getConsensus();
    void getMapRange( ConMap* cm, int& start, int& end );
//    void getSeq( BubbleMap& map );
//    vector<pair<string,vector<Bubble*>>> getAlignableSeqs();
    bool isClash( Bubble* b );
    bool isFoldTarget( ConMap* cm, int (&limits)[2], bool finalised, bool drxn );
    bool remove( ConMap* cm );
    int score( unordered_set<ConMap*>& maps );
    static void setRemapped( vector<pair<int,int>>& remapped, unordered_map<Bubble*,int>& oldBubbles, vector<Bubble*>& curBubbles );
    static void test( vector<Bubble*>& bubs );
    vector<BubbleMap> maps_;
//    vector<IntraBubble*> bubbles_;
    vector<Bubble*> bubs_;
    string template_, consensus_;
    int start_, len_, type_, score_;
};


#endif /* BUBBLE_H */

