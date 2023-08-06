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

#ifndef GLIN_ALIGNMENT_H
#define GLIN_ALIGNMENT_H

#include "types.h"
#include "bubble.h"
#include "align_result.h"

struct SNPs;

class Alignment
{
    void score();
    void debug();
    bool isFreeStart( int i, int j );
//    vector< vector<int> > p_;
protected:
    string a_, b_;
    vector< vector<int> > m_, p_;
    int hit_, miss_, gapOpen_, gapExt_, iMax_, jMax_;
    bool freeEnds_[2], scored_;
public:
    Alignment( string& a, string& b, int aStart, int aLen, int bStart, int bLen );
    AlignResult align( int hit, int miss, int gapOpen, int gapExt, bool lAnchored, bool rAnchored );
    AlignResult align( bool lAnchored, bool rAnchored );
    static AlignResult alignBySeed( string& a, string& b, int aStart, int bStart, int len );
};

class SnpAlignment : Alignment
{
    struct GapPointer
    {
        void set( SnpAlignment* align, int i, int score )
        {
            align_ = align;
            i_ = i;
            score_ = score;
            len_ = del_ = 0;
        }
        int score( int gapOpen, int gapExt )
        {
            return score_ - ( len_ > del_ ? gapOpen + ( len_-del_ ) * gapExt : 0 );
        }
        SnpAlignment* align_;
        int i_, score_, len_, del_;
    };
    struct AlignPointer
    {
        bool set( SnpAlignment* align, int i, int j, int score )
        {
            align_ = align;
            i_ = i;
            j_ = j;
            score_ = score;
            return true;
        }
        bool update( SnpAlignment* align, int i, int j, int score )
        {
            if ( score > score_ ) return set( align, i, j, score );
            return false;
        }
        SnpAlignment* align_;
        int i_, j_, score_;
    };
    struct BubbleAlign
    {
        BubbleAlign( SnpAlignment* parent, Bubble* b, int coord, int parentLen  );
        BubbleAlign( SnpAlignment* parent, SNPs* snps, int i, int coord );
        bool isThis( SnpAlignResult::BubbleAlignCoords& bac );
        void set( vector<GapPointer> iMax, int ii );
        SnpAlignment* align_,* parent_;
        Bubble* bubble_;
        SNPs* snps_;
        SNP* snp_;
        vector<GapPointer> iMax_;
        int start_, end_;
    };
//    SnpAlignment( string& a, string& b, int aStart, int aLen, int bStart, int bLen );
    SnpAlignment( BubbleAlign* wrapper, string& a, int coord );
//    SnpAlignment( SnpAlignment* base, Bubble* bubble );
//    SnpAlignment( SnpAlignment* base, SNPs* snps, int i );
//    SnpAlignResult alignOld( bool lAnchored, bool rAnchored );
//    bool build( SnpAlignment& result, AlignPointer*& ap, int& i, int& j );
//    AlignPointer* build( SnpAlignResult& result, vector<BubbleAlign*> stack, int i, int& j );
    void add( SnpAlignResult& result, int i );
    void fill( SnpAlignResult& result, int& i, int iEnd, bool deletion, bool descending, bool ascending );
    vector<pair<SnpAlignment*, int>> findEnd();
    void finish( SnpAlignResult& result, AlignPointer& start );
    AlignPointer getEnd();
//    bool getStack( vector<BubbleAlign*>& stack, SnpAlignment* base );
    vector<BubbleAlign*> getStack( SnpAlignment* base );
//    void scoreOld();
    void score();
    void score( vector<GapPointer>& iMax );
    void setBubble( BubbleAlign* ba, vector<GapPointer> iMax, int ii );
    void getBubble( BubbleAlign* ba, vector<GapPointer>& iMax, int i, int j, bool& matched );
    vector<SNPs*> snps_;
    vector< vector<AlignPointer> > s_;
//    vector< vector< pair<int, int> > > snp_;
//    vector< pair<Bubble*, SnpAlignment> > bubble_;
    vector<BubbleAlign*> bubbles_;
    BubbleAlign* wrapper_;
    SnpAlignment* parent_;
    pair<int,int> ijMax_;
    int coord_[2];
public:
    SnpAlignment( string& a, string& b, int aStart, int aLen, int bStart, int bLen, vector<SNPs*> snps, vector<Bubble*>& bubbles );
    SnpAlignResult align( bool lAnchored, bool rAnchored );
};


#endif /* GLIN_ALIGNMENT_H */

