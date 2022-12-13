/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   glin_alignment.h
 * Author: glen
 *
 * Created on 24 October 2022, 4:06 PM
 */

#ifndef GLIN_ALIGNMENT_H
#define GLIN_ALIGNMENT_H

#include "types.h"

class Alignment
{
    void debug();
    bool isFreeStart( int i, int j );
    void score();
    string a_, b_;
    vector< std::vector<int> > m_, p_;
    int hit_, miss_, gapOpen_, gapExt_, iMax_, jMax_;
    bool freeEnds_[2], scored_;
public:
    Alignment( string& a, string& b, int aStart, int aLen, int bStart, int bLen );
    pair<int,int> align( bool lAnchored, bool rAnchored );
};


#endif /* GLIN_ALIGNMENT_H */

