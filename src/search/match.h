/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   glin_match.h
 * Author: glen
 *
 * Created on 18 May 2022, 11:28 AM
 */

#ifndef GLIN_MATCH_H
#define GLIN_MATCH_H

#include "types.h"

class Target;
struct MappedRead;

struct Match
{
    Match( Target* tar, MappedRead* read, string t, string r, vector<pair<int,int>> &anchors, int tarCoord );
    int getAnchor( int drxn );
    vector<pair<int,int>> getGaps();
    Target* tar_;
    MappedRead* read_;
    vector<int> coords_;
    int anchors_[2], score_;
};


#endif /* GLIN_MATCH_H */

