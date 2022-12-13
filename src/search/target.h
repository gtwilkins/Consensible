/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   glin_target.h
 * Author: glen
 *
 * Created on 18 May 2022, 11:27 AM
 */

#ifndef GLIN_TARGET_H
#define GLIN_TARGET_H

#include "types.h"
#include "match.h"
#include "read.h"

class Target
{
    vector<Match*> matches_;
    
    vector<pair<int,int>> getGaps();
    void sortMatches();
public:
    Target( string header, string seq ): header_( header ), seq_( seq ){};
    bool addMatch( MappedRead* read, int coord );
    void assemble();
    void print( ofstream& ofs );
    string seq_, header_;
};

#endif /* GLIN_TARGET_H */

