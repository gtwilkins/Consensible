/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   glin_result.h
 * Author: glen
 *
 * Created on 18 May 2022, 11:25 AM
 */

#ifndef RESULT_H
#define RESULT_H

#include "types.h"
#include "target.h"
#include "read.h"

class Result
{
    MappedRead* addRead( ReadId id, string seq );
    
    vector<Target*> targets_;
    unordered_map<ReadId, MappedRead*> reads_;
public:
    ~Result();
    Target* addTarget( string header, string seq );
    void addMatch( Target* tar, ReadId id, string seq, int coord );
};



#endif /* GLIN_RESULT_H */

