/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   glin_consensus.h
 * Author: glen
 *
 * Created on 16 September 2022, 3:11 PM
 */

#ifndef GLIN_CONSENSUS_H
#define GLIN_CONSENSUS_H

#include "match.h"

class Target;


class Consensus
{
    struct ConMap
    {
        Match* node_;
        int coord_[2]/*of consensus*/, range_[2]/*of read*/;
    };
    Target* tar_;
    int coord_[2];
    string seq_;
    vector<ConMap> maps_;
public:
    Consensus( vector<Match*>& matches, Target* tar );
    static vector<Consensus*> seed( Match* match, int i, unordered_set<Match*>& used );
};



#endif /* GLIN_CONSENSUS_H */

