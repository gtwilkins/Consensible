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

#ifndef CONSENSUS_H
#define CONSENSUS_H

#include "match.h"
#include "consensus_kmers.h"
#include "align_result.h"
#include "bubble.h"
#include "consensus_map.h"

class Target;

class Consensus
{
    struct Branch
    {
        string template_, consensus_;
        int coord_, drxn_;
    };
    void addBranch( ConMap* cm, vector<pair<ConMap*, AlignResult>>& hits, bool drxn );
    bool addBubble( Match* l, Match* r, AlignResult& result );
    void addMatch( Match* match );
    void addMatch( ConMap* cm, Match* m, int blockStart, int blockLast, bool preexisting );
    void addMismatch( ConMap* cm, Match* match, int lGood, int rGood );
    void addSNP( Match*, int start, int tarLen, string seq, int coord );
    bool foldEnd( ConMap* cm, vector<Bubble*>* branches, bool force, bool cleanup, bool drxn );
    void foldEnds();
    ConMap* getConMap( Match* m );
    void mapKmers();
    bool merge( Consensus* rhs );
    bool merge( Consensus* rhs, AlignResult& result, Match* l, Match* r );
    void resolveBranches();
    void resolveBubbles();
    void setBranches();
    void setBubbles();
    void updateMapped( ConMap* cm, SnpAlignResult& result, int base, bool drxn );
    Target* tar_;
    ConsensusKmers* kmers_;
    int coord_[2];
    string template_, consensus_;
    vector<ConMap*> maps_;
    vector<Bubble*> bubble_, branch_[2];
    vector<SNPs*> snps_;
public:
    Consensus( vector<Match*>& matches, Target* tar );
    ~Consensus();
    string resolve();
    static bool bridge( Consensus* lhs, Consensus* rhs );
    static void reSort( vector<Consensus*>& cs );
};



#endif /* GLIN_CONSENSUS_H */

