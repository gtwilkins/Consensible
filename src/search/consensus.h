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
    struct Bubble
    {
        string template_, consensus_;
        int start_, len_;
    };
    struct Branch
    {
        string template_, consensus_;
        int coord_, drxn_;
    };
    struct SNPs
    {
        struct SNP
        {
            vector< pair<Match*, int> > matches_;
            string seq_;
        };
        string resolve( string base );
        vector<SNP> snps_;
        int start_, len_;
    };
    void addBranch( Match* match, int len, bool drxn );
    void addMatch( Match* match );
    void addMismatch( Match* match, int lGood, int rGood );
    void resolveBranches();
    void resolveBubbles();
    Target* tar_;
    int coord_[2];
    string template_, consensus_;
    vector<ConMap> maps_;
    vector<Bubble> bubble_;
    vector<Branch> branch_[2];
    vector<SNPs> snps_;
public:
    Consensus( vector<Match*>& matches, Target* tar );
    void resolve();
};



#endif /* GLIN_CONSENSUS_H */

