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

#ifndef CONSENSUS_RESOLUTION_H
#define CONSENSUS_RESOLUTION_H

#include "bubble.h"
#include "consensus_map.h"
#include "types.h"

struct ConsensusResolution
{
    ConsensusResolution( Bubble* b ): bub_( b ), snps_( NULL ), start_( b->start_), end_( b->start_+b->len_ ), score_( b->score_ ){};
    ConsensusResolution( SNPs* snps ): bub_( NULL ), snps_( snps ), start_( snps->start_), end_( snps->start_+snps->len_ ), score_( snps->score_ ){};
    static void branch( vector<Bubble*>& branches, vector<ConMap*>& map, vector<ConsensusResolution>& resolves, bool drxn );
    static vector<ConsensusResolution> create( vector<Bubble*>& bubs, vector<SNPs*>& snps, vector<ConMap*>& maps );
    string getConsensus();
    bool isClash( ConsensusResolution& rhs );
    static vector<ConsensusResolution> resolve( vector<Bubble*> branches[2], vector<Bubble*>& bubs, vector<SNPs*>& snps, vector<ConMap*>& maps );
    Bubble* bub_;
    SNPs* snps_;
    int start_, end_, score_;
};

#endif /* CONSENSUS_RESOLUTION_H */

