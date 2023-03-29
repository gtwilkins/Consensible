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

#include "consensus_kmers.h"
#include "read.h"
#include "shared_functions.h"

void ConsensusKmers::addMatch( Match* match, int range[2] )
{
    auto addKmer = [&]( unordered_map<uint16_t, Kmer>& kmers, uint16_t k, int i ){
        auto it = kmers.find( k );
        if ( it != kmers.end() ) it->second.coords_.push_back( make_pair( match, i ) );
        else kmers.insert( make_pair( k, Kmer( match, i ) ) );
    };
    vector<uint16_t> kmers = get8mers( match->read_->seq_, 0, match->read_->seq_.size() );
    for ( int i = 0; i < kmers.size(); i++ )
    {
        addKmer( kmers_, kmers[i], i );
        if ( i < range[0] ) addKmer( ends_[0], kmers[i], i );
        if ( range[1] < i+8 ) addKmer( ends_[1], kmers[i], i );
    }
}

Kmer* ConsensusKmers::getEndKmer( uint16_t k, bool drxn )
{
    auto it = ends_[drxn].find( k );
    return it != ends_[drxn].end() ? &it->second : NULL;
}

Kmer* ConsensusKmers::getKmer( uint16_t k )
{
    auto it = kmers_.find( k );
    return it != kmers_.end() ? &it->second : NULL;
}
