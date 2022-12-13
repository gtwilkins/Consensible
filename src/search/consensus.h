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
    Target* tar_;
    int coord_[2];
    string seq_;
    vector<ConMap> maps_;
public:
    Consensus( vector<Match*>& matches, Target* tar );
    static vector<Consensus*> seed( Match* match, int i, unordered_set<Match*>& used );
};



#endif /* GLIN_CONSENSUS_H */

