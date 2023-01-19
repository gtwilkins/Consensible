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
    void assemble( string& outPrefix );
    void outputFullAlign( string& outPrefix );
};



#endif /* GLIN_RESULT_H */

