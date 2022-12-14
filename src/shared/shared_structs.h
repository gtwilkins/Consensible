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

#ifndef SHARED_STRUCTS_H
#define SHARED_STRUCTS_H

#include <vector>
#include "types.h"

struct Coords
{
    Coords():Coords( 0, 0 ){};
    Coords( int32_t i, int32_t j ){ coords_[0] = i; coords_[1] = j; };
    int32_t& operator[]( int i ){ return coords_[i]; };
    void operator+=( int i ){ coords_[0] += i; coords_[1] += i; };
    static vector< pair<Coords, Coords> > align( string a, string b, int minAlign );
    static bool conflict( pair<Coords, Coords>& a, pair<Coords, Coords>& b );
    bool conflict( Coords& alt, bool drxn );
    int len(){ return coords_[1] - coords_[0]; };
    int32_t coords_[2];
};

struct Mapped
{
    Mapped( ReadId id, Coords coords ): id_( id ), coords_( coords ){};
    static unordered_set<ReadId> getIds( vector<Mapped>& reads, int coord, bool coordDrxn, bool readDrxn );
    static void sort( vector<Mapped>& reads, bool ascending, bool coordDrxn );
    ReadId id_;
    Coords coords_;
};

struct Read
{
    Read( string seq, ReadId id, int i, int j ): seq_( seq ), id_( id ), coords_( Coords( i, j ) ){};
    static unordered_set<ReadId> getIds( vector<Read>& reads, int coord, bool coordDrxn, bool readDrxn );
    static void sort( vector<Read>& reads, bool ascending, bool coordDrxn );
    string seq_;
    ReadId id_;
    Coords coords_;
};


#endif /* SHARED_STRUCTS_H */

