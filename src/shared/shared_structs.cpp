/*
 * Copyright (C) 2017 Glen T. Wilkins <glen.t.wilkins@gmail.com>
 * Written by Glen T. Wilkins
 * 
 * This file is part of the LeanBWT software package <https://github.com/gtwilkins/LeanBWT>
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

#include "shared_structs.h"
#include "shared_functions.h"
#include "local_alignment.h"
#include "parameters.h"
#include <algorithm>
#include <cassert>
#include <iostream>
#include <sys/stat.h>
#include <chrono>
#include <iomanip>

vector< pair<Coords, Coords> > Coords::align( string a, string b, int minAlign )
{
    int len = minAlign / 2;
    vector< pair<Coords, Coords> > aligns;
    if ( len < a.size() ) for ( int i = 0; i < a.size(); )
    {
        int j = min( i, (int)a.size()-len );
        string q = a.substr( j, len );
        size_t it = b.find( q );
        while ( it != string::npos )
        {
            Coords x[2]{ Coords( j, j+len ), Coords( it, it+len ) };
            while ( x[0][0] && x[1][0] && a[ x[0][0]-1 ] == b[ x[1][0]-1 ] ){ x[0][0]--; x[1][0]--; };
            while ( x[0][1] < a.size() && x[1][1] < b.size() && a[ x[0][1] ] == b[ x[1][1] ] ){ x[0][1]++; x[1][1]++; };
            bool good = x[0][1]-x[0][0] >= minAlign;
            for ( pair<Coords, Coords>& align : aligns ) if ( align.first[0] == x[0][0] && align.first[1] == x[0][1] && align.second[0] == x[1][0] ) good = false;
            if ( good ) aligns.push_back( make_pair( x[0], x[1] ) );
            it = b.find( q, it+1 );
        }
        i += len;
    }
    
    return aligns;
}

bool Coords::conflict( pair<Coords, Coords>& a, pair<Coords, Coords>& b )
{
    return b.first[0] < a.first[1] || b.second[0] < a.second[1];
}

unordered_set<ReadId> Mapped::getIds( vector<Mapped>& reads, int coord, bool coordDrxn, bool readDrxn )
{
    unordered_set<ReadId> ids;
    for ( Mapped& read : reads ) if ( readDrxn ? coord < read.coords_[coordDrxn] : read.coords_[coordDrxn] < coord ) ids.insert( read.id_ );
    return ids;
}

void Mapped::sort( vector<Mapped>& reads, bool ascending, bool coordDrxn )
{
    if ( ascending ) std::sort( reads.begin(), reads.end(), [&]( Mapped& a, Mapped& b ){ return a.coords_[coordDrxn] < b.coords_[coordDrxn]; } );
    else std::sort( reads.begin(), reads.end(), [&]( Mapped& a, Mapped& b ){ return a.coords_[coordDrxn] > b.coords_[coordDrxn]; } );
}

unordered_set<ReadId> Read::getIds( vector<Read>& reads, int coord, bool coordDrxn, bool readDrxn )
{
    unordered_set<ReadId> ids;
    for ( Read& read : reads ) if ( readDrxn ? coord < read.coords_[coordDrxn] : read.coords_[coordDrxn] < coord ) ids.insert( read.id_ );
    return ids;
}

void Read::sort( vector<Read>& reads, bool ascending, bool coordDrxn )
{
    if ( ascending ) std::sort( reads.begin(), reads.end(), [&]( Read& a, Read& b ){ return a.coords_[coordDrxn] < b.coords_[coordDrxn]; } );
    else std::sort( reads.begin(), reads.end(), [&]( Read& a, Read& b ){ return a.coords_[coordDrxn] > b.coords_[coordDrxn]; } );
}
