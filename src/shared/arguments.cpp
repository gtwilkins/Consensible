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

#include "arguments.h"
#include <string.h>
#include <iostream>
#include <cassert>

using namespace std;

Arguments::Arguments( int argc, char** argv )
: bwtPrefix_( "./consensible-bwt" ), outPrefix_( "./consensible-out" ), reindex_( false ), cleanup_( false ), help_( false )
{
    for ( int i ( 1 ); i < argc; i++ )
    {
        if ( !strcmp( argv[i], "-i" ) )
        {
            if ( i+1 >= argc || argv[i+1][0] == '-') error( "No inputs given with -i flag" );
            while ( i+1 < argc && argv[i+1][0] != '-' ) inputs_.push_back( argv[++i] );
        }
        else if ( !strcmp( argv[i], "-p" ) )
        {
            if ( i+1 >= argc || argv[i+1][0] == '-') error( "No inputs given with -p flag" );
            bwtPrefix_ = argv[++i] ;
        }
        else if ( !strcmp( argv[i], "-o" ) )
        {
            if ( i+1 >= argc || argv[i+1][0] == '-') error( "No inputs given with -o flag" );
            outPrefix_ = argv[++i] ;
        }
        else if ( !strcmp( argv[i], "-q" ) )
        {
            if ( i+1 >= argc || argv[i+1][0] == '-') error( "No inputs given with -q flag" );
            while ( i+1 < argc && argv[i+1][0] != '-' ) queries_.push_back( argv[++i] );
        }
        else if ( !strcmp( argv[i], "--reindex" ) ) reindex_ = true;
        else if ( !strcmp( argv[i], "--cleanup" ) ) cleanup_ = true;
        else if ( !strcmp( argv[i], "-h" ) || !strcmp( argv[i], "--help" ) ) help_ = true;
        else error( "Unrecognised argument: \"" + string( argv[i] ) + "\"" );
    }
}

void Arguments::error( string msg )
{
    assert( false );
    cerr << msg << endl;
    exit( EXIT_FAILURE );
}