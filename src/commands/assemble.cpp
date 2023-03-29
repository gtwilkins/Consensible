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

#include "assemble.h"
#include "filenames.h"
#include "match_query.h"
#include "parameters.h"
#include "timer.h"
#include "shared_functions.h"
#include <ctime> 
#include <iostream>
#include <string.h>
#include <unistd.h>
#include <cassert>
#include <algorithm>
#include "result.h"
#include "arguments.h"
#include <sys/stat.h>
#include <chrono>
#include <iomanip>

extern Parameters params;

Assemble::Assemble( Arguments& args )
:ir_( NULL ), qb_( NULL )
{
    string ifn, ofn, header, seq;
    int errors = 15;
    Filenames* fns = new Filenames( args.bwtPrefix_ );
    ir_ = new IndexReader( fns );
    qb_ = new QueryBinaries( fns );
    
    Result result;
    for ( string& ifn : args.queries_ )
    {
        ifstream ifs( ifn );
        while ( getSeq( ifs, header, seq ) )
        {
            Target* tar = result.addTarget( header, seq );
            for ( Read r : MatchQuery( seq, ir_, errors ).yield( qb_ ) ) result.addMatch( tar, r.id_, r.seq_, r.coords_[0] );
            break;
        }
        break;
    }
    result.assemble( args.outPrefix_ );
//    result.outputFullAlign( args.outPrefix_ );
    delete fns;
}

Assemble::~Assemble()
{
    delete ir_;
    delete qb_;
}

void Assemble::printUsage()
{
    cout << endl << "LeanBWT version " << LEANBWT_VERSION << endl;
    cout << endl << "Command:" << endl;
    cout << "    match" << endl;
    cout << endl << "Synopsis:" << endl;
    cout << "    Matches reads to one or more query sequences. By default, matches are exact. To find inexact matches, set an allowed mismatch rate (up to 15%)." << endl;
    cout << endl << "Usage:" << endl;
    cout << "    locass index [args]" << endl;
    cout << endl << "Required arguments:" << endl;
    cout << "    -p    Prefix for BWT data files." << endl;
    cout << endl << "Optional arguments:" << endl;
    cout << "    -o    Output file name (default: ./match_result.fa)." << endl;
    cout << "    -i    Input sequence query file (mutually exclusive with -s)." << endl;
    cout << "    -s    Input sequence query (mutually exclusive with -i)." << endl;
    cout << "    -e    Allowed mismatches per 100 bases for inexact matching (default: 0, maximum: 15)." << endl;
}
