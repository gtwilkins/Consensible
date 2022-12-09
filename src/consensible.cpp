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

#include <cstdlib>
#include <cassert>
#include <string.h>
#include <iostream>
#include "types.h"
#include "index.h"
#include "shared_functions.h"
#include "parameters.h"

Parameters params;

void printUsage()
{
    cout << endl << "Consensible version " << LEANBWT_VERSION << endl;
    cout << endl << "Usage:" << endl;
    cout << "\tconsensible [args]" << endl;
    cout << endl << "Arguments:" << endl;
    cout << "\t-i\t(Optional) Input sequence file(s)." << endl;
    cout << "\t-p\t(Required) Prefix for indexed sequence files. See notes for details." << endl;
    cout << "\t-q\t(Required) Query file containing one or more query sequences." << endl;
    cout << "\t-o\t(Optional) Output filename prefix." << endl;
    cout << endl << "Notes:" << endl;
    cout << "\t- Accepted read file formats are fasta, fastq or a list of sequences, one per line." << endl;
    cout << "\t- Input read libraries can be either paired or single." << endl;
    cout << "\t- Each paired read library can be input as either two separated files or one interleaved file." << endl;
    cout << "\t- Each line of the input text file is expected in one of the following forms:" << endl;
    cout << "\t\tpaired [separated_pair_file_1] [separated_pair_file_2]" << endl;
    cout << "\t\tpaired [interleaved_pair_file]" << endl;
    cout << "\t\tsingle [unpaired_file]" << endl;
    cout << "\t- At least one paired read library is required." << endl;
}

int main( int argc, char** argv )
{
    
    if ( argc > 1 )
    {
        if ( !strcmp( argv[1], "-h" ) || !strcmp( argv[1], "--help" ) || !strcmp( argv[1], "-help" ) )
        {
            printUsage();
        }
        else if ( !strcmp( argv[1], "index" ) )
        {
            Index( argc, argv );
        }
//        else if ( !strcmp( argv[1], "match" ) )
//        {
//            Match( argc, argv );
//        }
//        else if ( !strcmp( argv[1], "overlap" ) )
//        {
//            Overlap( argc, argv );
//        }
//        else if ( !strcmp( argv[1], "coverage" ) )
//        {
//            Coverage( argc, argv );
//        }
//        else if ( !strcmp( argv[1], "test" ) )
//        {
//            Test( argc, argv );
//        }
        else
        {
            cerr << "Unrecognised command: \"" << argv[1] << "\"" << endl << endl;
            printUsage();
        }
    }
    else
    {
        printUsage();
    }
    return 0;
}

