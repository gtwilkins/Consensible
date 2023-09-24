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
#include "arguments.h"
#include "assemble.h"

Parameters params;

void printUsage()
{
    cout << endl << "Consensible version " << LEANBWT_VERSION << endl;
    cout << endl << "Usage:" << endl;
    cout << "\tconsensible [args]" << endl;
    cout << endl << "Arguments:" << endl;
    cout << "\t-i\t(Required) Input shotgun sequence file(s)." << endl;
    cout << "\t-q\t(Required) Query file containing one or more query sequences." << endl;
    cout << "\t-w\t(Required) Working directory for temporary index files." << endl;
    cout << "\t-o\t(Optional) Output directory for results." << endl;
    cout << endl << "Example command:" << endl;
    cout << "\tconsensible -i /myinputs/project101_data.fa -q /myinputs/interesting_gene.fa -w /mytempdata -o /myoutput" << endl;
    cout << endl << "Explanation:" << endl;
    cout << "\tConsensible searches for reads among shotgun sequencing data (-i) that match" << endl;
    cout << "\tone or more query sequences (-q). It will then attempt to assmble overlapping" << endl;
    cout << "\treads among the matches into one or more consensus sequences that are output" << endl;
    cout << "\tas a fasta file (-o). However, consensible must first create an index of the" << endl;
    cout << "\tshotgun sequencing data (-w)." << endl;
    cout << endl << "Notes:" << endl;
    cout << "\t- Accepted read file formats are fasta, fastq or a list of sequences, one per line." << endl;
    cout << "\t- Once shotgun sequencing data (-i) has been indexed, only the index prefix (-p) need be provided on subsequent queries (-q)." << endl;
}

int main( int argc, char** argv )
{
    
    if ( argc > 1 )
    {
        Arguments arguments( argc, argv );
        if ( arguments.help_)
        {
            printUsage();
        }
        else
        {
            int i = 0;
            while ( arguments.setBwtPrefix() )
            {
                Index idx( arguments );
                arguments.updateFileIndex();
                Assemble ass( arguments );
                i++;
            }
        }
//        if ( !strcmp( argv[1], "-h" ) || !strcmp( argv[1], "--help" ) || !strcmp( argv[1], "-help" ) )
//        {
//            printUsage();
//        }
//        else if ( !strcmp( argv[1], "index" ) )
//        {
//            Index( argc, argv );
//        }
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
//        else
//        {
//            cerr << "Unrecognised command: \"" << argv[1] << "\"" << endl << endl;
//            printUsage();
//        }
    }
    else
    {
        printUsage();
    }
    return 0;
}

