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

#include "index.h"
#include "transform_structs.h"
#include <iostream>
#include <sys/stat.h>
#include <string.h>
#include <cassert>
#include <algorithm>
#include <fstream>
#include <unistd.h>

Index::Index( Arguments& args )
{
    PreprocessFiles* fns = new PreprocessFiles( args.bwtPrefix_, args.reindex_ );
    bool isComplete = false, canResume = false, isIndexed = false, doRevComp = true;
    int minScore = 0;
    fns->getState( isComplete, canResume, isIndexed );
    if ( args.reindex_ || ( !canResume && !isComplete ) )
    {
        Transform::load( fns, args.inputs_, doRevComp );
        Transform::run( fns );
        IndexWriter idx( fns, 1024, 20000 );
    }
    else if ( canResume  )
    {
        Transform::run( fns );
        IndexWriter idx( fns, 1024, 20000 );
    }
    else if ( !isIndexed  )
    {
        IndexWriter idx( fns, 1024, 20000 );
    }
    else
    {
        cout << "Index files already exist. Proceeding..." << endl << endl;
    }
    
    delete fns;
    
//    for ( int i ( 2 ); i < argc; i++ )
//    {
//        if ( !strcmp( argv[i], "-h" ) )
//        {
//            printUsage();
//            exit( EXIT_SUCCESS );
//        }
//        else if ( !strcmp( argv[i], "-i" ) )
//        {
//            if ( didInput )
//            {
//                cerr << "Error: multiple inputs provided." << endl;
//                exit( EXIT_FAILURE );
//            }
//            didInput = true;
//            infile.open( argv[++i] );
//            if ( !infile.good() )
//            {
//                cerr << "Failed to open file: \"" << argv[i] << "\"" << endl;
//                exit( EXIT_FAILURE );
//            }
//        }
//        else if ( !strcmp( argv[i], "-p" ) )
//        {
//            if ( fns )
//            {
//                cerr << "Error: more than one output prefix provided." << endl;
//                exit( EXIT_FAILURE );
//            }
//            prefix = argv[++i];
//            if ( prefix[0] != '/' )
//            {
//                string curr = getcwd( NULL, 0 );
//                if ( !prefix.empty() && prefix[0] == '.' ) prefix = prefix.substr( 1 );
//                if ( !prefix.empty() && prefix[0] != '/' ) prefix = "/" + prefix;
//                if ( prefix.empty() || curr.empty() )
//                {
//                    cerr << "Invalid file prefix. Please use the absolute path." << endl;
//                }
//                prefix = curr + prefix;
//            }
//        }
//        else if ( !strcmp( argv[i], "-s" ) ) minScore = stoi( argv[++i] );
//        else if ( !strcmp( argv[i], "--resume" ) ) isResume = true;
//        else if ( !strcmp( argv[i], "--no-rev-comp" ) ) doRevComp = false;
//        else
//        {
//            cerr << "Unrecognised argument: \"" << argv[i] << "\"" << endl << endl;
//            printUsage();
//            exit( EXIT_FAILURE );
//        }
//    }
    
//    if ( argc <= 2 )
//    {
//        cerr << "Error: no arguements supplied, see usage:" << endl << endl;
//        printUsage();
//        exit( EXIT_FAILURE );
//    }
//    
//    double preprocessStartTime = clock();
//    
//    if ( prefix.empty() )
//    {
//        cerr << "Error: no output prefix supplied." << endl;
//        exit( EXIT_FAILURE );
//    }
//    
//    fns = new PreprocessFiles( prefix, true );
//    
//    if ( isResume && didInput )
//    {
//        cerr << "Error: resume (--resume) and input (-i) are mutually exclusive arguments." << endl;
//        cerr << "If you wish to resume (--resume), simply specify the output prefix (-p) previously used." << endl;
//        exit( EXIT_FAILURE );
//    }
//    else if ( didInput )
//    {
//        newTransform( fns, minScore, infile, doRevComp );
//    }
//    else if ( isResume )
//    {
//        resumeTransform( fns );
//    }
//    else
//    {
//        cerr << "Error: did not specify either an input (-i), or the resume flag (--resume)." << endl << endl;
//        printUsage();
//        exit( EXIT_FAILURE );
//    }
//    
//    cout << "Preprocessing step 3 of 3: indexing transformed data..." << endl;
//    IndexWriter idx( fns, 1024, 20000 );
//    
//    cout << endl << "Preprocessing completed!" << endl;
//    cout << "Total time taken: " << getDuration( preprocessStartTime ) << endl;
}

//void Index::newTransform( PreprocessFiles* fns, vector<string>& infilenames, bool revComp )
//{
//    if ( infilenames.size() > 6 )
//    {
//        cerr << "Error: Excessive input read file count of " << (int)infilenames.size() << ". Maximum supported is 6." << endl;
//        exit( EXIT_FAILURE );
//    }
//    int minScore = 0;
//    vector<ReadFile*> infiles;
//    for ( string& ifn : infilenames ) infiles.push_back( new ReadFile( ifn, 0, minScore ) );
//    Transform::load( fns, infiles, revComp );
//    Transform::run( fns );
//}

//void Index::newTransform( PreprocessFiles* fns, int minScore, ifstream &infile, bool revComp )
//{
//    uint8_t fileCount = 0, pairedLibCount = 0;
//    
//    if ( infile.is_open() && infile.good() )
//    {
//        vector< vector<ReadFile*> > libs;
//        string line;
//        int lineNum = 1;
//        ReadFile* readFile = NULL;
//        while( getline( infile, line ) )
//        {
//	    if ( line.empty() || line[0] == '#' ) continue;
//            vector<string> args;
//            
//            size_t it = line.find( ' ' );
//            while ( it != line.npos && args.size() < 4 )
//            {
//                args.push_back( line.substr( 0, it ) );
//                line = line.substr( it + 1 );
//                it = line.find( ' ' );
//            }
//            
//            args.push_back( line );
//            
//            if ( ( args.size() == 2 || args.size() == 3 ) && args[0] == "paired" )
//            {
//                vector<ReadFile*> lib;
//                readFile = new ReadFile( args[1], 0, minScore );
//                fileCount++;
//                lib.push_back( readFile );
//                if ( args.size() == 3 )
//                {
//                    readFile = new ReadFile( args[2], 0, minScore );
//                    fileCount++;
//                }
//                lib.push_back( readFile );
//                libs.push_back( lib );
//                pairedLibCount++;
//            }
//            else if ( args.size() == 2 && args[0] == "single" )
//            {
//                readFile = new ReadFile( args[1], 0, minScore );
//                fileCount++;
//                vector<ReadFile*> lib = { readFile };
//                libs.push_back( lib );
//            }
//            else
//            {
//                cerr << "Error: invalid input form in line " << lineNum << " of input file." << endl << endl;
//                cerr << "Expecting inputs to be in one of the following forms:" << endl;
//                cerr << "\tpaired [separated_pair_file_1] [separated_pair_file_2]" << endl;
//                cerr << "\tpaired [interleaved_pair_file]" << endl;
//                cerr << "\tsingle [unpaired_file]" << endl;
//                exit( EXIT_FAILURE );
//            }
//            
//            lineNum++;
//        }
//        
//        if ( pairedLibCount > 6 )
//        {
//            cerr << "Error: Excessive library count of " << (int)pairedLibCount << ". Maximum supported is 6." << endl;
//            exit( EXIT_FAILURE );
//        }
//        
//        cout << "Preprocessing step 1 of 3: reading input files..." << endl << endl;
//        Transform::load( fns, libs, pairedLibCount, revComp );
//        Transform::run( fns );
//    }
//    else
//    {
//        cerr << "Error: No input file provided." << endl;
//        exit( EXIT_FAILURE );
//    }
//}
//
//void Index::resumeTransform( PreprocessFiles* fns )
//{
//    cout << "Resuming preprocessing..." << endl << endl;
//    Transform::run( fns );
//}

void Index::printUsage()
{
    
}
