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
#include <sys/stat.h>
#include <chrono>
#include <iomanip>

extern Parameters params;

Assemble::Assemble( int argc, char** argv )
:ir_( NULL ), qb_( NULL )
{
    string ifn, ofn, header, seq;
    int errors = 0;
    bool collapse = false, mismatches = false;
    Filenames* fns = NULL;
    
    for ( int i ( 2 ); i < argc; i++ )
    {
        if ( !strcmp( argv[i], "-h" ) )
        {
            printUsage();
            exit( EXIT_SUCCESS );
        }
        else if ( !strcmp( argv[i], "-i" ) )
        {
            if ( !ifn.empty() || !seq.empty() )
            {
                cerr << "Error: multiple inputs provided." << endl;
                exit( EXIT_FAILURE );
            }
            ifn = argv[++i];
        }
        else if ( !strcmp( argv[i], "-o" ) )
        {
            if ( !ofn.empty() )
            {
                cerr << "Error: multiple outputs provided." << endl;
                exit( EXIT_FAILURE );
            }
            ofn = argv[++i];
        }
        else if ( !strcmp( argv[i], "-p" ) )
        {
            if ( fns )
            {
                cerr << "Error: more than one output prefix provided." << endl;
                exit( EXIT_FAILURE );
            }
            string prefix = argv[++i];
            if ( prefix[0] != '/' )
            {
                string curr = getcwd( NULL, 0 );
                if ( !prefix.empty() && prefix[0] == '.' ) prefix = prefix.substr( 1 );
                if ( !prefix.empty() && prefix[0] != '/' ) prefix = "/" + prefix;
                if ( prefix.empty() || curr.empty() )
                {
                    cerr << "Invalid file prefix. Please use the absolute path." << endl;
                }
                prefix = curr + prefix;
            }
            
            fns = new Filenames( prefix );
        }
        else if ( !strcmp( argv[i], "-e" ) )
        {
            errors = stoi( argv[++i] );
            if ( errors > 15 || errors < 0 )
            {
                cerr << "Error: invalid mismatch rate set at: " << errors << "%, must be between 0-15%" << endl;
                exit( EXIT_FAILURE );
            }
        }
        else if ( !strcmp( argv[i], "-s" ) ) seq = argv[++i];
        else if ( !strcmp( argv[i], "--mismatches" ) ) mismatches = true;
        else if ( !strcmp( argv[i], "--collapse" ) ) collapse = true;
    }
    
    ir_ = new IndexReader( fns );
    qb_ = new QueryBinaries( fns );
    
    if ( ofn.empty() ) ofn = "./match_result.fa";
//    errors = 10;
    if ( !ifn.empty() )
    {
//        ifstream tfs( "/home/glen/Genomes/sediment/query_results_dump.fa" );
//        string line[2];
//        vector<string> tests;
//        tests.push_back("TACCTGGTTGATCCTGCCAGTAGTCATATGCTTGTCTCAAAGATTAAGCCATGCATGTCTCAGTATAAGCTTTTGATTTGTGAAACTGCAAATGGCTCATTAAAACAGTTATAATGCATTTGACAATTGCTCCTACATGGATATCTGTGGAAATTCTATAGTTAATACATGCACTGAAAACCTATCTTGGGGGAAAGGTTGTGGTCGTTGGTGCAGAACCAATTCAAGCATTGCTTGTATATTTGATTGAATCGCAATGACAAAATGAATTACATGGCAACAGCTGGTAATAATTCATGCCAGTTTCTGACCTATCAGCTTTCGACGGTAAGGTATTGGCTTACCGTGGCAATGACAGGTGACGGAGGATTAGGGTTTGATTCCGGAGAGGGAGCTTGAGAAATGGCTACCACATCTAAGGAAGGCAGCAGGCGCGCAAATTACCCAATCCTGGCACAGGGAGGTAGTGACAAGAAATAACAATACAAGGCATCCATGTCTTGTACTTGGAATGATTGGACATCAAACCTTTCTATAAGTATCAATTGGAGGGCAAGTCTGGTGCCAGCAGCCGCGGTAATTCCAGCTCCAATAGCGTATATTAAAGTTGTTGCGGTTAAAAAGCTCGTAGTTGGATTTCTGCTGAGGATGGCTGGTCCGCCCTCTGGGTGAGTATTTGGCACAGCCTGAGCATTTATCTTGAAAGCACAACTGCACTTGACTGTGTGGTGTGTTATTGAGAACATTTACTTTGAGGAAATCAGAGTGTTTCAAGCAGGTGTTTGGCCTTGAATACATTAGCATGGAATAATAATTAAGATCATGATTCTTTTTGTTAGTTTCTAGAATTGAGGTAATGATTAATAGGGATAGTTGGGGGCATTCGTATTTGATTGTCAGAGGTGAAATTCTTGGATTTGTTAAAGACGGACTACTGCGAAAGCATTTGCCAAGGATGTTTTCATTGATCAAGAACGAAAGTTAAGGGATCGAAGACGATCAGATACCGTCCTAGTCTTAACCATAAACCATGCCAACTAGAGATTGAAGGTTGTTACTCGTATGACTTTTTCAGCACCCTATGAGAAATCAAAGTGTTTGGGTTCCGGGGGGAGTATGGTCGCAAGGCTGAAACTTAAAGGAATTGACGGAAGGGCACCACCAGGAGTGGAGCCTGCGGCTTAATTTGACTCAACACGGGGAAACTTACCAGGTCCAGACATAATGAGGATTGACAGATTGATAGCTTTTTCTTGATTCTATGGGTGGTGGTGCATGGCCGTTCTTAGTTGGTGGAGTGATTTGTCTGGTTAATTCCGTTAACGAACGAGACCTTAACCTGCTAAATAGTTACATGTAATTTCGATTATGTGGGCAACTTCTTAGAGGGACTTTGTGTGCATAATGCAAGGAAGTTTGAGGCAATAACAGGTCTGTGATGCCCTTAGATGTTCTGGGCTGCACGCGCGCTACACTGATGTGTTCAACGAGTTTTTAACCTTGCCTGGAAAGGTTTGGTAATCTTGAACAGGCATCGTGATGGGGATTGTTTATTGTAATTATTAATCTTAAACGAGGAATTCCTAGTAAGCTTGAGTCATCAGCTTGTGCTGATTATGTCCCTGCCCTTTGTACACACCGCCCGTCGCTCCTACCGATTGAGTGATCCGGTGAATAATTTGGACTGTAGCAATGTTCAGTTCTTGAACAACGCAATGGAAAGTTTAATGAACCTTATCACTTAGAGGAAGGAGAAGTCGTAACAAGGTTTCTGTAGGTGAACCTGCAGAAGGATCATTT");
//        tests.push_back("CACCTGGTTGATCCTGCCAGTAGTCATATGCTTGTCTCAAAGATTAAGCCATGCATGTCTCAGTATAAGCTTCTACACGGCGAAACTGCGAATGGCTCATTAAAACAGTTATCGTTTATTTGGTGGTCATTCAGTACATGGATAACCATGGTAATTCTAGAGCTAATACATGCGCCCAGACCCAACTTTGTGGAAGGGTTGTGTTTATTAGGTACAGAACCAACCCAGGCTCTGCCTGGTCTTTTGGTGATTCATGATAACCGAACGAATCGCATGGCATCAGCTGGCGATGAATCATTCAAGTTTCTGACCTATCAGCTTCCGACGGTAGGGTATTGGCCTACCGTGGCAATGACGGGTAACGGAGAATTAGGGTTCGATTCCGGAGAGGGAGCCTGAGAAACGGCTACCACATCTAAGGAAGGCAGCAGGCGCGCAAATTACCCAATCCTGACACAGGGAGGTAGTGACAAGAAATAACAATACAGGGCATCCATGTCTTGTAATTGGAATGAATAGAATTTAAATCCCTTTATGAGTATCAATTGGAGGGCAAGTCTGGTGCCAGCAGCCGCGGTAATTCCAGCTCCAATAGCGTATATTAAAGTTGTTGCGGTTAAAAAGCTCGTAGTTGGATTTCTGTTGAGGACGATCGGTCCGCCCTCTGGGTGAGTATCTGGTTCGGCCTTGGCATCTTCTTGGAGAACGTATCTGCACTTGACTGTGTGGCGCGGTATCCAGGACTTTTACTTTGAGGAAATTAGAGTGTTTCAAGCAGGCGCACGCCTTGAATACATTAGCATGGAATAATAAGATAGGACCTTGGTTCTATTTTGTTGGTTTCTAGAGCTGAGGTAATGATTAATAGGGATAGTTGGGGGCATTCGTATTTAACTGTCAGAGGTGAAATTCTTGGATTTGTTAAAGACGGACTACTGCGAAAGCATTTGCCAAGGATGTTTTCATTGATCAAGAACGAAAGTTAGGGGATCGAAGACGATCAGATACCGTCCTAGTCTTAACCATAAACCATGCCGACTAGAGATTGGAGGCCGTTAACTATACGACTCCTTCAGCACCTTATGAGAAATCACAAGTTTCTGGGTTCCGGGGGGAGTATGGTCGCAAGGCTGAAACTTAAAGGAATTGACGGAAGGGCACCACCAGGAGTGGAGCCTGCGGCTTAATTTGACTCAACACGGGGAAACTTACCAGGTCCAGACATAGTAAGGATTGACAGATTGATAGCTCTTTCTTGATTCTATGGGTGGTGGTGCATGGCCGTTCTTAGTTGGTGGAGTGATTTGTCTGGTTAATTCCGTTAACGAACGAGACCTTAGCCTGCTAAATAGTTACGCATAACATCTGTTATGCGGGCAACTTCTTAGAGGGACTTTGTGTGTCTAACGCAAGGAAGTTTGAGGCAATAACAGGTCTGTGATGCCCTTAGATGTTCTGGGCTGCACGCGCGCTACACTGATGCACTCAACGAGTTTACGACCTTGTCCGACAGGATTGGGTAATCTTTTAAAAATGCATCGTGATGGGGATAGATTATTGCAATTATTAATCTTCAACGAGGAATTCCTAGTAAGCGCGAGTCATCAGCTCGTGCTGATTACGTCCCTGCCCTTTGTACACACCGCCCGTCGCTCCTACCGATTGAGTGATTCGGTGAATAATTCGGACTGCAGCAATGTTTGGATCCCGAACGTTGCAGCGGAAAGTTTAGTGAACCTCATCACTTAGAGGAAGGAGAAGTCGTAACAAGGTTTCCGTAGGTGAACCTGCAGAAGGATCAAGC");
//        tests.push_back("CACCTGGTTGATCCTGCCAGTAGTCATATGCTTGTCTCAAAGATTAAGCCATGCATGTCTCAGTATCAACTTTTACACGGTAACACTGCGAATGGCTCATTAAAACAGTTATAGTTTATTTGATGGCATTCATTATATGGATACCTTTGGTAATTCTAGAGCTAATACATGCATACAGACCCGACTCATGAAGGGTTGTGTTTATTAGATTCCAAACCTTCCAACTTCGGTTGTGTTCTTGGTGATTCATGATAAAATACTGATCGCCTTTGGCGACGAGTCTTTCAAGTTTCTGACCTATCAGCTTTCGATGGTAGGGTATTGGCCTACCATGGCAGTGACGGGTAACAGAGAATTAGGGTTCGATTCTGGAGAGGGAGCCTGAGAAACGGCTACCACTTCTAAGGAAGGCAGCAGGCGCGTAAATTACCCAATCCTGACACAGGGAGGTAGTGACAAGAAATAACAATACAGGGCATCCATGTCTTGTAATTGGAATGAGTAAAATTTAAATCTCTCTACGAGTATCAATTGGAGGGCAAGTCTGGTGCCAGCAGCCGCGGTAATTCCAGCTCCAATAGCGTATATTAAAGTTGTTGCGGTTAAAAAGCTCGTAGTTGGATTTCTGCCGAGGACAGTCGGTCCGCCCTTGCGGGTGAGCATCTGAATGGCCCTGGGCATCTTTTCCGAAGACTGTATCTGCACTTCATTGTGTGGAGCGGTATTCGGAACCTTTACTTTGAGGAAATTAGAGTGTTTCAGGCATGCGCTTGTGGTGAATACATTAGCATGGAATAATAAGATAGGACTCGTGGCCTACTTTGTTGGTTTGTGTGCCATAGGTAATGATTAATAGGGATAGTTGGGGGCATTCGTATTTAATTGTCAGAGGTGAAATTCTTGGATTTGTTAAAGACGAACTACTGCGAAAGCATTTGCCAAGGATGTTTTCATTGATCAAGAACGAAAGTTAGGGGATCGAAGACGATCAGATACCGTCGTAGTCTTAACCATAAACCATGCCGACTAGAGATTGGGGGTCGTTAGCTTTATGACTCCTTCAGCACCTTATGAGAAATCAAAGTTTTTGGGTTCTGGGGGGAGTATGGTCGCAAGGCTGAAACTTAAAGGAATTGACGGAAGGGCACCACCAGGAGTGGAGCCTGCGGCTTAATTTGACTCAACACGGGAAAACTCACCAGGTCCAGACATAGTAAGGATTGACAGATTGATAGCTCTTTCTTGATTCTATGGGTGGTGGTGCATGGCCGTTCTTAGTTGGTGGAGTGATTTGTCTGGTTAATTCCGTTAACGAACGAGACCTTAGCCTGCTAAATAGTGTCAAGTATGATTCTTATTTGATTACTTCTTAGAGGGACTTTGTGTGTCTAACGCAAGGAAGTTTGAGGCAATAACAGGTCTGTGATGCCCTTAGATGTTCTGGGCTGCACGCGCGCTACACTGATGCATTCAACGAGTTTATAACCTTGCCTGAAAAGGTTGGGTAATCTGCAATGTGCATCGTGATGGGGATAGATTATTGCAATTATTAATCTTCAACGAGGAATTCCTAGTAAGCGCGAGTCATCAGCTCGTGCTGATTACGTCCCTGCCCTTTGTACACACCGCCCGTCGCTCCTACCGATTGAGTGATTCGGTGAATAATTCGGAGATTGATTTGTTCAGCTTTGCAGAACATGTCTTTCGAAGTCTAGTGAACCTCATCACTTAGAGGAAGGAGAAGTCGTAACAAGGTTTCCGTAGGTGAACCTGCAGAAGGATCAAGC");
//        tests.push_back("CTACCTGGTTGATCCTGCCAGTAGTCATATGCTTGTCTCAAAGATTAAGCCATGCATGTCTAAGTATAAGCGACTATACTGTGAAACTGCGAATGGCTCATTAAATCAGTTATGGTTTATTTGATGGTACCTTGCTACTTGGATAACCGTAGTAATTCTAGAGCTAATACATGCAGGAGTTCCCGACTCACGGAGGGATGTATTTATTAGATAAGAAACCAAACCGGTCTCCGGTTGCGTGCTGAGTCATAATAACTGCTCGAATCGCACGGCTCTACGCCGGCGATGGTTCATTCAAATTTCTGCCCTATCAGCTTTCGATGGTAGGATAGAGGCCTACCATGGCGTTAACGGGTAACGGAGAATTAGGGTTCGATTCCGGAGAGGGAGCCTGAGAAATGGCTACCACATCCAAGGAAGGCAGCAGGCGCGTAAATTGCCCGAATCCTGACACAGGGAGGTAGTGACAAGAAATAACAATACAGGGCTATTTTAGTCTTGTAATTGGAATGAGTACAATTTACATCTCTTCACGAGGATCAATTGGAGGGCAAGTCTGGTGCCAGCAGCCGCGGTAATTCCAGCTCCAATAGCGTATATTAAAGTTGTTGCAGTTAAAACGCTCGTAGTCGGATTTCGGGGCGGGCCGACCGGTCTGCCGATGGGTATGCACTGGCCGGCGCGTCCTTCCACCCGGAGACCGCGCCTACTCTTAACTGAGCGGGCGCGGGAGACGGGTCTTTTACTTTGAAAAAATCAGAGTGTTTCAAGCAGGCAGTCGCTCTTGCATGGATTAGCATGGGATAATGAAATAGGACTCTGGTGCTATTTTGTTGGTTTCGAACACCGGAGTAATGATTAACAGGGACAGTCAGGGGCACTCGTATTCCGCCGAGAGAGGTGAAATTCTCAGACCAGCGGAAGACGAACCACTGCGAAAGCATTTGCCAGGGATGTTTTCACTGATCAAGAACGAAAGTTAGGGGATCGAAGACGATCAGATACCGTCGTAGTCTTAACCATAAACCATGCCGACTAGGGATTGGAGGATGTTCCATTTGTGACTCCTTCAGCACCTTTCGGGAAACTAAAGTCTTTGGGTTCCGGGGGGAGTATGGTCGCAAGGCTGAAACTTAAAGGAATTGACGGAAGGGCACCACCAGGAGTGGAGCCTGCGGCTTAATTTGACTCAACACGGGGAAACTTACCAGGTCCAGACATTGTGAGGATTGACAGATTGAGAGCTCTTTCTTGATTCGATGGGTGGTGGTGCATGGCCGTTCTTAGTTGGTGGAGTGATTTGTCTGGTTAATTCCGTTAACGAACGAGACCGCAGCCTGCTAAATAGCGACGCGAACCCTCCGTTCGCTGGAGCTTCTTAGAGGGACAACTTGTCTTCAACAAGTGGAAGTTCGCGGCAATAACAGGTCTGTGATGCCCTTAGATGTTCTGGGCCGCACGCGCGCTACACTGATGCACTCAACGAGTCTATCACCTTGACCGAGAGGTCCGGGTAATCTTTTGAAATTGCATCGTGATGGGGATAGATTATTGCAACTATTAATCTTCAACGAGGAATTCCTAGTAAGCGTGTGTCATCAGCGCACGTTGATTACGTCCCTGCCCTTTGTACACACCGCCCGTCGCTCCTACCGATTGAATGATCCGGTGAGGCCCCCGGACTGCGGCGCCGCAGCTGGTTCTCCAGCCGCGACGCCGCGGGAAGCTGTCCGAACCTTATCATTTAGAGGAAGGAGAAGTCGTAACAAGGTTTCCGTAGGTGAACCTGCGG");
//        tests.push_back("CTGCCAGTAGTCATATGCTTGTCTCAAAGATTAAGCCATGCATGTCTAAGTACAGGCCTAAGCAAGGCGAAACCGCGAATGGCTCATTAAATCAGTTATGGTTCATTGGATCTTACCCTTTCTTACTTGGATAACTGTGGTAATTCTAGAGCTAATACATGCCAACCAGTAGCTCCGACATACTCCCGTATGGAAGAGCGCTTTTATTAGATCAAAACCAGTCCTCTTTCGGGAGGTCAAACCCCTTTGTGGTGACTCTGAATAACTTTTAGCTGAGCGCACGGTCTTCGTGCCGGCGCCGCATCTTTCAAGTGTCTGCCTTATCAGCTTTCGATTGTAGGTTATGCGCCTACAATGGCTATAACGGGTAACGGGGAATCAGGGTTCGATTCCGGAGAGGGAGCCTGAGAAACGGCTACCACATCTAAGGAAGGCAGCAGGCACGCAAATTACCCACTCCCGGCTCGGGGAGGTAGTGACGAAAAATAACGATGCGGGACTCATCCGAGGCCTCGCAATCGGAATGAGTACACTTTAAATCCTTTAACGAGGATCTATTGGAGGGCAAGTCTGGTGCCAGCAGCCGCGGTAATTCCAGCTCCAATAGCGTATATTAAAGTTGTTGCGGTTAAAAAGCTCGTAGTTGGATCTCAGTTCCGGACCGACGGTACGCCGCCCGGCGTCTACTGTCGGCTTCCGAACGCTAAAGCCCGTCGGCTCGGTCGGGGTGCTCTTCATCGAGTGTCCCGGGCGGCCGGCATGTTTACTTTGAAAAAATTAGGGTGCTCAGAGCAGGCCAAGGCGCCTGAATACTTGTGCATGGAATAATGGAATAGGACCTCGGTTCTATTTTGTTGGTTTTAAGAGCCCGAGGTAATGATTAATAGGAACAGACGGGGGCATTCGTATTGCGACGTTAGAGGTGAAATTCTTGGATCGTCGCAAGACGAACTACTGCGAAAGCATTTGCCAAGAATGTTTTCATTAATCAAGAACGAAAGTTAGAGGTTCGAAGGCGATCAGATACCGCCCTAGTTCTAACCATAAACGATGCTGACTAGCGATCCGCCGGCGTTAATCCCATGACTCGGCGGGCAGCTTCCCGGGAAACCAAAGTCTCTGAGTTCCGGGGGAAGTATGGTTGCAAAGCTGAAACTTAAAGGAATTGACGGAAGGGCACCACCAGGAGTGGAGCCTGCGGCTTAATTTGACTCAACACGGGGAACCTCACCAGGCCCGGACACCGGAAGGATTGACAGATTGAGAGCTCTTTCTCGATTCGGTGGGTGGTGGTGCATGGCCGTTCTTAGTTGGTGGAGCGATTTGTCTGGTTAATTCCGATAACGAACGAGACTCTGGCCTACTAACTAGTCGTCGGGTTGTCTTGTACCGGCTTCGGCCGGTCGCCCGTCGCAACTTCTTCTTAGAGGGATCCGCGGCATCTAGCCGCAAGAGACTGAGCAATAACAGGTCTGTGATGCCCTTAGATGTTCTGGGCCGCACGCGCGCTACACTGAAGGGATCAGCGTGTCCTTCCCTCTCCGAGAGGAGCGGGTAACCCGTTGAAACCCCTTCGTGATAGGGATTGGGGCTTGAAATTTTTCCCCATCAACGAGGAATTCCCAGTAAGCGCGAGTCATCAGCTCGCGTTGATTACGTCCCTGCCCTTTGTACACACCGCCCGTCGCTACTACCGATTGAATGATTTAGTGAGGCCTCGGGACTGGCGTTCCTGGCGGCCTTGGTCGCCTCGAGCTGCCGGAAACATGTCCAAACTTGATCATTTAGAGGAAGTAAAAGTCGTA");
//        tests.push_back("CTGGTTGATCCTGCCAGTAGTCATACGCTCGTCTCAAAGATTAAGCCATGCATGTCTAAGTATAAATATTTTACTTTGAAACTGCGAACGGCTCATTATATCAGTTATAGTTTATTTGATAGTCCCTTACTACTTGGATACCCGTAGTAATTCTAGAGCTAATACATGCGTCAATACCCTTCTGGGGTAGTATTTATTAGATTGAAACCAACCCCTTCGGGGTGATGTGGTGATTCATAATAAGCTTGCGGATCGCATGCCTTTGGCGGCGATGGATCATTCAAGTTTCTGCCCTATCAGCTTTGGATGGTAGGGTATTGGCCTACCATGGCTTTAACGGGTAACGGGAAATTAGGGTTTGATTCCGGAGAGGGAGCCTGAGAGACGGCTACCACATCCAAGGAAGGCAGCAGGCGCGTAAATTACCCAATCCTGACACAGGGAGGTAGTGACAATAAATAACAATGCCGGGCCTTCTTAGGTCTGGCAATTGGAATGAGAACAATTTAAACCCCTTATCGAGTATCAATTGGAGGGCAAGTCTGGTGCCAGCAGCCGCGGTAATTCCAGCTCCAATAGCGTATATTAAAGTTGTTGCAGTTAAAAAGCTCGTAGTTGGATTTGTGGTGTGTCCAGTTGGCCTTTGCTCTTTGAGTGATAGTGCTTTACTGGTCTGCCATGTTTGGGTGGAATCTGTGTGGCATTAAGTTGTCGTGCAGGGGATGCCCATCGTTTACTGTGAAAAAATTAGAGTGTTCAAAGCAGGCTTATGCCTCTGAATATATTAGCATGGAATAATGATATAGGACCTTGGTACTATTTTGTTGGTTTGCGCACTGAGGTAATGATTAAGAGGGACAGTTGGGGGTATTTGTATTCCATTGTCAGAGGTGAAATTCTTGGATTTTTGGAAGACAAACTACTGCGAAAGCATTTACCAAGGATGTTTTCATTAATCAAGAACGAAAGTTAGGGGATCGAAGATGATTAGATACCATCGTAGTCTTAACCATAAACTATGCCGACAAGGGATTGGTGGAGTTTCGTTTCGTCTCCATCAGCACCTTGTGAGAAATCATAAGTCTTTGGGTTCCGGGGGGAGTATGGTCGCAAGGCTGAAACTTAAAGAAATTGACGGAAGGGCACCACCAGGAGTGGAGCCTGCGGCTTAATTTGACTCAACACGGGAAAACTTACCAGGTCCAGACATAGTGAGGATTGACAGATTGAGAGCTCTTTCTTGATTCTATGGGTGGTGGTGCATGGCCGTTCTTAGTTGGTGGAGTGATTTGTCTGGTTAATTCCGTTAACGAACGAGACCCCTGCCTGCTAAATAGCACGTGTAGTGTTTATCACTGCATTTGTGCTTCTTAGAGGGACGTGCGTTCTATTAGACGCAGGAAGATAGGGGCAATAACAGGTCTGTGATGCCCTTAGATGTTCTGGGCCGCACGCGCGCTACACTGATGCATTCAACGAGTTCTTACCTTGGCCGAGAGGCCTGGGCAATCTTTTGAACTTGCATCGTGATAGGGATAGATTATTGCAATTATTAATCTTGAACGAGGAATTCCTAGTAAACGCAGATCATCAATCTGCATTGATTACGTCCCTGCCCTTTGTACACACCGCCCGTCGCACCTACCGATTGAATGGTCCGGTGAGGCCTCGGGATTGTGGTTAGTTTCCTTTATTGGAAGTTAGTCGCGAGAACTTGTCCAAACCTTATCATTTAGAGGAAGGTGAAGTCGTAACAAGGTTTCCGTGGTGAACCTGCAG");
//        KmerScaffold ks( 12 );
//        for ( int i = 0; i < tests.size(); i++ ) ks.add( i, tests[i] );
//        while ( getline( tfs, line[0] ) && getline( tfs, line[1] ) )
//        {
//            ks.addMapped( stoll( line[0].substr( 1 ) ), line[1] );
//        }
//        ks.test();
//        assert( false );
        ifstream ifs( ifn );
//        vector<MatchedQuery> queries;
        int i = 0;
        Result result;
//        ofstream ofs( ofn );
        while ( getSeq( ifs, header, seq ) ){
            Target* tar = result.addTarget( header, seq );
            for ( Read r : MatchQuery( seq, ir_, errors ).yield( qb_ ) ) result.addMatch( tar, r.id_, r.seq_, r.coords_[0] );
            tar->assemble();
//            tar->print( ofs );
            assert( false );
//            queries.push_back( MatchedQuery( header, seq, ir_, qb_, errors ) );
        }
        assert( false );
//        MatchedQuery::compete( queries );
        
//        {
//            unordered_set<ReadId> blocked, printed;
//            AmpliconReduction ar;
//            for ( MatchedQuery& mq : queries ) for ( Read r : mq.exact_ ) ar.add( r.id_, r.seq_ );
//            for ( MatchedQuery& mq : queries ) for ( MatchRead mr : mq.inexact_ ) ar.add( mr.id_, mr.seq_ );
//            blocked = ar.getBlocked();
//            for ( MatchedQuery& mq : queries ) for ( int i = 0; i < mq.exact_.size(); i++ ) if ( blocked.find( mq.exact_[i].id_ ) != blocked.end() ) mq.exact_.erase( mq.exact_.begin() + i-- );
//            for ( MatchedQuery& mq : queries ) for ( int i = 0; i < mq.inexact_.size(); i++ ) if ( blocked.find( mq.inexact_[i].id_ ) != blocked.end() ) mq.inexact_.erase( mq.inexact_.begin() + i-- );
//            
//            
//        }
        
//        KmerGraph kg( 12 );
//        for ( Read r : queries[0].exact_ ) kg.add( r.id_, r.seq_ );
//        for ( MatchRead mr : queries[0].inexact_ ) kg.add( mr.id_, mr.seq_ );
//        kg.assemble();
//        output( ofn, queries, true, true );
        assert( false );
    }
    else if ( !seq.empty() )
    {
    
        ofstream ofs( ofn );
        header = "query";
        match( seq, header, ofs.good() ? &ofs : NULL, min( 15, errors ) );
        if ( ofs.good() ) ofs.close();
    }
    else
    {
        cerr << "Error: no query sequence(s) provided." << endl;
        exit( EXIT_FAILURE );
    }
}

void Assemble::match( string& q, string& header, ofstream* ofs, int errors )
{
    vector<Read> reads = MatchQuery( q, ir_, errors ).yield( qb_ );
    Read::sort( reads, true, 0 );
    int base = !reads.empty() ? max( -reads[0].coords_[0], 0 ) : 0;
    if ( ofs ) ( *ofs ) << ">" << header << "|matched:" << reads.size() << endl << string( base, '-' ) << "reads" << q << endl;
    else cout << ">" << header << "|matched:" << reads.size() << endl << string( base, '-' ) << "reads" << q << endl;
    Lib* lib;
    unordered_set<ReadId> used;
    bool addPairs = false;
    for ( Read& r : reads ) if ( used.find( r.id_ ) == used.end() )
    {
        string header = ">read_" + to_string( r.id_ );
        string seq = r.seq_;
        int d, coord = r.coords_[0]+base;
        if ( addPairs )
        {
            lib = params.getLib( r.id_ );
            if ( lib && lib->isPe )
            {
                ReadId id = r.id_;
                lib->getPair( id, d );
                string alt = qb_->getSequence( id );
                int ol = mapSeqOverlap( seq, alt, 15, d );
                if ( ol )
                {
                    if ( !d ) coord -= alt.size()-ol;
                    seq = (d?seq:alt) + (d?alt:seq).substr( ol );
                }
                else
                {
                    seq = (d?seq:alt) + string( 100, '-' ) + (d?alt:seq);
                    if ( !d ) coord -= 250;
                }
                used.insert( id );
            }
        }
        assert( coord >= 0 );
        seq = string( coord, '-' ) + seq;
        if ( ofs ) ( *ofs ) << header << endl << seq << endl;
        else cout << header << endl << seq << endl;
    }
}

//void Assemble::output( string ofn, vector<MatchedQuery>& queries, bool exact, bool inexact )
//{
//    ofstream ofs( ofn );
//    unordered_set<ReadId> used;
//    for ( MatchedQuery& mq : queries )
//    {
//        ofs << ">" + mq.header_ << "\n" << string( params.readLen, '-' ) << mq.seq_ << "\n";
//        if ( exact ) for ( Read& r : mq.exact_ ) if ( used.insert( r.id_ ).second ) ofs << ">Exact_match_" << r.id_ << "\n" << string( max( 0, r.coords_[0]+params.readLen ), '-' ) << r.seq_ << "\n";
//        if ( inexact ) for ( MatchRead& r : mq.inexact_ ) if ( used.insert( r.id_ ).second ) ofs << ">Inexact_match_" << r.id_ << "\n" << string( max( 0, r.query_[0]-r.read_[0]+params.readLen ), '-' ) << r.seq_ << "\n";
//        if ( inexact ) for ( Read& r : mq.unmatched_ ) if ( used.insert( r.id_ ).second ) ofs << ">Inexact_match_" << r.id_ << "\n" << string( max( 0, r.coords_[0]+params.readLen ), '-' ) << r.seq_ << "\n";
////        if ( exact ) for ( Read& r : mq.exact_ ) cout << "Exact_match_" << r.id_ << "\n" << string( max( 0, r.coords_[0]+params.readLen ), '-' ) << mq.seq_ << "\n";
////        if ( inexact ) for ( MatchRead& r : mq.inexact_ ) cout << "Inexact_match_" << r.id_ << "\n" << string( max( 0, r.query_[0]-r.read_[0]+params.readLen ), '-' ) << mq.seq_ << "\n";
////        if ( inexact ) for ( Read& r : mq.unmatched_ ) cout << "Inexact_match_" << r.id_ << "\n" << string( max( 0, r.coords_[0]+params.readLen ), '-' ) << mq.seq_ << "\n";
//    }
//    ofs.close();
//}

void Assemble::test( int tests, int errors )
{
    srand( time(NULL) );
    vector<MatchQuery> matches;
    vector< pair<ReadId, string> > queries;
    while ( queries.size() < tests )
    {
        ReadId id = ( ( rand() & 65535 ) << 16 | ( rand() & 65535 ) ) % params.seqCount;
        string seq = qb_->getSequence( id );
        if ( seq.size() == params.readLen ) queries.push_back( make_pair( id, seq ) );
    }
    
    double startTime = clock();
    for ( int i = 0; i < queries.size(); i++ ) matches.push_back( MatchQuery( queries[i].second, ir_, errors ) );
    cout << "Errors: " << errors << ", queryies: " << tests << ", time taken: " << getDuration( startTime ) << endl;
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