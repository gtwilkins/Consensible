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

#include "result.h"
#include "fstream"
#include <iostream>

Result::~Result()
{
    for ( auto read : reads_ ) delete read.second;
}

MappedRead* Result::addRead( ReadId id, string seq )
{
    auto it = reads_.find( id );
    if ( it != reads_.end() ) return it->second;
    MappedRead* read = new MappedRead( id, seq );
    reads_.insert( make_pair( id, read ) );
    return read;
}

void Result::addMatch( Target* tar, ReadId id, string seq, int coord )
{
    tar->addMatch( addRead( id, seq ), coord );
}

Target* Result::addTarget( string header, string seq )
{
    Target* tar = new Target( header, seq );
    targets_.push_back( tar );
    return tar;
}

void Result::assemble( string& outPrefix )
{
    vector<Consensus*> consensus;
    for ( Target* tar : targets_ )
    {
        for ( Consensus* c : tar->assemble() ) consensus.push_back( c );
    }
    string ofn = outPrefix + "_consensus.fa";
    size_t it = outPrefix.find_last_of( '.' );
    if ( it != string::npos )
    {
        string excess = outPrefix.substr( it+1 );
        if ( excess == "fasta" || excess == "fa" ) ofn = outPrefix;
    }
    
    cout << "    " << to_string( reads_.size() ) << " reads were found to match the query sequence." << endl;
    cout << "    " << to_string( consensus.size() ) << " consensus sequences were assembled!" << endl;
    cout << endl << "Writing results to: " << ofn << endl;
    ofstream ofs( ofn );
    for ( int i = 0; i < consensus.size(); i++ )
    {
        string seq = consensus[i]->resolve();
        ofs << ">consensus_" + to_string( i+1 ) + "\n";
        ofs << seq + "\n";
    }
    ofs.close();
}

void Result::outputFullAlign( string& outPrefix )
{
    for ( int i = 0; i < targets_.size(); i++ )
    {
        string ofn = outPrefix + "-fullalign" + ( targets_.size() > 1 ? "-" + to_string( i+1 ) : "" ) + ".fa";
        ofstream ofs( ofn );
        targets_[i]->print( ofs );
        ofs.close();
    }
}

