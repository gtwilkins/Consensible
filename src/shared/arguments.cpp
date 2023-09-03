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
#include "filenames.h"
#include <string.h>
#include <iostream>
#include <cassert>
#include <dirent.h>
#include <unistd.h>

using namespace std;

Arguments::Arguments( int argc, char** argv )
: workDir_( "./consensible_index" ), outFolder_( "./" ), outPrefix_( "./consensible-out" ), pInput_( 0 ), reindex_( false ), cleanup_( false ), help_( false ), finished_( false )
{
    for ( int i ( 1 ); i < argc; i++ )
    {
        if ( !strcmp( argv[i], "-i" ) )
        {
            if ( i+1 >= argc || argv[i+1][0] == '-') error( "No inputs given with -i flag" );
            while ( i+1 < argc && argv[i+1][0] != '-' ) for ( string fn : getFilenames( argv[++i] ) ) addInput( fn );
        }
        else if ( !strcmp( argv[i], "-w" ) )
        {
            if ( i+1 >= argc || argv[i+1][0] == '-') error( "No inputs given with -w flag" );
            workDir_ = argv[++i] ;
        }
        else if ( !strcmp( argv[i], "-p" ) )
        {
            if ( i+1 >= argc || argv[i+1][0] == '-') error( "No inputs given with -p flag" );
            bwtPrefix_ = argv[++i] ;
        }
        else if ( !strcmp( argv[i], "-o" ) )
        {
            if ( i+1 >= argc || argv[i+1][0] == '-') error( "No inputs given with -o flag" );
            setOutFolder( argv[++i] );
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
    if ( help_ ) return;
    checkWorkingDir();
    setOutputs();
    finished_ = inputs_.empty();
}

void Arguments::addInput( string fn )
{
    for ( string ifn : inputs_ ) if ( fn == ifn ) return;
    inputs_.push_back( fn );
}

void Arguments::checkWorkingDir()
{
    if ( workDir_.empty() && bwtPrefix_.empty() ) error( "Either -w or -p arguments must be provided." );
    if ( !workDir_.empty() && !bwtPrefix_.empty() ) error( "-p and -w are mutually exclusive arguments." );
    if ( !workDir_.empty() )
    {
        if ( workDir_[0] != '/' )
        {
            string curr = getcwd( NULL, 0 );
            if ( !workDir_.empty() && workDir_[0] == '.' ) workDir_ = workDir_.substr( 1 );
            if ( !workDir_.empty() && workDir_[0] != '/' ) workDir_ = "/" + workDir_;
            if ( workDir_.empty() || curr.empty() )  error( "Invalid working directory. Please try using the absolute path." );
            workDir_ = curr + workDir_;
        }
        if ( workDir_.back() == '/' ) workDir_.pop_back();
        size_t it = workDir_.rfind( "/" );
        string folder = workDir_.substr( it < workDir_.size() ? it+1 : workDir_.size() );
        if ( folder != "consensible_index" ) workDir_ += "/consensible_index";
        Filenames::makeFolder( workDir_ );
        string line;
        ifstream dataIndex( workDir_ + "/dataset_index.txt" );
        if ( dataIndex ) while ( getline( dataIndex, line ) ) if ( !line.empty() && line[0] != '#' )
        {
            it = line.find( '=' );
            if ( it == string::npos || !it ||  it+1 >= line.size() ) continue;
            string indexed = line.substr( 0, it ), dataset = line.substr( it+1 );
            fileIndex_.push_back( make_pair( indexed, dataset ) );
        }
    }
}

string Arguments::getCode( int num )
{
    if ( num <= 0 ) error( "Unexpectedly tried to create a negatively-keyed index prefix?!" );
    
    string codestr;
    for ( int i = 0; i < 4 || num; i++ )
    {
        codestr += to_string( num % 10 );
        num /= 10;
    }
    return "cons" + string( codestr.rbegin(), codestr.rend() );
}

std::vector<std::string> Arguments::getFilenameParts( std::string filestr )
{
    size_t it = filestr.find( "/" );
    vector<std::string> parts{ filestr.substr( 0, it ) };
    for ( size_t pos = it; pos != string::npos && pos < filestr.size(); )
    {
        it = filestr.find( "/", pos+1 );
        parts.push_back( filestr.substr( pos+1, it-pos-1 ) );
        pos = it;
    }
    return parts;
}

vector<string> Arguments::getFilenames( string filestr )
{
    vector<std::string> filenames, parts = getFilenameParts( filestr );
    getFilenames( "", parts, filenames, 0 );
    return filenames;
}

void Arguments::getFilenames( string filestr, vector<string>& parts, vector<string>& filenames, int i )
{
    if ( i == parts.size() )
    {
        if ( !filestr.empty() ) filenames.push_back( filestr );
        return;
    }
    if ( i ) filestr += "/";
    int wild[2]{ !parts[i].empty() && parts[i][0] == '*'?1:0, !parts[i].size()>1 && parts[i].back() == '*'?1:0 };
    if ( wild[0] || wild[1] )
    {
        struct dirent *ent;
        DIR *dir = opendir( filestr.c_str() );
        if ( dir == NULL ) error( "Failed to open directory: " + filestr );
        string needle = parts[i].substr( wild[0], parts[i].size()-wild[0]-wild[1] );
        while ( ( ent = readdir(dir) ) != NULL )
        {
            string fn( ent->d_name );
            size_t pos = wild[0] ? fn.rfind( needle ) : fn.find( needle );
            if ( pos != string::npos && pos == ( wild[0] ? fn.size()-needle.size() : 0 ) ) getFilenames( filestr+fn, parts, filenames, i+1 );
        }
    }
    else
    {
        getFilenames( filestr+parts[i], parts, filenames, i+1 );
    }
}

bool Arguments::setBwtPrefix()
{
    if ( finished_ ) return false;
    
    input_.clear();
    if ( inputs_.empty() ) finished_ = true;
    
    if ( pInput_ < inputs_.size() )
    {
        outPrefix_ = outputs_[pInput_];
        int highest = 0;
        bool exists = false;
        for ( pair<string, string>& indexed : fileIndex_ )
        {
            if ( indexed.second == inputs_[pInput_] )
            {
                bwtPrefix_ = workDir_ + "/" + indexed.first + "/" + indexed.first;
                exists = true;
                break;
            }
            if ( indexed.first.find( "cons" ) == 0 )
            {
                string num = indexed.first.substr( 4 );
                if ( num.find_first_not_of( "0123456789" ) == string::npos ) highest = max( highest, stoi( num ) );
            }
        }
        if ( !exists )
        {
            string code = getCode( ++highest );
            while ( Filenames::isFolder( workDir_ + "/" + code ) ) code = getCode( ++highest );
            fileIndex_.push_back( make_pair( code, inputs_[pInput_] ) );
            bwtPrefix_ = workDir_ + "/" + code + "/" + code;
        }
        input_.push_back( inputs_[pInput_] );
        finished_ = ++pInput_ >= inputs_.size();
    }
    
    return true;
}

void Arguments::setOutFolder( std::string folder )
{
    size_t slash = folder.rfind( "/" ), dot = folder.rfind( "." );
    if ( dot != string::npos && ( slash == string::npos || slash < dot ) ) folder = folder.substr( 0, dot );
    while ( !folder.empty() && folder.back() == '/' ) folder.pop_back();
    if ( folder.empty() ) error( "Invalid output folder provided" );
    Filenames::makeFolder( folder );
    outFolder_ = folder + "/";
}

void Arguments::setOutputs()
{
    for ( string ifn : inputs_ )
    {
        size_t pos = ifn.rfind( "/" );
        string base = ifn.substr( pos == string::npos ? 0 : pos+1 );
        pos = base.rfind( "." );
        string ext = base.substr( pos );
        if ( ext == ".fa" || ext == ".fasta" || ext == ".fq" || ext == ".fastq" || ext == ".seq" ) base = base.substr( 0, pos );
        base = outFolder_ + base;
        string ofn = base;
        int dupe = 0;
        for ( int i = 0; i < outputs_.size(); i++ ) if ( outputs_[i] == ofn )
        {
            ofn = base + "-" + to_string( ++dupe );
            i = -1;
        }
        outputs_.push_back( ofn );
    }
}

void Arguments::updateFileIndex()
{
    if ( !pInput_ || fileIndex_.empty() ) return;
    ofstream ofs( workDir_ + "/dataset_index.txt" );
    ofs << "# This file lists datasets indexed by consensible\n";
    ofs << "# Retaining these index folders will make any subsequent queries of the same dataset much quicker\n";
    ofs << "# However, once you are done with a given dataset, you may safely delete its folder to free up disk space\n\n";
    for ( pair<string,string>& indexed : fileIndex_ ) ofs << indexed.first + "=" + indexed.second + "\n";
    ofs.close();
}

void Arguments::error( string msg )
{
    cerr << msg << endl;
    exit( EXIT_FAILURE );
}