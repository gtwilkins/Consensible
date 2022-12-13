/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "result.h"

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

