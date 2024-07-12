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

#ifndef LONG_QUERY_H
#define LONG_QUERY_H

#include "types.h"
#include "index_reader.h"
#include "query_binary.h"
#include "shared_structs.h"

struct ReadMatch
{
    struct Align
    {
        Align( size_t tar=0, size_t read=0, size_t len=0 ): tar_( tar ), read_( read ), len_( len ){};
        size_t gapLen( Align& rhs );
        size_t readEnd();
        size_t tarEnd();
        size_t totalLen( Align& rhs );
        size_t tar_, read_, len_;
    };
    ReadMatch( string seq, string& base, ReadId id, size_t coord, size_t len, bool drxn );
    bool assess( uint hitLen, uint errorRate );
    int score();
    char seq( int i );
    string& seq();
    size_t size();
    string seq_;
    ReadId id_;
    vector<Align> align_;
    size_t read_[2], tar_[2];
};

struct LongResult
{
    LongResult( vector<ReadMatch*>& reads ): reads_( reads ){};
    ~LongResult(){ for ( ReadMatch* rm : reads_ ) delete rm; };
    vector<ReadMatch*> reads_;
};

class LongQuery
{
    struct LQ
    {
        LQ( bool drxn ): drxn_( drxn ){};
//        LQ( size_t blockEnd, bool noRev, bool drxn ): blockEnd_( blockEnd ), noRev_( noRev ), drxn_( drxn ){};
//        LQ( string& seq, size_t p, size_t blockEnd, size_t maxLen, bool noRev, bool drxn );
        vector<uint8_t> q_, t_;
        vector<uint> blocks_, errBlocks_;
        size_t coords_[2], blockEnd_;
        bool drxn_;
    };
    struct LQState
    {
        LQState(): len_( 1 ), pBlock_( 0 ), pErr_( 0 ), errCount_( 0 ), errAllowed_( 0 ), retrieve_( false ){};
        size_t len_;
        uint pBlock_, pErr_, errCount_, errAllowed_;
        bool retrieve_;
    };
    struct EndRank
    {
        int score(){ return coords_[1] - coords_[0] - errors_; };
        size_t coords_[2], errors_;
    };
    void add( LQ& lq, LQState& lqe, CharId rank, CharId count );
    bool advance( LQ& lq, LQState& lqe, uint8_t c );
    ReadMatch* create( ReadId id, string seq, EndRank& rank, bool drxn );
    bool flush( QueryBinaries* qb, bool force=false );
    LQ initiate( CharId& rank, CharId& count, bool drxn );
    bool package( bool force );
    bool query();
    void query( LQ& lq, LQState lqe, CharId rank, CharId count, uint8_t c );
//    void query( LQ& lq, CharId rank, CharId count, uint8_t c, int len, int errCount, int errAllowed );
    void retrieve( CharId rank, CharId count );
    bool reverse( LQ& lq, int errLeft );
    bool setChar( LQ& lq, LQState& lqe, uint8_t c );
    bool setChar( LQ& lq, uint8_t c, int len, int& errCount, int& errAllowed );
    
    string seq_;
    IndexReader* ir_;
    vector<ReadMatch*> matches_, result_;
    unordered_map<CharId, EndRank> endRanks_[2];
    vector<uint> blocks_;
    size_t p_, coords_[2], flushLen_;
    uint hitLen_, maxLen_, errorRate_, tarErrors_, tarLen_, pBlock_;
public:
    LongQuery( string seq, IndexReader* ir, uint hitLen, uint errorRate );
    LongResult* yield( QueryBinaries* qb );
};

#endif /* LONG_QUERY_H */

