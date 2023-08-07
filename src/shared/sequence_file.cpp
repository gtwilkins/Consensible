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

#include "sequence_file.h"
#include <cassert>
#include <string.h>

SequenceFile::SequenceFile( string ifn )
: ifs_( ifn ), ifn_( ifn )
{
    ended_ = !getline( ifs_, line_ );
}

bool SequenceFile::getSeq( InputSequence& seq )
{
    if ( line_.empty() || ended_ ) return false;
    if ( line_[0] == '>' )
    {
        seq.header_ = line_.substr( 1 );
        seq.seq_ = "";
        assert( getline( ifs_, line_ ) && isSequence( line_ ) );
//        if ( !getline( ifs_, line_ ) || !isSequence( line_ ) ) return false;
        while ( isSequence( line_ ) )
        {
            seq.seq_ += line_;
            if ( ended_ = ( !getline( ifs_, line_ ) ) ) break;
        }
        
    }
    else if ( line_[0] == '@' )
    {
        seq.header_ = line_.substr( 1 );
        assert( getline( ifs_, line_ ) && isSequence( line_ ) );
//        if ( !getline( ifs_, line_ ) || !isSequence( line_ ) ) return false;
        seq.seq_ = line_;
        assert( getline( ifs_, line_ ) && line_[0] == '+' );
        assert( getline( ifs_, line_ ) && line_.size() == seq.seq_.size() );
        ended_ = !getline( ifs_, line_ );
    }
    else if ( !isSequence( line_ ) ) return false;
    else
    {
        seq.header_ = "";
        seq.seq_ = line_;
    }
    
//    if ( seq.header_.empty() && seq.seq_.empty() )
//    {
//        
//    }
    
    return true;
    
}

bool SequenceFile::isSequence( string &s )
{
    for ( char c : s )  if ( !strchr( "ACGTN", c ) ) return false;
    return true;
}
