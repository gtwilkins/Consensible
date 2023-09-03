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

#include "transform.h"
#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include <string.h>
#include <cassert>
#include <algorithm>
//#include <chrono>
//#include <iomanip>

void Transform::load( PreprocessFiles* fns, vector<string>& infilenames, bool revComp )
{
    assert( infilenames.size() == 1 );
    int minScore = 0;
    uint8_t readLen = 0, minLen = 25;
    vector<ReadFile*> infiles;
    for ( string& ifn : infilenames ) infiles.push_back( new ReadFile( ifn, 0, minScore ) );
    for ( ReadFile* rf : infiles ) readLen = max( readLen, rf->readLen );
    
    cout << "    Read length set to " << to_string( readLen ) << "." << endl;
    cout << "    Min length set to " << to_string( minLen ) << "." << endl;
    
    ReadId readCount = 0, discardCount = 0;
    double readStartTime = clock();
    
    BinaryWriter* binWrite = new BinaryWriter( fns, 0, readLen, revComp );
    
    for ( ReadFile* rf : infiles )
    {
        ReadId thisReadCount = 0, thisDiscardCount = 0;
        string read;
        while ( rf->getNext( read ) )
        {
            if ( read.length() >= minLen ) binWrite->write( read );
            else discardCount++;
            thisReadCount++;
        }
        delete rf;
        readCount += thisReadCount;
        readCount += thisDiscardCount;
        
        cout << "    Found " << to_string( thisReadCount-discardCount ) << " useable reads in file, discarded " << to_string( discardCount ) << " short reads." << endl;
    }
    
    binWrite->close();
    delete binWrite;
}

void Transform::run( PreprocessFiles* fns )
{
    cout << "Preprocessing step 2 of 3: transforming sequence data..." << endl << endl;
    
    BinaryReader* bin = new BinaryReader( fns );
    BwtCycler* cycler = new BwtCycler( fns );
    double totalStart = clock();
//    auto t_start = std::chrono::high_resolution_clock::now();
    
    while ( bin->cycle < bin->readLen )
    {
        double cycleStart = clock();
        cout << "    Cycle " << to_string( bin->cycle ) << " of " << to_string( bin->readLen ) << "... " << flush;
        
        bin->read();
        cycler->run( bin->chars, ( bin->anyEnds ? bin->ends : NULL ), bin->cycle );
        bin->update();
        
        cout << " completed in " << getDuration( cycleStart ) << endl;
    }
    
    double finalStart = clock();
    cout << "    Cycle " << to_string( bin->cycle ) << " of " << to_string( bin->readLen ) << "... " << flush;
    cycler->finish( bin->cycle + 1 );
    cout << " completed in " << getDuration( finalStart ) << endl;
    bin->finish();
    fns->clean();
    delete bin;
    delete cycler;
    
    cout << endl << "Transforming sequence data... completed!" << endl;
    cout << "Time taken: " << getDuration( totalStart );
//    cout << "   " << std::fixed << std::setprecision(2) << ( clock() - totalStart ) / CLOCKS_PER_SEC << " vs " << ( ( std::chrono::high_resolution_clock::now() - t_start ).count() / 1000.0 ) / CLOCKS_PER_SEC << endl << endl;
    cout << endl << endl;
}
