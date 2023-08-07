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

#ifndef SEQUENCE_FILE_H
#define SEQUENCE_FILE_H

#include "types.h"
#include <fstream>

struct InputSequence
{
    string header_, seq_;
};

struct SequenceFile
{
    SequenceFile( string ifn );
    bool getSeq( InputSequence& seq );
    vector<InputSequence> getSeqs();
    bool isSequence( string &s );
    
    ifstream ifs_;
    string ifn_, line_;
    bool ended_;
};

#endif /* SEQUENCE_FILE_H */

