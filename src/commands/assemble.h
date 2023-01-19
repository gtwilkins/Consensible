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

#ifndef ASSEMBLE_H
#define ASSEMBLE_H

#include "types.h"
#include "index_reader.h"
#include "query_binary.h"
#include "match_query.h"
#include "arguments.h"

class Assemble
{
public:
    Assemble( Arguments& args );
    ~Assemble();
    
private:
    void printUsage();
    IndexReader* ir_;
    QueryBinaries* qb_;
};


#endif /* MATCH_H */

