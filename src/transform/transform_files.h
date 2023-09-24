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

#ifndef TRANSFORM_FILES_H
#define TRANSFORM_FILES_H

#include "types.h"
#include "transform_constants.h"
#include <fstream>

struct TransformFile
{
    TransformFile(): fp_( NULL ), pos_( 0 ), p_( 0 ), open_( false ){};
    void checkFile( string mode );
    void open();
    void close();
//    virtual void read();
//    virtual void write();
    std::string fn_;
    FILE* fp_;
    uint64_t pos_;
    uint32_t p_, buffSize_, totalSize_, curSize_;
    size_t bytes_;
    bool read_, open_;
};

//struct TransFileSmall : TransformFile
//{
//    void set();
//    void read();
//    void write();
//    uint8_t* buff_;
//};

struct TransFileLarge : TransformFile
{
    TransFileLarge():TransformFile(), buff_( NULL ){ buffSize_ = 64*1024; bytes_ = 4; };
    ~TransFileLarge();
    void flush();
    void set( string fn, bool reader );
    void read();
    void write();
    uint32_t readNext();
    void writeNext( uint32_t b );
    uint32_t* buff_;
};

#endif /* TRANSFORM_FILES_H */

