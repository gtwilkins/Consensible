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

#include "transform_files.h"
#include <iostream>

void TransformFile::checkFile( string mode )
{
    if ( fp_ == NULL )
    {
        if ( mode == "wb" ) cerr << "Error creating output file \"" << fn_ << "\" for writing." << endl;
        else if ( mode == "ab" ) cerr << "Error opening output file \"" << fn_ << "\" for appending." << endl;
        else if ( mode == "rb+" ) cerr << "Error opening output file \"" << fn_ << "\" for editing." << endl;
        else if ( mode == "rb" ) cerr << "Error opening input file \"" << fn_ << "\" for reading." << endl;
        else  cerr << "Error opening file \"" << fn_ << "\"." << endl;
        exit( EXIT_FAILURE );
    }
}

//void TransformFile::open()
//{
//    fp_ = fopen( fn_.c_str(), ( read_ ? "rb" : ( open_ ? "ab" : "wb" ) ) );
//    if ( !open_ )
//    {
//        if ( read_ ) fread( &totalSize_, 4, 1, fp_ );
//        else fwrite( &totalSize_, 4, 1, fp_ );
//        pos_ = 4;
//        open_ = true;
//    }
//    else if ( read_ ) fseek( fp_, pos_, SEEK_SET );
//    
//}
//
//void TransformFile::close()
//{
//    if ( !read_ )
//    {
//        
//    }
//    fclose( fp_ );
//}


TransFileLarge::~TransFileLarge()
{
    if ( buff_ ) delete buff_;
}

void TransFileLarge::flush()
{
    if ( p_ ) write();
    fp_ = fopen( fn_.c_str(), "rb+" );
    checkFile( "rb+" );
    fwrite( &totalSize_, 4, 1, fp_ );
    fclose( fp_ );
}

void TransFileLarge::set( std::string fn, bool reader )
{
    fn_ = fn;
    p_ = totalSize_ = curSize_ = 0;
    pos_ = 0;
    read_ = reader;
    open_ = false;
    if ( !buff_ ) buff_ = new uint32_t[buffSize_];
    
}

void TransFileLarge::read()
{
    fp_ = fopen( fn_.c_str(), "rb" );
    checkFile( "rb" );
    if ( !open_ && ( open_ = true ) ) fread( &totalSize_, 4, 1, fp_ );
    else fseek( fp_, pos_, SEEK_SET );
    curSize_ = min( totalSize_, buffSize_ );
    totalSize_ -= curSize_;
    fread( buff_, 4, curSize_, fp_ );
    pos_ = ftell( fp_ );
    fclose( fp_ );
    p_ = 0;
}

uint32_t TransFileLarge::readNext()
{
    if ( p_ == curSize_ ) read();
    return buff_[p_++];
}

void TransFileLarge::write()
{
    fp_ = fopen( fn_.c_str(), "rb+" );
    checkFile( "rb+" );
    if ( !open_ && ( open_ = true ) ) fwrite( &totalSize_, 4, 1, fp_ );
    else fseek( fp_, ( totalSize_*4 ) + 4, SEEK_SET );
    
    totalSize_ += p_;
    fwrite( buff_, 4, p_, fp_ );
    fclose( fp_ );
    p_ = 0;
}

void TransFileLarge::writeNext( uint32_t b )
{
    if ( p_ == buffSize_ ) write();
    buff_[p_++] = b;
}

