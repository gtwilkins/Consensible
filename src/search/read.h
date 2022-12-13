/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   glin_read.h
 * Author: glen
 *
 * Created on 18 May 2022, 11:27 AM
 */

#ifndef GLIN_READ_H
#define GLIN_READ_H

#include "types.h"

struct MappedRead
{
    MappedRead( ReadId id, string seq ):id_( id ), seq_( seq ){};
    string seq_;
    vector<Match*> matches_;
    ReadId id_;
};


#endif /* GLIN_READ_H */

