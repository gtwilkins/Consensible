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

#include "bubble_resolution.h"
#include "alignment.h"
#include <cassert>

BubbleResolution::BubbleResolution( Bubble* b, string& base )
{
    vector<SNPs*> snps;
    SnpAlignment align( b->template_, base, 0, b->template_.size(), b->start_, b->len_, snps, b->bubs_ );
    SnpAlignResult result = align.align( true, true );
    assert( false );
}