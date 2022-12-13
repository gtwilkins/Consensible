/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "consensus.h"
#include "target.h"
#include <cassert>
#include <algorithm>

Consensus::Consensus( vector<Match*>& matches, Target* tar )
: tar_( tar )
{
    coord_[0] = matches[0]->coords_[matches[0]->anchors_[0]];
    coord_[1] = matches[0]->coords_[matches[0]->anchors_[1]];
    for ( Match* gm : matches )
    {
        ConMap cm;
        cm.node_ = gm;
        cm.coord_[0] = gm->coords_[gm->anchors_[0]];
        cm.coord_[1] = gm->coords_[gm->anchors_[1]];
        cm.range_[0] = gm->anchors_[0];
        cm.range_[1] = gm->anchors_[1];
        coord_[0] = min( coord_[0], cm.coord_[0] );
        coord_[1] = max( coord_[1], cm.coord_[1] );
        maps_.push_back( cm );
    }
    sort( maps_.begin(), maps_.end(), []( ConMap& a, ConMap& b ){ return a.coord_[0] == b.coord_[0] ? a.coord_[1] > b.coord_[1] : a.coord_[0] < b.coord_[0]; } );
}
