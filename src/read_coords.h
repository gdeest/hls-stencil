#ifndef READ_COORDS_H
#define READ_COORDS_H

#include "hls_types.h"

#define COORDS_BSIZE 64

void read_coords(int ntiles,
                 int *coords,
                 hls::stream<tile_coord_t> &fi_ctrl2read,
                 hls::stream<tile_coord_t> &fi_ctrl2write);

#endif
