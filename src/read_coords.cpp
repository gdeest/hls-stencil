#include "read_coords.h"

void read_coords
( int ntiles,
  int *coords,
  hls::stream<tile_coord_t> &fi_ctrl2read,
  hls::stream<tile_coord_t> &fi_ctrl2write)
{
  int buff[3 * COORDS_BSIZE];

  int read      = 0;
  int remaining = ntiles;

  #ifdef DEBUG
  int n = 0;
  #endif

  // Continue while there are still coordinates to fetch.
  while (remaining > 0) {
    // Read min(remaining, COORDS_BSIZE) coordinates.
    int to_read = remaining > COORDS_BSIZE ? COORDS_BSIZE : remaining;

    // Perform burst. Start address is determined from the number of read coordinates.
    memcpy(buff, coords + 3*read, 3*to_read*sizeof(int));

    // Update counters.
    read      += to_read;
    remaining -= to_read;

    // Send coordinates to other actors.
    for (int i=0; i<to_read; i++) {
      #pragma HLS PIPELINE II=1
      tile_coord_t coord;
      coord.t  = buff[3*i];
      coord.x0 = buff[3*i + 1];
      coord.x1 = buff[3*i + 2];

      #ifdef DEBUG
      global_coords[n] = coord;
      n++;
      #endif

      fi_ctrl2read.write(coord);
      fi_ctrl2write.write(coord);
    }
  }
}
