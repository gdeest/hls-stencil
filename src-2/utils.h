#ifndef UTILS_H
#define UTILS_H

int in_array(int T, int N0, int N1, int t, int x0, int x1);
void unskew_tile(int t, int x0, int x1, int *_t, int *_x0, int *_x1);
void unskew_point(int t, int x0, int x1, int *_t, int *_x0, int *_x1);

#endif
