#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define RAND() ((double)rand() / (RAND_MAX - 1) - 0.5)
#define SCRAND(scale) ((scale)*RAND())
#define NOISE(scale) SCRAND((scale))

