
#ifndef __MA_H__
#define __MA_H__
#include <stdio.h>
#include "common.h"
extern "C" void maInitialization(int textureSize);

extern "C" void maDeInitialization();

extern "C" void ma(short *input,short *output,int STEP);

//#define MAMARKER -32768
#define INF 1800000000

#endif
