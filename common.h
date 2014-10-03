#ifndef __common_H__
#define __common_H__

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <map>
#include <vector>
#include <ctime>
#include <cmath>
#include <cstring>
#include <memory>
#include <omp.h>

using namespace std;

#define MAX_LINE 50000
#define LIMIT 50000
#define LIMITB 5000
#define MIN(a, b) (((a) < (b))? (a):(b))
#define MAX(a, b) (((a) > (b))? (a):(b))
#define MAX3(a, b, c) (((a) > (b))? MAX((a), (c)):MAX((b), (c)))
#define MIN3(a, b, c) (((a) > (b))? MIN((b), (c)):MIN((a), (c)))
#define ABS(a,b) (((a) > (b))? ((a)-(b)):((b)-(a)))
#define SQUARE(a,b) (ABS((a),(b))*ABS((a),(b)))
#define DELETE(a) \
{\
    while((a)!=NULL)\
    {\
        delete[] (a);\
        (a) = NULL;\
    }\
}

#endif
