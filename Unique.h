#ifndef __Unique_H__
#define __Unique_H__

#include "common.h"

using namespace std;

typedef struct _UniqueSeqs{
        string Name;
        string Members;
        string Sequences;
        int Size;        
}UniqueSeqs;

int ExtractUnique(const char*);
void StringUpper(string&);
#endif //ifndef _Unique_H
