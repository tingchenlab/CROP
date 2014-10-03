#ifndef __CROP_H__
#define __CROP_H__

#include "common.h"
#include "bayesianclustering.h"

using namespace std;

const int ReadSize = 5000;

typedef struct _TempBlockFile{
        int SeqNum;
        string* SeqID;
        string* Seq;
        double* ClusterStandardDeviation;
        int* ClusterSize;
        string** EigenSeqs;
        int* Assignment;
}TempBlockFile;

typedef struct _CROPParameters{
        int Block[20];
        int step;
        double Lower;
        double Upper;
        int MaxiterBlockSize;
        int RareCut;
        BayesianClusteringResult** BlockResult;
        int SeqNum;
        int UniqueSeqNum;
        int PrevIterk;
        int Finalk;
        int Maxiter;
        string* FinalCenterID;
        string* FinalCenterSeq;
        double* FinalSD;
        int* FinalSize;
        map<string,string> List;
        gsl_rng *r;
}CROPParameters;

void InitCROPParameters(CROPParameters&);
bool ChompFile(const char*, int, CROPParameters&);
bool ChompInterFile(const char*, int, CROPParameters&);

//************* 2011-5-4 For Diverse Data ********************
void CopyClusteringResult(const BayesianClusteringResult&, BayesianClusteringResult&, bool);
void MergingClusteringResult(BayesianClusteringResult&, BayesianClusteringResult&, BayesianClusteringResult&);
void ClearClusteringResult(BayesianClusteringResult&);
void MappingRareToAbundant(string, BayesianClusteringResult&, BayesianClusteringResult&, double);
void LoadBlockFile(string, TempBlockFile&);
void TempBlockFileDelete(TempBlockFile&);
void SummarizingResult(string,BayesianClusteringResult&, int, CROPParameters&);

//void MergeResult(int);
int LoadingList(const char*, CROPParameters&);
void MergeUsingBC(const char*, string, CROPParameters&);
void InterResultRead(const char*, int, CROPParameters&);
void InterResultReadForBC(string, int, CROPParameters&); //if cluster size < Rarecut, it will be stored separately
int CROP(string, string, int, int, int, float, float, int ,int, int, CROPParameters&);
void DeleteTmpFiles(string, CROPParameters&);


#endif //ifndef _CROP_H
