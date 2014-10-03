#ifndef __BAYESIANCLUSTERING_H__
#define __BAYESIANCLUSTERING_H__

#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <iomanip>

#include "common.h"
#include "alignment.h"

using namespace std;

//Global Variables
//*******************************
#define MAXSize 20000


/*
//Golden Standard Store
//*******************************
extern int StandardZ[MAXSize];
extern int Standardk;
extern int StandardNN[MAXSize];
extern double Purity,Pu;
extern double VI,V;
extern int EditDistance;
extern double ED;
*/


//Parameters
//*******************************
typedef struct _Parameters{
        int N;
        int Step;
        int k;
        double* Distrk;
        double* Time;
//extern int* Indicator;
         double* Pi;
         double* PiPrev;
         int* Center;
         double* Sigma;
         bool* SigmaChange;
         double* Likelihood;
         double* LikelihoodRatio; //for UpdateLikelihood;
//************************************
         double** LikelihoodTemp;
//************************************
         char** SeqID;     //for Sequence ID
         char** seq;  //for reading Fasta
         char*** EigenSeqs; //for further iteration use
         float** X;         //Distance Matrix
//double X2[MAXSize][MAXSize];
         int* Z;       //missing data
         int* NN;       //NN[i-1]=# of samples in i-th cluster
         int* IterNN;   //Record # of seqs in a cluster in this round (everyone is a center seq in previous round)
         int* TempNN;
         int* Abundance; //for further iteration use
         double* Centerlikelihood;  //square sum of current center
         gsl_rng *r;
         bool* Weight;  //for Birth-Death control;
         unsigned long int WeightSum;
//*******************************

         //Result Process
//*******************************
         int Resultk;
         double* ResultPi;
         int* ResultCenter;
         double* ResultSigma;
         double* ResultLikelihood;
//*******************************
         //Hyperparameters
//*******************************
//For Pi
         double* Gamma;
//For Sigma
         double Alpha;
         double Beta;
         double g;
         double h;
         double SigmaSum;
         double Lower, Upper;
//For New-born Pi
         double Gamma1;
//For Birth-rate
         double Lamdab;
//For Stationary Time
         double Lamda;
//*******************************
//Time Inspection
//*******************************
         double TimeZ;
         double TimeCenter;
         double TimeSigma;
         double TimePi;
         double TimeBirth;
         double TimeDeath;
         double TimeLikelihood;       
}Parameters;

//Result Storage
//********************************
typedef struct _BayesianClusteringResult{
        int ClusterNumber;
        char** ClusterCenterID;
        char** ClusterCenterSeq;
        double* ClusterStandardDeviation;
        int* ClusterSize;
        int* ClusterIterationSize;
        string** ClusterMember;
        string** ClusterMemberSeq;
        bool** ClusterMemberFlag;        
}BayesianClusteringResult;

//Functions
//*******************************
void Initialize(Parameters&);
void InitRandomGenerator(Parameters&);
char** readBlockSeq(const char*, unsigned int*, int, Parameters&);
void LoadDistanceMatrix(const char*, int, Parameters&);
int GenerateInitial(int, Parameters&);
int UpdateZ(Parameters&);
int UpdatePi(Parameters&);
int UpdateCenter(int, bool, Parameters&);
int UpdateSigma(Parameters&);
int UpdateLikelihoodRatio(Parameters&);
int Birth(int,Parameters&);
int Death(int,int,Parameters&);
int ProcessSampling(int, int,Parameters&);
int ResultProcess(Parameters&);
//int CriteriaCalc(int);
double factorial(int,int);
double gaussian(double, double, Parameters&);
void NewMem(int,Parameters&);
void DeleteMem(int,int,Parameters&);
BayesianClusteringResult BayesianClustering(const char*,int, double, double, BayesianClusteringResult&, int, Parameters&);
extern "C"
{
    char** readFastASeq(char*, unsigned int*, Parameters&);
}
//*******************************



#endif//__BAYESIANCLUSTERING_H__
