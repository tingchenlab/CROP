#include <iostream>
#include <getopt.h>

#include "CROP.h"
#include "Unique.h"

using namespace std;

char* program = NULL;

int main(int argc, char* argv[])
{
	int next_option;
	/* A string listing valid short options letters. */
	const char* const short_options = "i:o:b:gse:m:z:l:u:r:";
	/* An array describing valid long options. */
	const struct option long_options[] = {
		{"input", 0, NULL, 'i'},
		{"output", 0, NULL, 'o'},
		{"block", 0, NULL, 'b'},
		{"genus", 0, NULL, 'g'},
		{"species", 0, NULL, 's'},
		{"step", 0, NULL, 'e'},
		{"maxiter", 0, NULL, 'm'},
		{"maxiterblocksize", 0, NULL, 'z'},
		{"threshold lower", 0, NULL, 'l'},
		{"threshold upper", 0, NULL, 'u'},
		{"rare cut", 0, NULL, 'r'},
		{NULL, 0, NULL, 0}
	};
	program = argv[0];
    char* fname=NULL;
	char* oname=NULL;
    int TempBlockSize=50;
    bool BlockSizeFlag=false;
    int TempStep=2000;
    int TempCtrl=0;
    int TempMaxiterBlockSize=500;
    int TempMaxiter=20;
    int NumSeq=0;
    int tRareCut=2;
    float Lower=0.0;
    float Upper=0.0;
	/* Process the arguments */
	do{
		next_option = getopt_long(argc, argv, short_options, long_options, NULL);
		switch(next_option){
			case 'i': 
                fname=new char[strlen(optarg)+1];
                memset(fname,0x00,sizeof(char)*(strlen(optarg)+1));
				strcpy(fname,optarg);
				break;
			case 'o':
                oname=new char[strlen(optarg)+1];
                memset(oname,0x00,sizeof(char)*(strlen(optarg)+1));
                strcpy(oname,optarg);
                break;
            case 'b':
                 TempBlockSize=atoi(optarg);
                 BlockSizeFlag=true;
                 break;
            case 'g':
                 TempCtrl=0;
                 break;
            case 's':
                 TempCtrl=1;
                 break;
            case 'e':
                 TempStep=atoi(optarg);
                 break;
            case 'm':
                 TempMaxiter=atoi(optarg);
                 break;
            case 'z':
                 TempMaxiterBlockSize=atoi(optarg);
                 break;
            case 'u':
                 TempCtrl=2;
                 Upper=atof(optarg);
                 break;
            case 'l':
                 TempCtrl=2;
                 Lower=atof(optarg);
                 break;
            case 'r':
                 tRareCut=atoi(optarg);
                 break;
			case -1: /* Done with options */
				 break;
			default: /* Unexpected */
				cout<<"Unexpected argument was taken."<<endl;
		}
	}while(next_option != -1);
	if((TempCtrl==2)&&(Lower*Upper==0))
	{
        cout<<"invalid threshold setting (lower or upper = 0), using default -g option"<<endl;
        TempCtrl=0;
    }
	CROPParameters CParameter;   
    NumSeq=ExtractUnique(fname);
    cout<<NumSeq<<endl;
	if(NumSeq==0)
    {
        return 0;
    }
    if(BlockSizeFlag==false)  //if -b parameter is not specified
        TempBlockSize=(NumSeq/TempMaxiterBlockSize)+1;
    if(oname==NULL)           //if -o parameter is not specified
    {
        string onametemp(fname);
        string fnametemp(fname);
        fnametemp+=".unique";
        CROP(fnametemp, onametemp, TempBlockSize, TempStep, TempCtrl, Lower, Upper, TempMaxiter, TempMaxiterBlockSize, tRareCut, CParameter);
    }
    else
    {
        string onametemp(oname);
        string fnametemp(fname);
        fnametemp+=".unique";
        CROP(fnametemp, onametemp, TempBlockSize, TempStep, TempCtrl, Lower, Upper,TempMaxiter, TempMaxiterBlockSize, tRareCut, CParameter);
    }
	DELETE(fname);
	DELETE(oname);
	system("pause");
    return 1;	
}
