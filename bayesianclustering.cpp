#include "bayesianclustering.h"

void Initialize(Parameters &Parameter)
{
//Global Variables
//*******************************
    Parameter.k=1;
    Parameter.Distrk=NULL;
    Parameter.Time=NULL;
    Parameter.Pi=NULL;
    Parameter.PiPrev=NULL;
    Parameter.Center=NULL;
    Parameter.Abundance=NULL;
    Parameter.Sigma=NULL;
    Parameter.SigmaChange=NULL;
    Parameter.Likelihood=NULL;
    Parameter.LikelihoodRatio=NULL; //for UpdateLikelihood;
//************************************
    Parameter.LikelihoodTemp=NULL;
//************************************
    Parameter.SeqID=NULL;     //for Sequence ID
    Parameter.seq=NULL;  //for reading Fasta
    Parameter.EigenSeqs=NULL; //for further iteration use
    Parameter.X=NULL;         //Distance Matrix
//double X2[MAXSize][MAXSize];
    Parameter.Z=NULL;       //missing data
    Parameter.NN=NULL;       //NN[i-1]=# of samples in i-th cluster
    Parameter.IterNN=NULL;
    Parameter.TempNN=NULL;
    Parameter.Centerlikelihood=NULL;  //square sum of current center
    Parameter.r=NULL;
    Parameter.Weight=NULL;  //for Birth-Death control;
    Parameter.WeightSum=0;
//*******************************

//Result Process
//*******************************
    Parameter.Resultk=0;
    Parameter.ResultPi=NULL;
    Parameter.ResultCenter=NULL;
    Parameter.ResultSigma=NULL;
    Parameter.ResultLikelihood=NULL;
//*******************************

//Hyperparameters
//*******************************
//For Pi
    Parameter.Gamma=NULL;
//For Sigma
    Parameter.Alpha=2;
    Parameter.Beta=3;
    Parameter.g=0.2;
    Parameter.h=0;
    Parameter.SigmaSum=0;
    Parameter.Lower=0; 
    Parameter.Upper=0;
//For New-born Pi
    Parameter.Gamma1=1;
//For Birth-rate
    Parameter.Lamdab=1.0;
//For Stationary Time
    Parameter.Lamda=2;
//Time Inspection
//*******************************
    Parameter.TimeZ=0;
    Parameter.TimeCenter=0;
    Parameter.TimeSigma=0;
    Parameter.TimePi=0;
    Parameter.TimeBirth=0;
    Parameter.TimeDeath=0;
    Parameter.TimeLikelihood=0;
    return;
}

void InitRandomGenerator(Parameters &Parameter)
{
     Parameter.r = gsl_rng_alloc (gsl_rng_taus);
     gsl_rng_set(Parameter.r,(unsigned int)time((time_t *)NULL));
     return;   
}

char** readBlockSeq(const char *fname, unsigned int *num_seq, int iter, Parameters &Parameter)
{
	char dummy[MAX_LINE], **seq = NULL;
	unsigned int len = 0, lentotal = 0;
	*num_seq=0;
	fstream BlockSeq;
	BlockSeq.open(fname);
	if(BlockSeq.fail())
	    cout<<"Cannot open block sequence file"<<endl;
    while(!BlockSeq.eof())
    {
        BlockSeq.getline(dummy, MAX_LINE);
        if(dummy[0]=='>')
            (*num_seq)++;
    }
    if(*num_seq==0)
    {
        cout<<"Not A Valid Fasta File!"<<endl;
        return NULL;
    }    
    BlockSeq.clear();   
    BlockSeq.seekg(0,ios::beg);
	//cout<<"num_seq = "<<*num_seq<<endl;

    NewMem(*num_seq,Parameter);
	MALLOC(Parameter.seq, char*, *num_seq);
	MALLOC(Parameter.SeqID, char*, *num_seq);
	Parameter.EigenSeqs=new char** [*num_seq];
    if(iter!=0)
    {
    for(int i=0;i<*num_seq; i++)
    {
        //**********Read Center Seq ID***************  
        BlockSeq.getline(dummy, MAX_LINE);
        len=strlen(dummy);
        len+=(dummy[len-1]=='\n')? -1:0;
        MALLOC(Parameter.SeqID[i], char, len);
        memcpy(Parameter.SeqID[i],dummy+1, sizeof(char)*(len-1)); //erase ">"
        //**************************************
        
        //**********Read Center Seq********************
        BlockSeq.getline(dummy, MAX_LINE);
        len=strlen(dummy);
        len+=(dummy[len-1]=='\n')? -1:0;
        MALLOC(Parameter.seq[i], char, len+1);
        memcpy(Parameter.seq[i],dummy, sizeof(char)*len);
        //**************************************
        
        BlockSeq>>Parameter.Abundance[i];
        BlockSeq>>Parameter.Sigma[i];
        BlockSeq.getline(dummy, MAX_LINE);
        int NumberOfSeqs=MIN(20,Parameter.Abundance[i]);
        if(NumberOfSeqs==1)
        {
            Parameter.EigenSeqs[i]=new char*[NumberOfSeqs+1];
            Parameter.EigenSeqs[i][NumberOfSeqs]=NULL;
        }
        else
            Parameter.EigenSeqs[i]=new char*[NumberOfSeqs];
        //***********Read Eigen Seqs*************
	    for(int j=0;j<NumberOfSeqs;j++)
	    {
            BlockSeq.getline(dummy, MAX_LINE);
            len=strlen(dummy);
            len+=(dummy[len-1]=='\n')? -1:0;
            Parameter.EigenSeqs[i][j]=new char[len+1];
            memset(Parameter.EigenSeqs[i][j], 0x00, sizeof(char)*(len+1));
            memcpy(Parameter.EigenSeqs[i][j], dummy, sizeof(char)*len);
		}
		//*************************************************
	}
    }//if(iter!=0)
    else
    {
        for(int i=0;i<*num_seq; i++)
        {
        //**********Read Center Seq ID***************  
        BlockSeq.getline(dummy, MAX_LINE);
        len=strlen(dummy);
        len+=(dummy[len-1]=='\n')? -1:0;
        MALLOC(Parameter.SeqID[i], char, len);
        memcpy(Parameter.SeqID[i],dummy+1, sizeof(char)*(len-1));
        //**************************************
        
        //**********Read Center Seq********************
        BlockSeq.getline(dummy, MAX_LINE);
        len=strlen(dummy);
        len+=(dummy[len-1]=='\n')? -1:0;
        MALLOC(Parameter.seq[i], char, len+1);
        memcpy(Parameter.seq[i],dummy, sizeof(char)*len);
        //**************************************        
        BlockSeq>>Parameter.Abundance[i];
        BlockSeq.getline(dummy, MAX_LINE); //to pass this line
        }
    } 
    BlockSeq.close();
	return Parameter.seq;
}

void LoadDistanceMatrix(const char* fname, int iter, Parameters &Parameter)
{
    //fstream tempout1,tempout2;
    float maxtemp=0.0, mintemp=100.0;
    int matrix = EDNAFULL;
	unsigned int seq_id = 0, tar_id = 0;
	alignment (*alignTool) (const char *, const char *, alignment, int, int) = NULL, align = {0, 0, 0};
	unsigned int UnN=(unsigned int)Parameter.N;    
    if(iter==0)
    {
    Parameter.seq = readBlockSeq(fname, &UnN, iter, Parameter);
    //Parameter.seq = readFastASeq(ffname, &UnN, Parameter);
    Parameter.N=(int)UnN;
    //NewMem(Parameter.N,Parameter);
    //for(int i=0;i<Parameter.N;i++)
        //Parameter.Abundance[i]=1; 
    alignTool = NW_alignmentEx;
	for(seq_id = 0; seq_id < Parameter.N-1; seq_id++){
        Parameter.X[seq_id][seq_id]=0.0;
		for(tar_id = seq_id+1; tar_id < Parameter.N; tar_id++){
            int tempint=seq_id*Parameter.N-(seq_id)*(seq_id-1)/2;
            //cout<<"\rLoading Distance Matrix "<<(tempint+(tar_id-seq_id)+1)/(Parameter.N*(Parameter.N+1)/2-1)<<"%";
			align = alignTool(Parameter.seq[seq_id], Parameter.seq[tar_id], align, matrix, 1);
			Parameter.X[seq_id][tar_id]=(float)(1.0-(double)align.matches/(double)SHORT_ALIGN(align))*100.0;
		}
	}
	//cout<<endl;
	Parameter.X[Parameter.N-1][Parameter.N-1]=0.0;
    for(seq_id=0; seq_id < Parameter.N; seq_id++)
        for(tar_id = 0; tar_id < seq_id; tar_id++)
            Parameter.X[seq_id][tar_id]=Parameter.X[tar_id][seq_id];
    }//if iter==0

    else if(iter==20)
    {   
    Parameter.seq = readBlockSeq(fname, &UnN, iter, Parameter);
    Parameter.N=(int)UnN;
    if(Parameter.N>LIMIT||Parameter.N==0)     //2011-1-19  When merging using BC, if size too large, skip
        return;
    //k=N;   
    alignTool = NW_alignmentEx;
	for(seq_id = 0; seq_id < Parameter.N; seq_id++){
        Parameter.X[seq_id][seq_id]=Parameter.Sigma[seq_id];
        #pragma omp parallel for
	    for(tar_id = seq_id+1; tar_id < Parameter.N; tar_id++){
            int tempint=seq_id*Parameter.N-(seq_id)*(seq_id-1)/2;
            //#pragma omp critical
            //{
                //cout<<"\rLoading Distance Matrix "<<(tempint+(tar_id-seq_id)+1)/(Parameter.N*(Parameter.N+1)/2-1)<<"%";
            //}
			align = alignTool(Parameter.seq[seq_id], Parameter.seq[tar_id], align, matrix, 1);
			Parameter.X[seq_id][tar_id]=(float)(1.0-(double)align.matches/(double)SHORT_ALIGN(align))*100.0;
            Parameter.X[tar_id][seq_id]=Parameter.X[seq_id][tar_id];
			if((Parameter.X[seq_id][tar_id]<15)&&(Parameter.X[seq_id][tar_id]>(float)Parameter.Lower/2))
			{
                #pragma omp parallel for                                                                                        
                for(int i=0;i<MIN(20,Parameter.Abundance[seq_id]);i++)
                {
                    align = alignTool(Parameter.seq[tar_id], Parameter.EigenSeqs[seq_id][i], align, matrix, 1);
                    Parameter.X[seq_id][tar_id]+=(float)(1.0-(double)align.matches/(double)SHORT_ALIGN(align))*100.0;
                }
                Parameter.X[seq_id][tar_id]/=(MIN(20,Parameter.Abundance[seq_id])+1);
                #pragma omp parallel for
                for(int i=0;i<MIN(20,Parameter.Abundance[tar_id]);i++)
                {
                    align = alignTool(Parameter.seq[seq_id], Parameter.EigenSeqs[tar_id][i], align, matrix, 1);
                    Parameter.X[tar_id][seq_id]+=(float)(1.0-(double)align.matches/(double)SHORT_ALIGN(align))*100.0;
                }
                Parameter.X[tar_id][seq_id]/=(MIN(20,Parameter.Abundance[tar_id])+1);
            }
		}//tar_id
	}//seq_id
	cout<<endl;
	Parameter.X[Parameter.N-1][Parameter.N-1]=Parameter.Sigma[Parameter.N-1];
	cout<<"Loading Distance Complete"<<endl;
    }//else if
        
    else
    {   
    Parameter.seq = readBlockSeq(fname, &UnN, iter, Parameter);
    Parameter.N=(int)UnN;
    //k=N;   
    alignTool = NW_alignmentEx;
	for(seq_id = 0; seq_id < Parameter.N; seq_id++){
        Parameter.X[seq_id][seq_id]=Parameter.Sigma[seq_id];
	    for(tar_id = seq_id+1; tar_id < Parameter.N; tar_id++){
            int tempint=seq_id*Parameter.N-(seq_id)*(seq_id-1)/2;
            //cout<<"\rLoading Distance Matrix "<<(tempint+(tar_id-seq_id)+1)/(Parameter.N*(Parameter.N+1)/2-1)<<"%";
			align = alignTool(Parameter.seq[seq_id], Parameter.seq[tar_id], align, matrix, 1);
			Parameter.X[seq_id][tar_id]=(float)(1.0-(double)align.matches/(double)SHORT_ALIGN(align))*100.0;
            Parameter.X[tar_id][seq_id]=Parameter.X[seq_id][tar_id];
			if((Parameter.X[seq_id][tar_id]<15)&&(Parameter.X[seq_id][tar_id]>(float)Parameter.Lower/2))
			{
                for(int i=0;i<MIN(20,Parameter.Abundance[seq_id]);i++)
                {
                    align = alignTool(Parameter.seq[tar_id], Parameter.EigenSeqs[seq_id][i], align, matrix, 1);
                    Parameter.X[seq_id][tar_id]+=(float)(1.0-(double)align.matches/(double)SHORT_ALIGN(align))*100.0;
                }
                Parameter.X[seq_id][tar_id]/=(MIN(20,Parameter.Abundance[seq_id])+1);
                for(int i=0;i<MIN(20,Parameter.Abundance[tar_id]);i++)
                {
                    align = alignTool(Parameter.seq[seq_id], Parameter.EigenSeqs[tar_id][i], align, matrix, 1);
                    Parameter.X[tar_id][seq_id]+=(float)(1.0-(double)align.matches/(double)SHORT_ALIGN(align))*100.0;
                }
                Parameter.X[tar_id][seq_id]/=(MIN(20,Parameter.Abundance[tar_id])+1);
            }
		}//tar_id
	}//seq_id
	//cout<<endl;
	Parameter.X[Parameter.N-1][Parameter.N-1]=Parameter.Sigma[Parameter.N-1];
    }//else
    for(int i=0;i<Parameter.N;i++)
    {   
        for(int j=0;j<Parameter.N;j++)
        {
            if(Parameter.X[i][j]>=maxtemp)
                 maxtemp=Parameter.X[i][j];
            if((Parameter.X[i][j]<=mintemp)&&(Parameter.X[i][j]!=0))
                 mintemp=Parameter.X[i][j];
        }
    }
    Parameter.h=(double)100*Parameter.g/(maxtemp-mintemp)/(maxtemp-mintemp)/Parameter.Alpha;
    //if(iter==20)
    //{
        //tempout1.close();
        //tempout2.close();
    //}
    //cout<<"Loading Distance Complete"<<endl;
    return;    
}

int GenerateInitial(int iter, Parameters &Parameter)
{
    //if(iter==0)
    //{    
    //Init Weight & WeightSum    
    for(int i=0;i<Parameter.N;i++)
            Parameter.Weight[i]=true; //false
    memset(Parameter.PiPrev,0,sizeof(double)*Parameter.N);
    memset(Parameter.Pi,0,sizeof(double)*Parameter.N);
    Parameter.WeightSum=Parameter.N; //0
    //Init Pi
    for(int i=0;i<Parameter.N;i++)
            Parameter.Gamma[i]=1.0;
            //Pi[i]=0.2;
    gsl_ran_dirichlet(Parameter.r, Parameter.k, Parameter.Gamma, Parameter.Pi);
    memcpy(Parameter.PiPrev, Parameter.Pi, sizeof(double)*Parameter.k);
    memcpy(Parameter.ResultPi, Parameter.Pi, sizeof(double)*Parameter.k);
    //Init Center
    for(int i=0;i<Parameter.k;i++)
    {
            //Center[i]=(int)gsl_ran_flat(r,100*i,100*i+99);
            Parameter.Center[i]=(int)gsl_ran_flat(Parameter.r,0,Parameter.N);
            Parameter.ResultCenter[i]=Parameter.Center[i];
            Parameter.Weight[Parameter.Center[i]]=false;
            Parameter.WeightSum--;
    }
    //Init Sigma
    Parameter.SigmaSum=0;
    Parameter.Beta=gsl_ran_gamma(Parameter.r,Parameter.g,(double)1/Parameter.h);
    for(int i=0;i<Parameter.k;i++)
    {
            for(int j=0;j<10;j++)
            {
                Parameter.Sigma[i]=sqrt(1/gsl_ran_gamma(Parameter.r, Parameter.Alpha, (double)1/Parameter.Beta));
                if((Parameter.Sigma[i]<=Parameter.Upper)&&(Parameter.Sigma[i]>=Parameter.Lower))
                    break;
            }
            while((Parameter.Sigma[i]>Parameter.Upper)||(Parameter.Sigma[i]<Parameter.Lower))
            {
                Parameter.Sigma[i]=sqrt(1/gsl_ran_gamma(Parameter.r, Parameter.Alpha+1, (double)1/(Parameter.Upper+Parameter.Lower)));
            }
            Parameter.ResultSigma[i]=Parameter.Sigma[i];
            Parameter.SigmaSum+=(double)1/Parameter.Sigma[i]/Parameter.Sigma[i];
    }
    memset(Parameter.SigmaChange,0,sizeof(bool)*Parameter.k);
    //Init LikelihoodRatio
    UpdateLikelihoodRatio(Parameter);
    for(int i=0;i<Parameter.N;i++)
        Parameter.ResultLikelihood[i]=Parameter.Likelihood[i];
    Parameter.Resultk=Parameter.k;
    //}//if(iter==0)
    /*
    else
    {
    //Init Weight & WeightSum    
    for(int i=0;i<N;i++)
            Weight[i]=false; //false
    WeightSum=0; //0
    //Init Pi
    for(int i=0;i<N;i++)
            Gamma[i]=1.0;
    double* GammaTemp=new double[k+1];
    memset(GammaTemp,0,sizeof(double)*(k+1));
    for(int i=0;i<k;i++)
        GammaTemp[i]=Gamma[i]+Abundance[i];
    memset(PiPrev,0,sizeof(double)*N);
    memset(Pi,0,sizeof(double)*N);
    gsl_ran_dirichlet(r, k, GammaTemp, Pi);      
    delete[] GammaTemp;
    memcpy(PiPrev, Pi, sizeof(double)*k);
    memcpy(ResultPi, Pi, sizeof(double)*k);
    //Init Center
    for(int i=0;i<k;i++)
    {
            Center[i]=i;
            ResultCenter[i]=Center[i];
    }
    //Init Sigma
    SigmaSum=0;
    Beta=gsl_ran_gamma(r,g,(double)1/h);
    for(int i=0;i<k;i++)
    {
            ResultSigma[i]=Sigma[i];
            SigmaSum+=(double)1/Sigma[i]/Sigma[i];
    }
    memset(SigmaChange,0,sizeof(bool)*k);
    //Init LikelihoodRatio
    UpdateLikelihoodRatio();
    for(int i=0;i<N;i++)
        ResultLikelihood[i]=Likelihood[i];
    delete[] GammaTemp;
    Resultk=k;        
    }//else
    */
    return 1;   
}

int UpdateZ(Parameters &Parameter)
{
     time_t start;
     start=time(NULL);
     memset(Parameter.NN,0,sizeof(int)*Parameter.k);
     double z,temp=0;
     bool flag;
     //******************************************
     for(int i=0;i<Parameter.N;i++)
     {
          flag=false;
          for(int j=0;j<Parameter.k;j++)
          {
                  if(i==Parameter.Center[j])
                  {
                       Parameter.Z[i]=j;
                       Parameter.NN[j]+=Parameter.Abundance[i];
                       flag=true;
                       break;
                  }
                       
          }
          if(flag==false)
          {   
             z=gsl_ran_flat(Parameter.r, 0, Parameter.Likelihood[i]);
             temp=0;
             for(int j=0;j<Parameter.k;j++)
             {
                     temp+=Parameter.LikelihoodTemp[i][j];
                     if(temp>=z)
                     {
                         Parameter.Z[i]=j;
                         Parameter.NN[j]+=Parameter.Abundance[i];
                         break;
                     }
             }
          }
     }
     /*
     if(k==1)
     {
         for(int i=0;i<N;i++)
             if(Z[i]!=0)
                 cout<<"Z["<<i<<"]";
     }
     */
     //*******************************************
     time_t stop;
     stop=time(NULL);
     Parameter.TimeZ+=stop-start;
     //cout<<"Z";
     return 1;
}

int UpdatePi(Parameters &Parameter)
{
     time_t start;
     start=time(NULL);
       int d=0;
       //Sampling Method
       memcpy(Parameter.PiPrev,Parameter.Pi,sizeof(double)*Parameter.k);
       double* GammaTemp=new double[Parameter.k+1];
       memset(GammaTemp,0,sizeof(double)*(Parameter.k+1));
       for(int i=0;i<Parameter.k;i++)
               GammaTemp[i]=Parameter.Gamma[i]+Parameter.NN[i];
       gsl_ran_dirichlet(Parameter.r, Parameter.k, GammaTemp, Parameter.Pi);      
       //Posterior Mean Method        
       /*
       for(int i=0;i<k;i++)
               d+=Gamma[i]+NN[i];
       for(int i=0;i<k;i++)
               Pi[i]=(double)(Gamma[i]+NN[i])/d;
       */
       DELETE(GammaTemp);
     time_t stop;
     stop=time(NULL);
     Parameter.TimePi+=stop-start;  
     //cout<<"P";            
     return 1;
}

int UpdateCenter(int m, bool f, Parameters &Parameter)
{
     time_t start;
     start=time(NULL);
     double *p=NULL,*rand=NULL,*temp=NULL;
     bool *flag=NULL;
     if(m%MAX(10,Parameter.N/100)==0) //2011-4-13 less update center if N is large & only update for those with size>10
     {
     memset(Parameter.Centerlikelihood,0,sizeof(double)*Parameter.k);
     for(int i=0;i<Parameter.N;i++)
     {
         double x=Parameter.X[i][Parameter.Center[Parameter.Z[i]]];
         Parameter.Centerlikelihood[Parameter.Z[i]]+=((x*x/2)*Parameter.Abundance[i]);
     }
     double* TempCenterLikelihood=new double[Parameter.N];      //for UpdateCenter;
     memset(TempCenterLikelihood, 0,sizeof(double)*Parameter.N);
     if(Parameter.k>1)
     {
         p=new double[Parameter.k];
         rand=new double[Parameter.k];
         temp=new double[Parameter.k];
         flag=new bool[Parameter.k];
         memset(rand,0,sizeof(double)*Parameter.k);
         memset(flag,0,sizeof(bool)*Parameter.k);
         memset(p,0,sizeof(double)*Parameter.k);
         memset(temp,0,sizeof(double)*Parameter.k);
     }
     else
     {
         p=new double[Parameter.k+1];
         rand=new double[Parameter.k+1];
         temp=new double[Parameter.k+1];
         flag=new bool[Parameter.k+1];
         memset(rand,0,sizeof(double)*(Parameter.k+1));
         memset(flag,0,sizeof(bool)*(Parameter.k+1));
         memset(p,0,sizeof(double)*(Parameter.k+1));
         memset(temp,0,sizeof(double)*(Parameter.k+1));
     }
     for(int i=0;i<Parameter.N;i++)
     {
             for(int j=0;j<Parameter.N;j++)
             {
                  if(Parameter.Z[j]==Parameter.Z[i]&&Parameter.NN[Parameter.Z[i]]>10)
                  {
                      double x=Parameter.X[j][i];
                      TempCenterLikelihood[i]+=((x*x/2)*Parameter.Abundance[j]);
                  }
             }
             if(Parameter.NN[Parameter.Z[i]]>10)
             {
                 TempCenterLikelihood[i]-=Parameter.Centerlikelihood[Parameter.Z[i]];
                 TempCenterLikelihood[i]=gsl_expm1(-TempCenterLikelihood[i]/Parameter.Sigma[Parameter.Z[i]]/Parameter.Sigma[Parameter.Z[i]])+1;
                 p[Parameter.Z[i]]+=MIN(10e300,TempCenterLikelihood[i]);
             }
     }
     for(int j=0;j<Parameter.k;j++)
             if(Parameter.NN[j]>10)
                  rand[j]=gsl_ran_flat(Parameter.r, 0, p[j]);
     for(int i=0;i<Parameter.N;i++)
     {
         if(Parameter.NN[Parameter.Z[i]]>10)
         {
             temp[Parameter.Z[i]]+=MIN(10e300,TempCenterLikelihood[i]);
             if((temp[Parameter.Z[i]]>=rand[Parameter.Z[i]])&&(flag[Parameter.Z[i]]==0))
             {
                     Parameter.Weight[Parameter.Center[Parameter.Z[i]]]=true;                                     
                     Parameter.Center[Parameter.Z[i]]=i;
                     Parameter.Weight[i]=false;
                     flag[Parameter.Z[i]]=1;
             }
         }
     }
    DELETE(temp);
    DELETE(rand);
    DELETE(p);
    DELETE(TempCenterLikelihood);
    DELETE(flag);
    time_t stop;
    stop=time(NULL);
    Parameter.TimeCenter+=stop-start;
    return 1;
    }
    else if(f==1)
    {
     Parameter.Centerlikelihood[Parameter.k-1]=0;
     for(int i=0;i<Parameter.N;i++)
         if(Parameter.Z[i]==Parameter.k-1)
         {
             double x=Parameter.X[i][Parameter.Center[Parameter.k-1]];
             Parameter.Centerlikelihood[Parameter.k-1]+=((x*x/2)*Parameter.Abundance[i]);
         }
     double* TempCenterLikelihood=new double[Parameter.N];      //for UpdateCenter;
     memset(TempCenterLikelihood, 0,sizeof(double)*Parameter.N);
     if(Parameter.k>1)
     {
         p=new double[Parameter.k];
         rand=new double[Parameter.k];
         temp=new double[Parameter.k];
         flag=new bool[Parameter.k];
         memset(rand,0,sizeof(double)*Parameter.k);
         memset(flag,0,sizeof(bool)*Parameter.k);
         memset(p,0,sizeof(double)*Parameter.k);
         memset(temp,0,sizeof(double)*Parameter.k);
     }
     else
     {
         p=new double[Parameter.k+1];
         rand=new double[Parameter.k+1];
         temp=new double[Parameter.k+1];
         flag=new bool[Parameter.k+1];
         memset(rand,0,sizeof(double)*(Parameter.k+1));
         memset(flag,0,sizeof(bool)*(Parameter.k+1));
         memset(p,0,sizeof(double)*(Parameter.k+1));
         memset(temp,0,sizeof(double)*(Parameter.k+1));
     }
     for(int i=0;i<Parameter.N;i++)
     {
         if(Parameter.Z[i]==Parameter.k-1)
         {
             for(int j=0;j<Parameter.N;j++)
             {
                  if(Parameter.Z[j]==Parameter.k-1)
                  {
                      double x=Parameter.X[j][i];
                      TempCenterLikelihood[i]+=((x*x/2)*Parameter.Abundance[j]);
                  }
             }
             TempCenterLikelihood[i]-=Parameter.Centerlikelihood[Parameter.k-1];
             TempCenterLikelihood[i]=gsl_expm1(-TempCenterLikelihood[i]/Parameter.Sigma[Parameter.k-1]/Parameter.Sigma[Parameter.k-1])+1;
             p[Parameter.k-1]+=TempCenterLikelihood[i];
         }
     }
     rand[Parameter.k-1]=gsl_ran_flat(Parameter.r, 0, p[Parameter.k-1]);
     for(int i=0;i<Parameter.N;i++)
     {
         if(Parameter.Z[i]==Parameter.k-1)
         {
             temp[Parameter.k-1]+=TempCenterLikelihood[i];
             if((temp[Parameter.k-1]>=rand[Parameter.k-1])&&(flag[Parameter.k-1]==0))
             {
                     Parameter.Weight[Parameter.Center[Parameter.k-1]]=true;                                     
                     Parameter.Center[Parameter.k-1]=i;
                     Parameter.Weight[i]=false;
                     flag[Parameter.k-1]=1;
             }
        }
    }
    DELETE(temp);
    DELETE(rand);
    DELETE(p);
    DELETE(TempCenterLikelihood);
    DELETE(flag);
    }
    //cout<<"C";
    time_t stop;
    stop=time(NULL);
    Parameter.TimeCenter+=stop-start;
    return 1;      
}

int UpdateSigma(Parameters &Parameter)
{
     time_t start;
     start=time(NULL);
     double temp;
    Parameter.Beta=gsl_ran_gamma(Parameter.r, (double)Parameter.g+Parameter.k*Parameter.Alpha, (double)1/(Parameter.h+Parameter.SigmaSum));
    Parameter.SigmaSum=0;
    memset(Parameter.SigmaChange,0,sizeof(bool)*Parameter.k);
    //double temp=0;
    //for(int i=0;i<k;i++)
        //temp+=Centerlikelihood[i]; 
    for(int i=0;i<Parameter.k;i++)
    {
            //Sigma[i]=sqrt(1/gsl_ran_gamma(r, Alpha+(double)NN[i]/2, (double)1/(Beta+temp)));   
            temp=sqrt(1/gsl_ran_gamma(Parameter.r, Parameter.Alpha+(double)Parameter.NN[i]/2, (double)1/(Parameter.Beta+Parameter.Centerlikelihood[i])));
            /*   
            if(Sigma[i]>Threshold)
                Sigma[i]=sqrt(1/gsl_ran_gamma(r, Alpha*10+1, (double)1/Threshold/10));
            if(Sigma[i]<1)
                Sigma[i]=sqrt(1/gsl_ran_gamma(r, Alpha*10+1, (double)1/Threshold/10));
            */
            /*
            while((Sigma[i]>Threshold)||(Sigma[i]<1))
                Sigma[i]=sqrt(1/gsl_ran_gamma(r, Alpha+1, (double)1/Threshold));
            */
           while((temp>Parameter.Upper)||(temp<Parameter.Lower))
                temp=sqrt(1/gsl_ran_gamma(Parameter.r, Parameter.Alpha+1, (double)1/(Parameter.Upper+Parameter.Lower)));
           double ratio=temp/Parameter.Sigma[i];
           if((ratio<1.2)&&(ratio>0.8))
                Parameter.SigmaChange[i]=true;
           else
               Parameter.Sigma[i]=temp;
           Parameter.SigmaSum+=1/Parameter.Sigma[i]/Parameter.Sigma[i];    
    }
    //cout<<"V";
     time_t stop;
     stop=time(NULL);
     Parameter.TimeSigma+=stop-start;
    return 1;
}

int Birth(int m,Parameters &Parameter)
{
     time_t start;
     start=time(NULL);
    //Init new state
    Parameter.Pi[Parameter.k-1]=gsl_ran_beta(Parameter.r, Parameter.Gamma1, (double)Parameter.k*Parameter.Gamma1);
    //**********************choose intial sigma*****************************
    for(int j=0;j<10;j++)
    {
        Parameter.Sigma[Parameter.k-1]=sqrt(1/gsl_ran_gamma(Parameter.r, Parameter.Alpha, (double)1/Parameter.Beta));
        if((Parameter.Sigma[Parameter.k-1]<=Parameter.Upper)&&(Parameter.Sigma[Parameter.k-1]>=Parameter.Lower))
            break;
    }
    while((Parameter.Sigma[Parameter.k-1]>Parameter.Upper)||(Parameter.Sigma[Parameter.k-1]<Parameter.Lower))
    {
        Parameter.Sigma[Parameter.k-1]=sqrt(1/gsl_ran_gamma(Parameter.r, Parameter.Alpha+1, (double)2/(Parameter.Upper+Parameter.Lower)));
    }
    //************************************************************************ 
    for(int i=0;i<Parameter.k-1;i++)
    {
        Parameter.Pi[i]*=(1-Parameter.Pi[Parameter.k-1]);
        for(int j=0;j<Parameter.N;j++)
            Parameter.LikelihoodTemp[j][i]*=(1-Parameter.Pi[Parameter.k-1]);
    } 
    for(int i=0;i<Parameter.N;i++)
    {
        Parameter.LikelihoodTemp[i][Parameter.k-1]=Parameter.Pi[Parameter.k-1]*gaussian(Parameter.X[i][Parameter.Center[Parameter.k-1]],Parameter.Sigma[Parameter.k-1], Parameter);
        Parameter.Likelihood[i]=Parameter.Likelihood[i]*(1-Parameter.Pi[Parameter.k-1])+Parameter.LikelihoodTemp[i][Parameter.k-1];
    }
    //Update whole parameters
    //cout<<"Birth "<<Center[k-1]<<endl;
     time_t stop;
     stop=time(NULL);
     Parameter.TimeBirth+=stop-start;
    //cout<<"B";
    UpdateZ(Parameter);
    UpdatePi(Parameter);
    UpdateCenter(m,1,Parameter);
    UpdateSigma(Parameter);
    UpdateLikelihoodRatio(Parameter);    
    return 1;       
}

int Death(int j,int m,Parameters &Parameter)
{
     time_t start;
     start=time(NULL);
    //Init new state
    for(int i=0;i<Parameter.N;i++)
            Parameter.Likelihood[i]=(Parameter.Likelihood[i]-Parameter.LikelihoodTemp[i][j])/(1-Parameter.Pi[j]);
    for(int i=0;i<Parameter.k;i++)
    {
            if(i!=j)
                Parameter.Pi[i]/=(1-Parameter.Pi[j]);
            for(int l=0;l<Parameter.N;l++)
            {
                if(i<j)
                    Parameter.LikelihoodTemp[l][i]/=(1-Parameter.Pi[j]);
                else
                   if(i<Parameter.k-1)
                    Parameter.LikelihoodTemp[l][i]=Parameter.LikelihoodTemp[l][i+1]/(1-Parameter.Pi[j]);
            }
    }
    Parameter.Weight[Parameter.Center[j]]=true;
    Parameter.WeightSum++;
    Parameter.k--;
    for(int i=j;i<Parameter.k;i++)
    {
            Parameter.Pi[i]=Parameter.Pi[i+1];
            Parameter.Center[i]=Parameter.Center[i+1];
            Parameter.Sigma[i]=Parameter.Sigma[i+1];
    }
     time_t stop;
     stop=time(NULL);
     Parameter.TimeDeath+=stop-start;
    //Update whole Parameters
    //cout<<"Death "<<j<<endl;
    //cout<<"D";
    UpdateZ(Parameter);
    UpdatePi(Parameter);
    UpdateCenter(m,0,Parameter);
    UpdateSigma(Parameter); 
    UpdateLikelihoodRatio(Parameter);    
    return 1;    
}

int ProcessSampling(int m, int iter,Parameters &Parameter)
{
    bool flag=false;
    double temp=0.0,z,Ti,L=0.0;
    flag=false;
    if(Parameter.k==1)
    {
           if(m==0)
               Ti=0;
           else
               Ti=Parameter.Time[m-1];
           Parameter.Time[m]=gsl_ran_exponential(Parameter.r,(double)1/Parameter.Lamdab);
           Parameter.Distrk[Parameter.k]+=Parameter.Time[m];
           for(int i=0;i<Parameter.N;i++)
               if(Parameter.Weight[i]==true)
               {         
                   //if(Likelihood[i]>=10.0e299||Likelihood[i]<=10.0e-300)
                       //cout<<Likelihood[i]<<" ";
                   temp+=1.0/MAX(10e-300,Parameter.Likelihood[i]);
                   //if(Likelihood[i]>=10.0e299||Likelihood[i]<=10.0e-300)
                       //cout<<temp<<" ";
               }
           //cout<<temp;      
           z=gsl_ran_flat(Parameter.r, 0, temp);
           temp=0;
           for(int i=0;i<Parameter.N;i++)
           {
               if(Parameter.Weight[i]==true)
               {
                   temp+=1.0/MAX(10.0e-300,Parameter.Likelihood[i]); 
                   if(temp>=z)
                   {
                       Parameter.Center[Parameter.k]=i;
                       Parameter.Weight[i]=false;
                       Parameter.WeightSum--;
                       Parameter.k++;
                       flag=true;
                       break;
                   }
               }
           }
           //cout<<"[1,"<<k<<"]";
           Birth(m,Parameter);
           return 1;
    }
    else
    {
        if(m==0)
            Ti=0;
        else
           Ti=Parameter.Time[m-1];
        for(int i=0;i<Parameter.k;i++)
                temp+=Parameter.LikelihoodRatio[i];
        if(Parameter.k<Parameter.N)
                temp+=Parameter.Lamdab;//(double)k;
        Parameter.Time[m]=gsl_ran_exponential(Parameter.r,(double)1/temp);
        Parameter.Distrk[Parameter.k]+=Parameter.Time[m];
        z=gsl_ran_flat(Parameter.r, 0, temp);
        temp=0;
        for(int i=0;i<Parameter.k;i++)
        {
                temp+=Parameter.LikelihoodRatio[i];
                if(temp>=z)
                {
                           //cout<<"2";
                           Death(i,m,Parameter);
                           flag=true;
                           break; 
                }
        }
        if(flag==false)
        {
           temp=0;
       	       //Randomly choose a new group center from data.
           for(int i=0;i<Parameter.N;i++)
               if(Parameter.Weight[i]==true)
                   temp+=1.0/MAX(10.0e-300,Parameter.Likelihood[i]);        
           z=gsl_ran_flat(Parameter.r, 0, temp);
           temp=0;
           for(int i=0;i<Parameter.N;i++)
           {
               if(Parameter.Weight[i]==true)
               {
                       temp+=1.0/MAX(10.0e-300,Parameter.Likelihood[i]);
                   if(temp>=z)
                   {
                       Parameter.Center[Parameter.k]=i;
                       Parameter.Weight[i]=false;
                       Parameter.WeightSum--;
                       Parameter.k++;
                       flag=true;
                       break;
                   }
               }
           }
           if(flag==true)
           {
              //cout<<"[3,"<<k<<"]";
              Birth(m,Parameter); 
           }
        }
    }
    for(int i=0;i<Parameter.N;i++)
    {
        if(Parameter.Likelihood[i]!=0&&Parameter.ResultLikelihood[i]==0)
            L+=log10(Parameter.Likelihood[i])+300;
        if(Parameter.Likelihood[i]==0&&Parameter.ResultLikelihood[i]==0)
            continue;
        if(Parameter.Likelihood[i]!=0&&Parameter.ResultLikelihood[i]!=0)
            L+=log10(Parameter.Likelihood[i])-log10(Parameter.ResultLikelihood[i]);
        if(Parameter.Likelihood[i]==0&&Parameter.ResultLikelihood[i]!=0)
            L-=300;
        /*
        if(Likelihood[i]==0)
            if(ResultLikelihood[i]!=0)
                L*=(10E-308/ResultLikelihood[i]);
                else
                    L*=1.0;
            else
                if(ResultLikelihood[i]!=0)
                    L*=(Likelihood[i]/ResultLikelihood[i]);
                else
                    L*=(Likelihood[i]/10E-308);
        */
    }
    if(L>=307)
        temp=10.0E307;
    if(L<=-308)
        temp=0;
    if(L<307&&L>-308)
        temp=pow((double)10,L);
    L=temp;
        //if((L>1)&&(temp<MINProb+1))
    if(L>pow(Parameter.Lamdab,Parameter.Resultk-Parameter.k)*factorial(Parameter.Resultk,Parameter.k))
    {
        if(iter==20)
        {
            fstream TempRatioMoniter;
            TempRatioMoniter.open("LikelihoodRatio.txt",ios_base::out|ios_base::app);
            if(TempRatioMoniter.fail())
                cout<<"Cannot write likelihoodratio info!"<<endl;
            for(int i=0;i<Parameter.k;i++)
            {
                TempRatioMoniter<<Parameter.LikelihoodRatio[i]<<" ";
                //cout<<LikelihoodRatio[i]<<" ";
            }
            TempRatioMoniter<<endl;
            //cout<<endl;
            TempRatioMoniter.close();
        }                                            
        Parameter.Resultk=Parameter.k;
        for(int i=0;i<Parameter.k;i++)
        {
            Parameter.ResultPi[i]=Parameter.Pi[i];
            Parameter.ResultCenter[i]=Parameter.Center[i];
            Parameter.ResultSigma[i]=Parameter.Sigma[i];
        }
        for(int i=0;i<Parameter.N;i++)
            Parameter.ResultLikelihood[i]=Parameter.Likelihood[i];
    }
    if(Parameter.X[0]==0x00)
        cout<<"63\n"; 
    //cout<<"P";
    return 1;       
}

int UpdateLikelihoodRatio(Parameters &Parameter)
{
    time_t start;
    start=time(NULL);
    memset(Parameter.Likelihood, 0, sizeof(double)*Parameter.N);
    for(int i=0;i<Parameter.N;i++)
            for(int j=0;j<Parameter.k;j++)
            {
                 if(Parameter.SigmaChange[j]==false)
                     Parameter.LikelihoodTemp[i][j]=Parameter.Pi[j]*gaussian(Parameter.X[i][Parameter.Center[j]],Parameter.Sigma[j],Parameter);
                 else
                     Parameter.LikelihoodTemp[i][j]*=(Parameter.Pi[j]/Parameter.PiPrev[j]);                  
                 Parameter.Likelihood[i]+=Parameter.LikelihoodTemp[i][j]; 
            }
/*
    for(int i=0;i<N;i++)
            if(Likelihood[i]==0)
                Likelihood[i]=4e-300;
*/
    if(Parameter.k==1)
        Parameter.LikelihoodRatio[0]=0;
    else
    {
        double* TempLikelihood=new double[Parameter.N];
        double* LogLikelihoodRatio=new double[Parameter.k+1];
        memset(LogLikelihoodRatio,0,sizeof(double)*(Parameter.k+1));     
        for(int j=0;j<Parameter.k;j++)
        {
            Parameter.LikelihoodRatio[j]=1.0;
            for(int i=0;i<Parameter.N;i++)
            {
                    if(Parameter.Likelihood[i]==0)
                         TempLikelihood[i]=1;
                    else
                         TempLikelihood[i]=(Parameter.Likelihood[i]-Parameter.LikelihoodTemp[i][j])/(1-Parameter.Pi[j])/Parameter.Likelihood[i];
                    if(TempLikelihood[i]==0) //if a cluster is far away from all the other clusters
                    {
                        LogLikelihoodRatio[j]=-308;
                        break;
                    }
                    LogLikelihoodRatio[j]+=Parameter.Abundance[i]*log10(TempLikelihood[i]);
                    //LikelihoodRatio[j]*=pow(TempLikelihood[i],Abundance[i]);
                    //if(LikelihoodRatio[j]<0)
                        //cout<<pow(TempLikelihood[i],Abundance[i])<<endl;
            }
            if(LogLikelihoodRatio[j]<=-308)
                Parameter.LikelihoodRatio[j]=0;
            if(LogLikelihoodRatio[j]>=307)
                Parameter.LikelihoodRatio[j]=10.0e300;
            if((LogLikelihoodRatio[j]<307)&&(LogLikelihoodRatio[j]>-308))
                Parameter.LikelihoodRatio[j]=pow((double)10,LogLikelihoodRatio[j]);
        }
        //cout<<"T1";
        DELETE(TempLikelihood);
        //cout<<"T2";
        DELETE(LogLikelihoodRatio);
    }
     time_t stop;
     stop=time(NULL);
     Parameter.TimeLikelihood+=stop-start;
    //cout<<"L";
    return 1;                       
}

int ResultProcess(Parameters &Parameter)
{
    
     if(Parameter.Resultk==1)
     {
         Parameter.Pi[0]=Parameter.ResultPi[0];
         Parameter.Center[0]=Parameter.ResultCenter[0];
         Parameter.Sigma[0]=Parameter.ResultSigma[0];
         Parameter.IterNN[0]=0;
         Parameter.NN[0]=0;
         Parameter.k=Parameter.Resultk;
         for(int i=0;i<Parameter.N;i++)
         {
             Parameter.Z[i]=0;
             Parameter.NN[0]+=Parameter.Abundance[i];
             Parameter.IterNN[0]++;
         }
     }
     
     
     else
     {
     memcpy(Parameter.Pi, Parameter.ResultPi, sizeof(double)*Parameter.Resultk);
     memcpy(Parameter.Center, Parameter.ResultCenter, sizeof(int)*Parameter.Resultk);
     memcpy(Parameter.Sigma, Parameter.ResultSigma, sizeof(double)*Parameter.Resultk);
     memset(Parameter.NN,0,sizeof(int)*Parameter.Resultk);
     memset(Parameter.IterNN,0,sizeof(int)*Parameter.Resultk);
     Parameter.k=Parameter.Resultk;
     double z,temp,maxtemp;
     int Zmax=-1;
     bool flag;
     for(int i=0;i<Parameter.N;i++)
     {
          maxtemp=0;
          
          flag=false;
          for(int j=0;j<Parameter.k;j++)
          {
                  if(i==Parameter.Center[j])
                  {
                       Parameter.Z[i]=j;
                       Parameter.NN[j]+=Parameter.Abundance[i];
                       Parameter.IterNN[j]++;
                       flag=true;
                       break;
                  }                       
          }
          if(flag==false)
          {   
             for(int j=0;j<Parameter.k;j++)
             {                     
                     temp=Parameter.Pi[j]*gaussian(Parameter.X[i][Parameter.Center[j]],Parameter.Sigma[j],Parameter);
                     if(temp>maxtemp)
                     {
                         maxtemp=temp;
                         Zmax=j;
                     }
             }
             if(Zmax!=-1)
             {
                 Parameter.Z[i]=Zmax;
                 Parameter.NN[Zmax]+=Parameter.Abundance[i];
                 Parameter.IterNN[Zmax]++;
             }
             else
             {
                 Parameter.k++;
                 Parameter.NN[Parameter.k-1]=Parameter.Abundance[i];
                 Parameter.IterNN[Parameter.k-1]=1;
                 Parameter.Z[i]=Parameter.k-1;
                 Parameter.Center[Parameter.k-1]=i;
                 Parameter.Sigma[Parameter.k-1]=Parameter.Lower;
                 UpdatePi(Parameter);
             }
          }
    }
    }//else
    //cout<<"R";        
    return 1;   
}

/*
int CriteriaCalc(int t)
{
    int temp,tmp;
    //Purity
    Purity=(double)t;
    int Candidate[MAXSize];
    memset(Candidate, -1, sizeof(Candidate));
    for(int i=0;i<t;i++)
        for(int j=0;j<N;j++)
        {
                if(Z[j]==i)
                {
                    if((Candidate[i]!=-1)&&(StandardZ[j]!=Candidate[i]))
                    {
                            Purity--;
                            break;
                    }
                    else
                         Candidate[i]=StandardZ[j];
                }
        }
    tmp=t;
    for(int i=0; i<t; i++)
        if(NN[i]==0)
        {
            Purity--;
            tmp--;
        }
    Purity/=(double)tmp;
    
    //VI
    double HResult=0.0;
    double HStandard=0.0;
    double HRandS=0.0;
    for(int i=0;i<t;i++)
    {
        if(NN[i]!=0)
            HResult-=(double)NN[i]/N*log((double)NN[i]/N);
    }
    for(int i=0;i<Standardk;i++)
        HStandard-=(double)StandardNN[i]/N*log((double)StandardNN[i]/N);
    for(int i=0;i<t;i++)
        for(int j=0;j<Standardk;j++)
        {
            temp=0;
            for(int l=0;l<N;l++)
            {
                    if((Z[l]==i)&&(StandardZ[l]==j))
                        temp++;
            }
            if(temp>0)
                HRandS+=(double)temp/N*log((double)temp/NN[i]/StandardNN[j]*N);
        }
    VI=HResult+HStandard-2*HRandS;
    //Edit Distance
    EditDistance=0;
    for(int i=0;i<t;i++)
        for(int j=0;j<Standardk;j++)
        {
            temp=0;
            for(int l=0;l<N;l++)
            {
                    if((Z[l]==i)&&(StandardZ[l]==j))
                    {
                        temp=1;
                        EditDistance++;
                        break;
                    }
            }
        }
    EditDistance*=2;
    EditDistance-=(tmp+Standardk);
    
    return 1;        
}
*/

double factorial(int t, int l)
{
    if(t<l)
        return l*factorial(t,l-1);
    if(t>l)
        return (double)1/t*factorial(t-1,l);
    if(t==l)
        return 1;
}

double gaussian(double x, double y, Parameters &Parameter)
{
    if(x>MAX(15,6*Parameter.Upper))
        return 0;
    else
        return gsl_ran_gaussian_pdf(x,y);
}

char** readFastASeq(char *fname, unsigned int *num_seq, Parameters &Parameter)
{
	char dummy[MAX_LINE], **seq = NULL, tempdummy[MAX_LINE];
	unsigned int i = 0, len = 0, lentotal = 0;
	FILE *fp = fopen(fname, "r");

	if(fp == NULL) printf("Cannot open `%s' file.\n", fname);
	
    *num_seq = 0;
	while(feof(fp) == 0){
		fgets(dummy, MAX_LINE, fp);
		if(dummy[0] == '>')
			(*num_seq)++;
	}
	printf("num_seq = %u \n", *num_seq);
	if(*num_seq==0)
	{
        printf("file emply!");
        return NULL;
    }
    
	rewind(fp);
    
	MALLOC(Parameter.seq, char*, *num_seq);
	MALLOC(Parameter.SeqID, char*, *num_seq);

	fgets(dummy, MAX_LINE, fp); /* Read first line */
	for(i = 0; i < *num_seq; i++){
		if(dummy[0] == '>'){
            len = strlen(dummy);
            len += (dummy[len-1] == '\n')? -1:0;
            MALLOC(Parameter.SeqID[i], char, len+1);
            memcpy(Parameter.SeqID[i], dummy, sizeof(char)*len);        
			if(fgets(dummy, MAX_LINE, fp) == NULL) printf("This is not a valid fasta file (no read for a header).\n");
			lentotal=0;   
			len = strlen(dummy);
			if(len == MAX_LINE) printf("There is very long line (> %u).\n", MAX_LINE);
			len += (dummy[len-1] == '\n')? -1:0;
            memcpy(tempdummy+lentotal, dummy, sizeof(char)*len);			
			lentotal+=len;
			if(fgets(dummy, MAX_LINE, fp) == NULL){ //end of file
			    MALLOC(Parameter.seq[i], char, lentotal+1);
			    memcpy(Parameter.seq[i], tempdummy, sizeof(char)*(lentotal));
			    fclose(fp);
			    return Parameter.seq;
            }
			while(dummy[0]!='>'){
            if(dummy[0] == '\n' || dummy[0] == '\0'){
			    MALLOC(Parameter.seq[i], char, lentotal+1);
			    memcpy(Parameter.seq[i], tempdummy, sizeof(char)*(lentotal));
			    fclose(fp);
			    return Parameter.seq;
            }                
			len = strlen(dummy);
			if(len == MAX_LINE) printf("There is very long line (> %u).\n", MAX_LINE);
			len += (dummy[len-1] == '\n')? -1:0;
			memcpy(tempdummy+lentotal, dummy, sizeof(char)*len);
			lentotal+=len;
            if(fgets(dummy, MAX_LINE, fp) == NULL){ /* read next fasta element */
			    MALLOC(Parameter.seq[i], char, lentotal+1);
			    memcpy(Parameter.seq[i], tempdummy, sizeof(char)*(lentotal));
			    fclose(fp);
			    return Parameter.seq;
            }
            }
            MALLOC(Parameter.seq[i], char, lentotal+1);
            memcpy(Parameter.seq[i], tempdummy, sizeof(char)*(lentotal));
		}
		else if(dummy[0] == '\n' || dummy[0] == '\0') break;
		else printf("This is not a valid fasta file.\n");
	}
	if(fclose(fp)!=0){
         printf("Fail to close fasta file.\n");
    }
	return Parameter.seq;
}

void NewMem(int nn,Parameters &Parameter)
{
     //Global Variables
     //*******************************
     Parameter.Distrk=new double[Parameter.Step];
     Parameter.Time=new double[Parameter.Step];
     Parameter.Pi=new double[nn];
     Parameter.PiPrev=new double[nn];
     Parameter.Center=new int[nn];
     Parameter.Sigma=new double[nn];
     Parameter.SigmaChange=new bool[nn];
     Parameter.Likelihood=new double[nn];
     Parameter.LikelihoodRatio=new double[nn]; //for UpdateLikelihood;
     //************************************
     Parameter.LikelihoodTemp=new double*[nn];
     //************************************
     if(nn<LIMIT) //2010-1-19 if sequence number compatible, generate distance matrix
     {
         Parameter.X=new float*[nn];         //Distance Matrix
         for(int i=0;i<nn;i++)
             Parameter.X[i]=new float[nn];
         for(int i=0;i<nn;i++)
             Parameter.LikelihoodTemp[i]=new double[nn];
     }
     Parameter.Z=new int[nn];       //missing data
     Parameter.NN=new int[nn];       //NN[i-1]=# of samples in i-th cluster
     Parameter.IterNN=new int[nn];
     Parameter.Abundance=new int[nn];
     Parameter.Centerlikelihood=new double[nn];  //square sum of current center
     Parameter.Weight=new bool[nn];  //for Birth-Death control;
     //*******************************

     //Result Process
     //*******************************
     Parameter.ResultPi=new double[nn];
     Parameter.ResultCenter=new int [nn];
     Parameter.ResultSigma=new double[nn];
     Parameter.ResultLikelihood=new double[nn];
    //*******************************

    //Hyperparameters
    //*******************************
    //For Pi
    Parameter.Gamma=new double[nn];   
    return;
}

void DeleteMem(int nn, int iter,Parameters &Parameter)
{
     //Global Variables
     //*******************************
     DELETE(Parameter.Distrk);
     DELETE(Parameter.Time);
     DELETE(Parameter.Pi);
     DELETE(Parameter.PiPrev);
     DELETE(Parameter.Center);
     DELETE(Parameter.Sigma);
     DELETE(Parameter.SigmaChange);
     DELETE(Parameter.Likelihood);
     DELETE(Parameter.LikelihoodRatio);
     for(int i=0;i<nn;i++)
     {
             FREE(Parameter.SeqID[i]);
             FREE(Parameter.seq[i]);
     }
     FREE(Parameter.SeqID);
     FREE(Parameter.seq);
     if(iter!=0)
     {
     for(int i=0;i<nn;i++)
     {    
        for(int j=0;j<MIN(20,Parameter.Abundance[i]);j++)
        {
            DELETE(Parameter.EigenSeqs[i][j]); //delete[]
        }
        DELETE(Parameter.EigenSeqs[i]); //delete[]
     }
     DELETE(Parameter.EigenSeqs); //delete[]
     }  
     //************************************
     if(nn<LIMIT)   //2010-1-19 Distance matrix not generated if sequence number is too large
     {
         for(int i=0;i<nn;i++)
             DELETE(Parameter.X[i]);  //delete[]
         DELETE(Parameter.X);         //Distance Matrix
         //for UpdateLikelihood;
         //************************************
         for(int i=0;i<nn;i++) 
             DELETE(Parameter.LikelihoodTemp[i]); //delete[]
         DELETE(Parameter.LikelihoodTemp);  //delete[]
     }
     DELETE(Parameter.Z);       //missing data
     DELETE(Parameter.NN);       //NN[i-1]=# of samples in i-th cluster
     DELETE(Parameter.IterNN);
     DELETE(Parameter.Abundance);
     DELETE(Parameter.Centerlikelihood);  //square sum of current center
     DELETE(Parameter.Weight);  
     //for Birth-Death control;
     //*******************************

     //Result Process
     //*******************************
     DELETE(Parameter.ResultPi);
     DELETE(Parameter.ResultCenter);
     DELETE(Parameter.ResultSigma);
     DELETE(Parameter.ResultLikelihood);
    //*******************************

    //Hyperparameters
    //*******************************
    //For Pi
    DELETE(Parameter.Gamma);
    return;
}

BayesianClusteringResult BayesianClustering(const char* fname, int step, double lower, double upper, BayesianClusteringResult &CResult, int iter, Parameters &Parameter)
{ 
    Initialize(Parameter);   
    int len;
   	time_t start,stop;
    Parameter.Step=step;
    Parameter.Lower=lower;
    Parameter.Upper=upper;
    InitRandomGenerator(Parameter);
    LoadDistanceMatrix(fname,iter,Parameter);        
    if((Parameter.N>=LIMIT&&iter==20)||(Parameter.N==0)) //2010-1-19  If sequence number not compatible for MergingUsingBC
    {
         if(Parameter.N==0)
         {
             CResult.ClusterNumber=0;
             return CResult;
         }
         Parameter.k=Parameter.N;
         for(int i=0;i<Parameter.N;i++)
         {
             Parameter.NN[i]=Parameter.Abundance[i];
             Parameter.Center[i]=i;
             Parameter.Sigma[i]=Parameter.Upper;
             Parameter.IterNN[i]=1;
             Parameter.Z[i]=i;
         }
    }    
    else    //2010-1-19 If sequence number compatible, continue MergingUsingBC
    {       
    GenerateInitial(iter,Parameter);    
    //*******for time counting***********
	start=time(NULL);
    memset(Parameter.Time,0,sizeof(double)*Parameter.Step);
    memset(Parameter.Distrk,0,sizeof(double)*Parameter.Step);
    //***********************************
    if(iter==20)
    {
        cout<<Parameter.Step<<endl;                
        cout<<"******Progress******"<<endl;
        cout<<"*";
    }
    for(int i=0;i<Parameter.Step;i++)
    {
            if(iter==20)
            {
                if((20*i/Parameter.Step)!=(20*(i-1)/Parameter.Step))
                    cout<<"*";
            } 
            ProcessSampling(i,iter,Parameter);
            /*             
            if(i>=2000)
            {
                //for result record
                CriteriaCalc(k);
                Pu=(double)(Pu*i+Purity)/(i+1);
                V=(double)(V*i+VI)/(i+1);
                ED=(double)(ED*i+EditDistance)/(i+1);           
            }
            */
    }
    if(iter==20)
        cout<<endl;
    //for time counting
	stop=time(NULL);
	//cout<<"Time Cost: "<<stop-start<<" sec"<<endl;
    //fs<<Pu<<","<<V<<","<<ED<<endl;
    //cout<<Parameter.Resultk<<endl;
    //cout<<Parameter.k<<endl;
    ResultProcess(Parameter); 
    }
    CResult.ClusterNumber=Parameter.k;
    int TempK=Parameter.k;
    bool TempKFlag=false;
    if(Parameter.k==1)
    {
        TempK++;
        TempKFlag=true;
    }  //for else 2011-1-19
    CResult.ClusterCenterID=new char*[TempK];
    CResult.ClusterCenterSeq=new char*[TempK];
    CResult.ClusterStandardDeviation=new double[TempK];
    CResult.ClusterSize=new int[TempK];
    CResult.ClusterIterationSize=new int[TempK];
    CResult.ClusterMember=new string*[TempK];
    CResult.ClusterMemberSeq=new string*[TempK];
    CResult.ClusterMemberFlag=new bool*[TempK];
    Parameter.TempNN=new int[TempK];
    if(TempKFlag==true)    //When Resultk=1
    {
        CResult.ClusterCenterID[TempK-1]=NULL;
        CResult.ClusterCenterSeq[TempK-1]=NULL;
        CResult.ClusterStandardDeviation[TempK-1]=0;
        CResult.ClusterSize[TempK-1]=0;
        CResult.ClusterIterationSize[TempK-1]=0;
        CResult.ClusterMember[TempK-1]=NULL;
        CResult.ClusterMemberSeq[TempK-1]=NULL;
        CResult.ClusterMemberFlag[TempK-1]=NULL;
        Parameter.TempNN[TempK-1]=0;
    }
    for(int i=0;i<Parameter.k;i++)
    {
        len=strlen(Parameter.SeqID[Parameter.Center[i]]);
        CResult.ClusterCenterID[i]=new char[len+1];
        memset(CResult.ClusterCenterID[i], 0x00, sizeof(char)*(len+1));
        memcpy(CResult.ClusterCenterID[i], Parameter.SeqID[Parameter.Center[i]], sizeof(char)*len);
        len=strlen(Parameter.seq[Parameter.Center[i]]);
        CResult.ClusterCenterSeq[i]=new char[len+1];
        memset(CResult.ClusterCenterSeq[i],0x00, sizeof(char)*(len+1));
        memcpy(CResult.ClusterCenterSeq[i], Parameter.seq[Parameter.Center[i]], sizeof(char)*len);
        CResult.ClusterStandardDeviation[i]=Parameter.Sigma[i];
        CResult.ClusterSize[i]=Parameter.NN[i];
        CResult.ClusterIterationSize[i]=Parameter.IterNN[i];
        CResult.ClusterMember[i]=new string[Parameter.IterNN[i]];
        CResult.ClusterMemberFlag[i]=new bool[Parameter.IterNN[i]];
        memset(CResult.ClusterMemberFlag[i],0,sizeof(bool)*Parameter.IterNN[i]);
        CResult.ClusterMemberSeq[i]=new string[MIN(20,Parameter.NN[i])];
    }
    memset(Parameter.TempNN, 0,  sizeof(int)*Parameter.k);   
    for(int i=0;i<Parameter.N;i++)
    {
        CResult.ClusterMember[Parameter.Z[i]][Parameter.TempNN[Parameter.Z[i]]].assign(Parameter.SeqID[i]);
        Parameter.TempNN[Parameter.Z[i]]++;
    }
    for(int i=0;i<Parameter.k;i++)
    {
        if(Parameter.TempNN[i]!=Parameter.IterNN[i])
            cout<<"Unexpected Result Error in Block"<<endl;
    }
    //cout<<"Basic Operations Done"<<endl;
    //*****************************************Generate Eigen Seqs*********************************************************
    memset(Parameter.TempNN,0,sizeof(int)*Parameter.k);
    int TempSeqNum, temp;
    for(int i=0;i<Parameter.N;i++)
    {
        Parameter.TempNN[Parameter.Z[i]]+=MIN(20,Parameter.Abundance[i]);
    }
    if(iter!=0)
    {
    for(int i=0;i<Parameter.k;i++)
    {

        unsigned long int Range=MIN(20,Parameter.TempNN[i]);
        if(Range==20)
        {
            for(int j=0;j<20;j++)
            {
                TempSeqNum=gsl_rng_uniform_int(Parameter.r, Parameter.TempNN[i])+1;
                temp=0;
                for(int l=0;l<Parameter.N;l++)
                {
                    if(Parameter.Z[l]==i)
                    {
                        temp+=MIN(20,Parameter.Abundance[l]);
                        if(temp>=TempSeqNum)
                        {
                            CResult.ClusterMemberSeq[i][j].assign(Parameter.EigenSeqs[l][TempSeqNum-1-(temp-MIN(20,Parameter.Abundance[l]))]);
                            break;
                        }
                    }
                }
            }
        }
        else
        {
            int TempCount=0;
            for(int j=0;j<Parameter.N;j++)
            {
                if(Parameter.Z[j]==i)
                {
                    for(int l=0;l<MIN(20,Parameter.Abundance[j]);l++)
                    {
                        CResult.ClusterMemberSeq[i][TempCount].assign(Parameter.EigenSeqs[j][l]);
                        TempCount++;
                    }
                }
            }
        }
    }
    }//if iter!=0
    
    else
    {
        for(int i=0;i<Parameter.k;i++)
        {
            for(int j=0;j<MIN(20,Parameter.TempNN[i]);j++)
            {
               TempSeqNum=gsl_rng_uniform_int(Parameter.r, Parameter.TempNN[i]);
                temp=0;
                for(int l=0;l<Parameter.N;l++)
                {
                    if(Parameter.Z[l]==i)
                    {
                        temp+=MIN(20,Parameter.Abundance[l]);
                        if(temp>=TempSeqNum)
                        {
                            CResult.ClusterMemberSeq[i][j].assign(Parameter.seq[l]);
                            break;
                        }//if temp>TempSeqNum
                    }//if Z[l]==i
                }//l
            }//j
        }//i
    }//if iter==0
    //cout<<"Choosing Seqs Compelete"<<endl;
    //********************************************************************************************************************
            
    DELETE(Parameter.TempNN);
    //cout<<TimeZ<<" "<<TimeCenter<<" "<<TimePi<<" "<<TimeSigma<<" "<<TimeBirth<<" "<<TimeDeath<<" "<<TimeLikelihood<<endl;
    //cout<<"Build Complete"<<endl;
    gsl_rng_free(Parameter.r);
    DeleteMem(Parameter.N,iter,Parameter);
    return CResult;
}
