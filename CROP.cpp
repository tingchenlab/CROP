#include "CROP.h"


using namespace std;

void InitCROPParameters(CROPParameters &CParameter)
{
    memset(CParameter.Block,1,sizeof(int)*20);
    CParameter.SeqNum=0;
    CParameter.UniqueSeqNum=0;
    CParameter.PrevIterk=0;
    CParameter.FinalCenterID=NULL;
    CParameter.FinalCenterSeq=NULL;
    CParameter.FinalSD=NULL;
    CParameter.FinalSize=NULL;
    CParameter.r=NULL;
}

bool ChompFile(const char* fname, int block, CROPParameters &CParameter)
{
    fstream in,out;
    string ResultName;
    string ResultNamePrefix(fname);
    ResultNamePrefix+=".BlockSeq";
    string ResultNameNum="1";
    string Line;
    char Buffer[ReadSize];
    char temp[ReadSize];
    in.open(fname);
    if(in.fail())
    {
        cout<<"Read Fasta File Error"<<endl;
        return 0;
    }
    while(!in.eof())
    {
        in.getline(temp, ReadSize);
        if(temp[0]=='>')
            CParameter.UniqueSeqNum++;
    }
    if(CParameter.UniqueSeqNum==0)
    {
        cout<<"Not A Valid Fasta File!"<<endl;
        return 0;
    }    
    in.clear();   
    in.seekg(0,ios::beg);
    cout<<CParameter.UniqueSeqNum<<endl;
    ResultName=ResultNamePrefix+ResultNameNum; 
    out.open(ResultName.c_str(),ios_base::out);
    if(out.fail())
    {
        cout<<"Write Fasta File Error"<<endl;
    }
    in.getline(temp,ReadSize);
    Line.assign(temp);
    if(temp[0]!='>')
    {
        cout<<"Bad Header"<<endl;
        return 0;
    }
    out<<Line<<endl;
    in.getline(temp,ReadSize);
    Line.assign(temp);
    int numblock=1;
    for(int i=0; i<CParameter.UniqueSeqNum; i++)
    {
        while(temp[0]!='>')
        {
            out<<Line<<endl;
            in.getline(temp,ReadSize);
            if(temp[0] == '\n' || temp[0] == '\0') break;
            Line.assign(temp);
        }
        if(floor((double)MAX(i-1,0)/CParameter.UniqueSeqNum*block)!=floor((double)i/CParameter.UniqueSeqNum*block))
        { 
            numblock++;
            out.close();
            memset(Buffer, 0, sizeof(char)*ReadSize);
            itoa(numblock, Buffer, 10);
            ResultNameNum.assign(Buffer);
            ResultName=ResultNamePrefix+ResultNameNum;
            out.open(ResultName.c_str(),ios_base::out);
            if(out.fail())
            {
                cout<<"Write Fasta File Error"<<endl;
            }
        }
        if(i<CParameter.UniqueSeqNum-1)
        {
            out<<Line<<endl;
            in.getline(temp,ReadSize);
            Line.assign(temp);
        }
        else
            out.close();
    }
    in.close();
    return 1;
}

bool ChompInterFile(const char* fname, int block, CROPParameters &CParameter)
{
    int tempstat=0;
    fstream in,out;
    string ResultName;
    string ResultNamePrefix(fname);
    ResultNamePrefix+=".BlockSeq";
    string ResultNameNum="1";
    string Line;
    int tempint;
    double tempdouble;
    char Buffer[ReadSize];
    char temp[ReadSize];
    string TempCenters(fname);
    TempCenters+=".TempCenters";
    in.open(TempCenters.c_str());
    if(in.fail())
    {
        cout<<"Read Inter Result Fasta File Error"<<endl;
        return 0;
    }
    while(!in.eof())
    {
        in.getline(temp, ReadSize);
        if(temp[0]=='>')
            tempstat++;
    }
    if(tempstat==0)
    {
        cout<<"Not A Valid Fasta File!"<<endl;
        return 0;
    }    
    in.clear();   
    in.seekg(0,ios::beg);
    cout<<"Number of Seqs in File Before Chomp = "<<tempstat<<endl;
    ResultName=ResultNamePrefix+ResultNameNum; 
    out.open(ResultName.c_str(),ios_base::out);
    int numblock=1;
    if(out.fail())
    {
        cout<<"Write Fasta File Error"<<endl;
    }
    for(int i=0;i<CParameter.Finalk;i++)
    {
        if(floor((double)MAX(i-1,0)/CParameter.Finalk*block)!=floor((double)i/CParameter.Finalk*block))
        {
            numblock++;
            out.close();
            itoa(numblock, Buffer, 10);
            ResultNameNum.assign(Buffer);
            ResultName=ResultNamePrefix+ResultNameNum;
            out.open(ResultName.c_str(),ios_base::out);
            if(out.fail())
            {
                cout<<"Write Fasta File Error"<<endl;
            }
        }
        //*******Center Seq ID**********
        in.getline(temp,ReadSize);
        Line.assign(temp);
        out<<Line<<endl;
        //*******Center Seq************
        in.getline(temp,ReadSize);
        Line.assign(temp);
        out<<Line<<endl;
        //********Abudance and Sigma********
        in>>tempint;
        out<<tempint<<endl;
        in>>tempdouble;
        out<<tempdouble<<endl;
        in.getline(temp,ReadSize);
        //***********Eigen Seqs**************
        for(int j=0;j<MIN(tempint,20);j++)
        {
            in.getline(temp,ReadSize);
            Line.assign(temp);
            out<<Line<<endl;
        }
        if(i==CParameter.Finalk-1)
            out.close();
    }
    in.close();
    return 1;
}

int LoadingList(const char* fname, CROPParameters &CParameter)
{
    int N;
    fstream SeqList;
    string SeqListName(fname);
    char dummy[MAX_LINE];
    SeqListName+=".list";
    cout<<"Loading "<<SeqListName<<endl;
    SeqList.open(SeqListName.c_str());
    if(SeqList.fail())
    {
        cout<<"Cannot open name file"<<endl;
        return 0;
    }
    SeqList>>N;
    SeqList.getline(dummy, MAX_LINE);
    for(int i=0;i<N;i++)
    {
        string temp1, temp2, temp;
        int len = 0;
        SeqList>>temp1;
        SeqList.getline(dummy, MAX_LINE);//read .Name
        //len = strlen(dummy);
        //len += (dummy[len-1] == '\n')? -1:0;
        //temp1.assign(dummy,len);
        SeqList>>temp2;
        SeqList.getline(dummy, MAX_LINE);//read .Members
        len = strlen(dummy);        
        CParameter.List[temp1]=temp2;  //build map
    }
    SeqList.close();
    return N;
}

void LoadBlockFile(string fname, TempBlockFile& Temp)
{
    fstream Rare;
    string Line;
    Rare.open(fname.c_str());
    Temp.SeqNum=0;
    if(Rare.fail()){cout<<"cannot open rare clusters file\n"; return;}
    while(!Rare.eof())
    {
        getline(Rare, Line);
        if(Line[0]=='>')
            Temp.SeqNum++;
    }
    if(Temp.SeqNum==0)
        return;
    Temp.SeqID=new string[Temp.SeqNum];
    Temp.Seq=new string[Temp.SeqNum];
    Temp.EigenSeqs=new string*[Temp.SeqNum];
    Temp.ClusterSize=new int[Temp.SeqNum];
    Temp.ClusterStandardDeviation=new double[Temp.SeqNum];
    Temp.Assignment=new int[Temp.SeqNum];
    memset(Temp.Assignment, -1, sizeof(int)*(Temp.SeqNum));
    Rare.clear();   
    Rare.seekg(0,ios::beg);
    for(int i=0;i<Temp.SeqNum;i++)
    {
        Rare>>Temp.SeqID[i];
        Temp.SeqID[i]=Temp.SeqID[i].substr(1);
        getline(Rare,Line);
        Rare>>Temp.Seq[i];
        getline(Rare,Line);
        Rare>>Temp.ClusterSize[i];
        getline(Rare,Line);  //for cluster size
        Rare>>Temp.ClusterStandardDeviation[i];
        getline(Rare,Line);  //for cluster standard deviation
        Temp.EigenSeqs[i]=new string[MIN(20,Temp.ClusterSize[i])];
        for(int j=0;j<MIN(20,Temp.ClusterSize[i]);j++)
        {
            Rare>>Temp.EigenSeqs[i][j];
            getline(Rare, Line);
        }
    } 
    Rare.close();
    cout<<"load rare clusters file complete\n";
    return;
}

void CopyClusteringResult(const BayesianClusteringResult& src, BayesianClusteringResult& dest, bool ctrl)
{
        dest.ClusterNumber=src.ClusterNumber;
        dest.ClusterCenterID=new char*[dest.ClusterNumber];
        dest.ClusterCenterSeq=new char*[dest.ClusterNumber];
        dest.ClusterStandardDeviation=new double[dest.ClusterNumber];
        dest.ClusterSize=new int[dest.ClusterNumber];
        dest.ClusterIterationSize=new int[dest.ClusterNumber];
        if(ctrl==true)
        {
            dest.ClusterMember=new string*[dest.ClusterNumber];
            dest.ClusterMemberFlag=new bool*[dest.ClusterNumber];
            dest.ClusterMemberSeq=new string*[dest.ClusterNumber];
        }
        else
        {
            dest.ClusterMember=NULL;
            dest.ClusterMemberFlag=NULL;
            dest.ClusterMemberSeq=NULL;
        }
        for(int i=0;i<dest.ClusterNumber;i++)
        {
            //center ID    
            int len=strlen(src.ClusterCenterID[i]);
            dest.ClusterCenterID[i]=new char[len+1];
            memset(dest.ClusterCenterID[i], 0x00, sizeof(char)*(len+1));
            memcpy(dest.ClusterCenterID[i], src.ClusterCenterID[i], sizeof(char)*len);
            //center seq
            len=strlen(src.ClusterCenterSeq[i]);
            dest.ClusterCenterSeq[i]=new char[len+1];
            memset(dest.ClusterCenterSeq[i], 0x00, sizeof(char)*(len+1));
            memcpy(dest.ClusterCenterSeq[i], src.ClusterCenterSeq[i], sizeof(char)*len);
            //sd
            dest.ClusterStandardDeviation[i]=src.ClusterStandardDeviation[i];
            //size
            dest.ClusterSize[i]=src.ClusterSize[i];
            dest.ClusterIterationSize[i]=src.ClusterIterationSize[i];
            //member
            if(ctrl==true)
            {
                dest.ClusterMember[i]=new string[dest.ClusterIterationSize[i]];
                dest.ClusterMemberFlag[i]=new bool[dest.ClusterIterationSize[i]]; 
                dest.ClusterMemberSeq[i]=new string[MIN(20,dest.ClusterSize[i])];       
                for(int j=0;j<dest.ClusterIterationSize[i];j++)
                {
                    dest.ClusterMember[i][j]=src.ClusterMember[i][j];
                    dest.ClusterMemberFlag[i][j]=src.ClusterMemberFlag[i][j];
                }
                for(int j=0;j<MIN(20,dest.ClusterSize[i]);j++)
                    dest.ClusterMemberSeq[i][j]=src.ClusterMemberSeq[i][j];
            }
        }
        cout<<"copied "<<dest.ClusterNumber<<" temp clustering results\n";
        return;
}

void MergingClusteringResult(BayesianClusteringResult& src1, BayesianClusteringResult& src2, BayesianClusteringResult& dest)
{
     cout<<"merging clustering results\n";
     dest.ClusterNumber=src1.ClusterNumber+src2.ClusterNumber;
     dest.ClusterCenterID=new char*[dest.ClusterNumber];
     dest.ClusterCenterSeq=new char*[dest.ClusterNumber];
     dest.ClusterStandardDeviation=new double[dest.ClusterNumber];
     dest.ClusterSize=new int[dest.ClusterNumber];
     dest.ClusterIterationSize=new int[dest.ClusterNumber];    
     dest.ClusterMember=new string*[dest.ClusterNumber];
     dest.ClusterMemberFlag=new bool*[dest.ClusterNumber];
     dest.ClusterMemberSeq=new string*[dest.ClusterNumber];
     for(int i=0;i<src1.ClusterNumber;i++)
     {
         //center ID
         //cout<<src1.ClusterCenterID[i]<<"\n";  
         int len=strlen(src1.ClusterCenterID[i]);
         dest.ClusterCenterID[i]=new char[len+1];
         memset(dest.ClusterCenterID[i], 0x00, sizeof(char)*(len+1));
         memcpy(dest.ClusterCenterID[i], src1.ClusterCenterID[i], sizeof(char)*len);
         //cout<<"ID["<<i<<"]\n";
         //center seq
         //cout<<src1.ClusterCenterSeq[i]<<"\n";
         len=strlen(src1.ClusterCenterSeq[i]);
         dest.ClusterCenterSeq[i]=new char[len+1];
         memset(dest.ClusterCenterSeq[i], 0x00, sizeof(char)*(len+1));
         memcpy(dest.ClusterCenterSeq[i], src1.ClusterCenterSeq[i], sizeof(char)*len);
         //cout<<"Seq["<<i<<"]\n";
         //sd
         //cout<<src1.ClusterStandardDeviation[i]<<"\n";
         dest.ClusterStandardDeviation[i]=src1.ClusterStandardDeviation[i];
         //cout<<"SD["<<i<<"]\n";
         //size
         //cout<<src1.ClusterSize[i]<<"\t"<<src1.ClusterIterationSize[i]<<"\n";
         dest.ClusterSize[i]=src1.ClusterSize[i];
         dest.ClusterIterationSize[i]=src1.ClusterIterationSize[i];
         //cout<<"Size["<<i<<"]\n";           
         //member
         dest.ClusterMember[i]=new string[dest.ClusterIterationSize[i]];
         dest.ClusterMemberFlag[i]=new bool[dest.ClusterIterationSize[i]]; 
         dest.ClusterMemberSeq[i]=new string[MIN(20,dest.ClusterSize[i])]; 
         //cout<<"MemberNew["<<i<<"]\n";      
         for(int j=0;j<dest.ClusterIterationSize[i];j++)
         {
              //cout<<src1.ClusterMember[i][j]<<"\t"<<src1.ClusterMemberFlag[i][j]<<"\n";
              dest.ClusterMember[i][j]=src1.ClusterMember[i][j];
              dest.ClusterMemberFlag[i][j]=src1.ClusterMemberFlag[i][j];
              //cout<<"Member["<<i<<"]["<<j<<"]\n";
         }
         for(int j=0;j<MIN(20,dest.ClusterSize[i]);j++)
         {
              //cout<<src1.ClusterMemberSeq[i][j]<<"\n";
              dest.ClusterMemberSeq[i][j]=src1.ClusterMemberSeq[i][j];
              //cout<<"MemberSeq["<<i<<"]["<<j<<"]\n";
         }
     }
     for(int i=src1.ClusterNumber;i<dest.ClusterNumber;i++)
     {
         int ii=i-src1.ClusterNumber;
         //center ID
         //cout<<src2.ClusterCenterID[ii]<<"\n";   
         int len=strlen(src2.ClusterCenterID[ii]);
         dest.ClusterCenterID[i]=new char[len+1];
         memset(dest.ClusterCenterID[i], 0x00, sizeof(char)*(len+1));
         memcpy(dest.ClusterCenterID[i], src2.ClusterCenterID[ii], sizeof(char)*len);
         //cout<<"ID["<<i<<"]\n";
         //center seq
         //cout<<src2.ClusterCenterSeq[ii]<<"\n";
         len=strlen(src2.ClusterCenterSeq[ii]);
         dest.ClusterCenterSeq[i]=new char[len+1];
         memset(dest.ClusterCenterSeq[i], 0x00, sizeof(char)*(len+1));
         memcpy(dest.ClusterCenterSeq[i], src2.ClusterCenterSeq[ii], sizeof(char)*len);
         //cout<<"Seq["<<i<<"]\n";
         //sd
         //cout<<src2.ClusterStandardDeviation[ii]<<"\n";
         dest.ClusterStandardDeviation[i]=src2.ClusterStandardDeviation[ii];
         //cout<<"SD["<<i<<"]\n";
         //size
         //cout<<src2.ClusterSize[ii]<<"\t"<<src2.ClusterIterationSize[ii]<<"\n";
         dest.ClusterSize[i]=src2.ClusterSize[ii];
         dest.ClusterIterationSize[i]=src2.ClusterIterationSize[ii];
         //cout<<"Size["<<i<<"]\n";         
         //member
         dest.ClusterMember[i]=new string[dest.ClusterIterationSize[i]];
         dest.ClusterMemberFlag[i]=new bool[dest.ClusterIterationSize[i]]; 
         dest.ClusterMemberSeq[i]=new string[MIN(20,dest.ClusterSize[i])]; 
         //cout<<"MemberNew["<<i<<"]\n";          
         for(int j=0;j<dest.ClusterIterationSize[i];j++)
         {
              //cout<<src2.ClusterMember[ii][j]<<"\t"<<src2.ClusterMemberFlag[ii][j]<<"\n";
              dest.ClusterMember[i][j]=src2.ClusterMember[ii][j];
              dest.ClusterMemberFlag[i][j]=src2.ClusterMemberFlag[ii][j];
              //cout<<"Member["<<i<<"]["<<j<<"]\n";
         }
         for(int j=0;j<MIN(20,dest.ClusterSize[i]);j++)
         {
              //cout<<src2.ClusterMemberSeq[ii][j]<<"\n";
              dest.ClusterMemberSeq[i][j]=src2.ClusterMemberSeq[ii][j];
              //cout<<"MemberSeq["<<i<<"]["<<j<<"]\n";
         }
     }          
     cout<<"merging clustering result complete\n";
     return;
}

void ClearClusteringResult(BayesianClusteringResult& trash)
{
     for(int i=0;i<trash.ClusterNumber;i++)
     {
         delete[] trash.ClusterCenterID[i];
         delete[] trash.ClusterCenterSeq[i];
         delete[] trash.ClusterMember[i];
         delete[] trash.ClusterMemberSeq[i];
         delete[] trash.ClusterMemberFlag[i];
     }
     delete[] trash.ClusterCenterID;
     delete[] trash.ClusterCenterSeq;
     delete[] trash.ClusterStandardDeviation;
     delete[] trash.ClusterSize;
     delete[] trash.ClusterIterationSize;     
     delete[] trash.ClusterMember;
     delete[] trash.ClusterMemberSeq;
     delete[] trash.ClusterMemberFlag;
     cout<<"cleared "<<trash.ClusterNumber<<" temp clustering results\n";
     return;
}

void TempBlockFileDelete(TempBlockFile& Temp)
{
    for(int i=0;i<Temp.SeqNum;i++)
        delete[] Temp.EigenSeqs[i];
    delete[] Temp.EigenSeqs;
    delete[] Temp.ClusterSize;
    delete[] Temp.ClusterStandardDeviation;
    delete[] Temp.SeqID;
    delete[] Temp.Seq;
    delete[] Temp.Assignment;
    cout<<"cleared temp block file\n";
    return;
}


void MappingRareToAbundant(string fname, BayesianClusteringResult& AbundantResult, BayesianClusteringResult& MappedResult, double Upper)
{
    cout<<"mapping rare clusters to abundant clusters\n";
    int matrix = EDNAFULL;
	unsigned int seq_id = 0, tar_id = 0;
	alignment (*alignTool) (const char *, const char *, alignment, int, int) = NULL, align = {0, 0, 0};
	alignTool = NW_alignmentEx;
    TempBlockFile RareCluster;
    LoadBlockFile(fname,RareCluster);
    if(RareCluster.SeqNum!=0)
    {
    CopyClusteringResult(AbundantResult,MappedResult, false);   
    int** RareAssignment=NULL;
    RareAssignment=new int*[MappedResult.ClusterNumber];   
    int* count=NULL;
    count=new int[MappedResult.ClusterNumber];
    memset(count,0,sizeof(int)*MappedResult.ClusterNumber);
    for(int i=0;i<MappedResult.ClusterNumber;i++)
    {
        RareAssignment[i]=new int[RareCluster.SeqNum];
        memset(RareAssignment[i],-1,sizeof(int)*(RareCluster.SeqNum));
    }
    //#pragma omp parallel for
    for(seq_id=0;seq_id<RareCluster.SeqNum;seq_id++)
    {
        for(tar_id=0;tar_id<MappedResult.ClusterNumber;tar_id++)
        {
			align = alignTool(RareCluster.Seq[seq_id].c_str(), MappedResult.ClusterCenterSeq[tar_id], align, matrix, 1);
			float dist=(float)(1.0-(double)align.matches/(double)SHORT_ALIGN(align))*100.0;
			if(dist<Upper*2)
			{
                //#pragma omp critical
                //{
                    MappedResult.ClusterIterationSize[tar_id]++;
                    MappedResult.ClusterSize[tar_id]+=RareCluster.ClusterSize[seq_id];
                    RareAssignment[tar_id][count[tar_id]]=seq_id;
                    count[tar_id]++;
                    RareCluster.Assignment[seq_id]=tar_id;
                //}
                break;
            }  
        }
    }
    cout<<"dist done\n";
    MappedResult.ClusterMember=new string*[MappedResult.ClusterNumber];
    MappedResult.ClusterMemberFlag=new bool*[MappedResult.ClusterNumber];
    MappedResult.ClusterMemberSeq=new string*[MappedResult.ClusterNumber];
    for(int i=0;i<MappedResult.ClusterNumber;i++)
    {
        MappedResult.ClusterMember[i]=new string[MappedResult.ClusterIterationSize[i]];
        MappedResult.ClusterMemberFlag[i]=new bool[MappedResult.ClusterIterationSize[i]];
        memset(MappedResult.ClusterMemberFlag[i],false,sizeof(bool)*(MappedResult.ClusterIterationSize[i]));
        MappedResult.ClusterMemberSeq[i]=new string[MIN(20,MappedResult.ClusterSize[i])];
        for(int j=0;j<AbundantResult.ClusterIterationSize[i];j++)
            MappedResult.ClusterMember[i][j]=AbundantResult.ClusterMember[i][j];
        int kmin=MIN(20,AbundantResult.ClusterSize[i]), kmax=kmin;
        for(int j=AbundantResult.ClusterIterationSize[i];j<MappedResult.ClusterIterationSize[i];j++)
        {
            kmax=MIN(20,kmax+RareCluster.ClusterSize[RareAssignment[i][j-AbundantResult.ClusterIterationSize[i]]]);
            MappedResult.ClusterMember[i][j]=RareCluster.SeqID[RareAssignment[i][j-AbundantResult.ClusterIterationSize[i]]];
            for(int k=kmin;k<kmax;k++)
                MappedResult.ClusterMemberSeq[i][k]=RareCluster.EigenSeqs[RareAssignment[i][j-AbundantResult.ClusterIterationSize[i]]][k-kmin];
            kmin=kmax;
        }
        for(int j=0;j<MIN(20,AbundantResult.ClusterSize[i]);j++)
            MappedResult.ClusterMemberSeq[i][j]=AbundantResult.ClusterMemberSeq[i][j];
    }
    //output remaining rare clusters
    cout<<"output new rare\n";
    fstream RareOut;
    RareOut.open(fname.c_str(),ios_base::out);
    if(RareOut.fail()){cout<<"cannot create new rare file\n";}
    for(int i=0;i<RareCluster.SeqNum;i++)
    {
        if(RareCluster.Assignment[i]==-1)
        {
            RareOut<<">"<<RareCluster.SeqID[i]<<"\n";
            RareOut<<RareCluster.Seq[i]<<"\n";
            RareOut<<RareCluster.ClusterSize[i]<<"\n";
            RareOut<<RareCluster.ClusterStandardDeviation[i]<<"\n";
            for(int j=0;j<MIN(20,RareCluster.ClusterSize[i]);j++)
            {
                RareOut<<RareCluster.EigenSeqs[i][j]<<"\n";
            }
        }
    }
    RareOut.close();
    //clear
    for(int i=0;i<MappedResult.ClusterNumber;i++)
        delete[] RareAssignment[i];
    delete[] RareAssignment;
    delete[] count;
    TempBlockFileDelete(RareCluster);
    }//for if(RareCluster.SeqNum!=0)
    else //else (RareCluster.SeqNum==0)
        CopyClusteringResult(AbundantResult,MappedResult,true);
    cout<<"done\n";
    return;
}

void SummarizingResult(string oname, BayesianClusteringResult& MergeCenterResult, int Cluster, CROPParameters& CParameter)
{
    vector<string> TempFinalMembership;
    for(int p=0;p<MergeCenterResult.ClusterIterationSize[Cluster];p++)
        for(int j=0;j<CParameter.Block[CParameter.Maxiter];j++)
            for(int l=0;l<CParameter.BlockResult[CParameter.Maxiter][j].ClusterNumber;l++)
            {
                if(strcmp(MergeCenterResult.ClusterMember[Cluster][p].c_str(),CParameter.BlockResult[CParameter.Maxiter][j].ClusterCenterID[l])==0)
                {                                                                                                                                                                
                    for(int q=0;q<CParameter.BlockResult[CParameter.Maxiter][j].ClusterIterationSize[l];q++)
                    {
                        if(CParameter.BlockResult[CParameter.Maxiter][j].ClusterMemberFlag[l][q]==true)
                        {
                            cout<<"Redundant Sequence Detected: "<<CParameter.BlockResult[CParameter.Maxiter][j].ClusterMember[l][q]<<endl;
                        }
                        else
                        {
                            CParameter.BlockResult[CParameter.Maxiter][j].ClusterMemberFlag[l][q]=true;
                            TempFinalMembership.push_back(CParameter.BlockResult[CParameter.Maxiter][j].ClusterMember[l][q]);                                  
                        }
                    }//q
                }//if
            }//l
    int TempMaxiter=CParameter.Maxiter;
    while(TempMaxiter>0)
    {
        TempMaxiter--;
        int TempFinalSize=(int)TempFinalMembership.size();
        for(int p=0;p<TempFinalSize;p++)
            for(int j=0;j<CParameter.Block[TempMaxiter];j++)
                for(int l=0;l<CParameter.BlockResult[TempMaxiter][j].ClusterNumber;l++)
                {
                    if(strcmp(TempFinalMembership[p].c_str(),CParameter.BlockResult[TempMaxiter][j].ClusterCenterID[l])==0)
                    {
                        for(int q=0;q<CParameter.BlockResult[TempMaxiter][j].ClusterIterationSize[l];q++)
                        {
                            if(strcmp(CParameter.BlockResult[TempMaxiter][j].ClusterMember[l][q].c_str(),CParameter.BlockResult[TempMaxiter][j].ClusterCenterID[l])!=0)
                            {
                                if(CParameter.BlockResult[TempMaxiter][j].ClusterMemberFlag[l][q]==true)
                                {
                                    cout<<"Redundant Sequence Detected: "<<CParameter.BlockResult[TempMaxiter][j].ClusterMember[l][q]<<endl;
                                }
                                else
                                {
                                    CParameter.BlockResult[TempMaxiter][j].ClusterMemberFlag[l][q]=true;
                                    TempFinalMembership.push_back(CParameter.BlockResult[TempMaxiter][j].ClusterMember[l][q]);
                                }
                            }//if
                        }//q
                    }//if
                }//l
    }//while   
    //output tempfinalmembership
    fstream finalout;
    if(Cluster==0)    
        finalout.open((oname+".cluster.list").c_str(), ios_base::out);
    else
        finalout.open((oname+".cluster.list").c_str(), ios_base::out|ios_base::app);
    if(finalout.fail()){cout<<"cannot open output file 3"<<endl; return;}
    finalout<<CParameter.FinalCenterID[Cluster]<<"\t";
    for(int j=0;j<(int)TempFinalMembership.size();j++)
    {
        finalout<<CParameter.List[TempFinalMembership[j]];
        if(j==(int)TempFinalMembership.size()-1)
            finalout<<"\n";
        else
            finalout<<",";
    }
    finalout.close();
    TempFinalMembership.clear();
    //**************************
    return;
}

void MergeUsingBC(const char* fname, string oname, CROPParameters &CParameter)
{
    cout<<"MergingUsingBC"<<endl;
    string TempCenters(fname);
    TempCenters+=".TempCenters";
    BayesianClusteringResult MergeCenterResult, AbundantResult, RareResult, MappedResult;
    Parameters Parameter;
    BayesianClustering(TempCenters.c_str(),MIN(CParameter.Finalk*10,LIMIT/5), CParameter.Lower, CParameter.Upper, AbundantResult, 20, Parameter);
    cout<<"clustering abundant cluster complete"<<endl;
    if(AbundantResult.ClusterNumber!=0)
        MappingRareToAbundant((TempCenters+".Rare"),AbundantResult, MappedResult, CParameter.Upper);
    BayesianClustering((TempCenters+".Rare").c_str(),MIN(CParameter.Finalk*10,LIMIT/5), CParameter.Lower, CParameter.Upper, RareResult, 20, Parameter); 
    if(RareResult.ClusterNumber!=0&&AbundantResult.ClusterNumber!=0)
    {
        MergingClusteringResult(MappedResult, RareResult, MergeCenterResult);
        ClearClusteringResult(RareResult);
        ClearClusteringResult(MappedResult);
        ClearClusteringResult(AbundantResult);
    }
    else if(RareResult.ClusterNumber==0)
    {
        CopyClusteringResult(MappedResult, MergeCenterResult, true);
        ClearClusteringResult(MappedResult);
        ClearClusteringResult(AbundantResult);
    }
    else if(AbundantResult.ClusterNumber==0)
    {
        CopyClusteringResult(RareResult, MergeCenterResult, true);
        ClearClusteringResult(RareResult);
    }  
    CParameter.Finalk=MergeCenterResult.ClusterNumber;
    CParameter.FinalCenterID=new string[CParameter.Finalk];
    CParameter.FinalCenterSeq=new string[CParameter.Finalk];
    CParameter.FinalSize=new int[CParameter.Finalk];
    CParameter.FinalSD=new double[CParameter.Finalk];
    memset(CParameter.FinalSize,0,sizeof(int)*CParameter.Finalk);
    
    for(int i=0;i<CParameter.Finalk;i++)
    {
        CParameter.FinalCenterID[i].assign(MergeCenterResult.ClusterCenterID[i]);
        CParameter.FinalCenterSeq[i].assign(MergeCenterResult.ClusterCenterSeq[i]);
        CParameter.FinalSD[i]=MergeCenterResult.ClusterStandardDeviation[i];
        CParameter.FinalSize[i]=MergeCenterResult.ClusterSize[i];
        SummarizingResult(oname, MergeCenterResult, i, CParameter);
    }//i
    /*
    for(int i=0;i<Block[Maxiter];i++)
        for(int j=0;j<BlockResult[Maxiter][i].ClusterNumber;j++)
            for(int l=0;l<BlockResult[Maxiter][i].ClusterIterationSize[j];l++)
                if(BlockResult[Maxiter][i].ClusterMemberFlag[j][l]==false)
                    cout<<"Missing Sequence Detected: "<<BlockResult[Maxiter][i].ClusterMember[j][l]<<endl;
    */       
   //*********Free Memory*****************
   ClearClusteringResult(MergeCenterResult);
   //**************************************
                                              
   return;                      
}

void InterResultRead(const char* fname, int iter, CROPParameters &CParameter)
{
    string TempCenters(fname);
    TempCenters+=".TempCenters";
    fstream BlockCenters;
    CParameter.Finalk=0;
    BlockCenters.open(TempCenters.c_str(),ios_base::out);
    if(BlockCenters.fail())
    {
        cout<<"Fail to Open Temp File"<<endl;
        return;
    }
    int totalcluster=0;
    for(int i=0;i<CParameter.Block[iter];i++) //for each block in this iteration
        totalcluster+=CParameter.BlockResult[iter][i].ClusterNumber; //count its cluster number
    int* Shuffle;
    Shuffle=new int[totalcluster];
    int Count=0;
    for(int i=0;i<CParameter.Block[iter];i++) //for each block in this iteration
        for(int j=0;j<CParameter.BlockResult[iter][i].ClusterNumber;j++) //for each cluster in one block
        {
            Shuffle[Count]=i;
            Count++;
        }
    CParameter.r = gsl_rng_alloc (gsl_rng_taus);
    gsl_rng_set(CParameter.r,(unsigned int)time((time_t *)NULL));
    gsl_ran_shuffle(CParameter.r,Shuffle,totalcluster,sizeof(int)); //shuffle the order
    int* TempCount=new int[CParameter.Block[iter]];
    memset(TempCount,0,sizeof(int)*CParameter.Block[iter]);
    for(int i=0;i<totalcluster;i++)
    {
        BlockCenters<<">"<<CParameter.BlockResult[iter][Shuffle[i]].ClusterCenterID[TempCount[Shuffle[i]]]<<"\n";
        BlockCenters<<CParameter.BlockResult[iter][Shuffle[i]].ClusterCenterSeq[TempCount[Shuffle[i]]]<<"\n";
        BlockCenters<<CParameter.BlockResult[iter][Shuffle[i]].ClusterSize[TempCount[Shuffle[i]]]<<"\n";
        BlockCenters<<CParameter.BlockResult[iter][Shuffle[i]].ClusterStandardDeviation[TempCount[Shuffle[i]]]<<"\n";
        for(int l=0;l<MIN(20,CParameter.BlockResult[iter][Shuffle[i]].ClusterSize[TempCount[Shuffle[i]]]);l++)
            BlockCenters<<CParameter.BlockResult[iter][Shuffle[i]].ClusterMemberSeq[TempCount[Shuffle[i]]][l]<<"\n";//l
        CParameter.Finalk++;        
        TempCount[Shuffle[i]]++;
    }//i
    BlockCenters.close();
    delete[] Shuffle;
    delete[] TempCount;
    gsl_rng_free(CParameter.r);
    cout<<"Number of Clusters in Previous Round = "<<CParameter.Finalk<<endl;
    return;
}

void InterResultReadForBC(string fname, int iter, CROPParameters &CParameter)
{
    fstream Abundant, Rare;
    Abundant.open((fname+".TempCenters").c_str(),ios_base::out);
    if(Abundant.fail()){cout<<"cannot create abundant temp file\n";}
    Rare.open((fname+".TempCenters"+".Rare").c_str(),ios_base::out);
    if(Rare.fail()){cout<<"cannot create rare temp file\n";}
    for(int i=0;i<CParameter.Block[iter];i++) //for each block in this iteration
        for(int j=0;j<CParameter.BlockResult[iter][i].ClusterNumber;j++)
        {
            if(CParameter.BlockResult[iter][i].ClusterSize[j]>CParameter.RareCut)
            {
                Abundant<<">"<<CParameter.BlockResult[iter][i].ClusterCenterID[j]<<"\n";
                Abundant<<CParameter.BlockResult[iter][i].ClusterCenterSeq[j]<<"\n";
                Abundant<<CParameter.BlockResult[iter][i].ClusterSize[j]<<"\n";
                Abundant<<CParameter.BlockResult[iter][i].ClusterStandardDeviation[j]<<"\n";
                for(int l=0;l<MIN(20,CParameter.BlockResult[iter][i].ClusterSize[j]);l++)
                    Abundant<<CParameter.BlockResult[iter][i].ClusterMemberSeq[j][l]<<"\n";
            }//if
            else
            {
                Rare<<">"<<CParameter.BlockResult[iter][i].ClusterCenterID[j]<<"\n";
                Rare<<CParameter.BlockResult[iter][i].ClusterCenterSeq[j]<<"\n";
                Rare<<CParameter.BlockResult[iter][i].ClusterSize[j]<<"\n";
                Rare<<CParameter.BlockResult[iter][i].ClusterStandardDeviation[j]<<"\n";
                for(int l=0;l<MIN(20,CParameter.BlockResult[iter][i].ClusterSize[j]);l++)
                    Rare<<CParameter.BlockResult[iter][i].ClusterMemberSeq[j][l]<<"\n";               
            }//else
        }
    Abundant.close();
    Rare.close();
    return;
}

int CROP(string fname, string oname, int TempBlockSize, int TempStep, int TempCtrl, float Lower, float Upper, int TempMaxiter, int TempMaxiterBlockSize, int tRareCut, CROPParameters &CParameter)
{
    time_t start,stop;
    fstream finalout;
    fstream logfile;
    string logfilename;
    logfilename=oname+".log";
    logfile.open(logfilename.c_str(), ios_base::out);
    if(logfile.fail())
    { 
        cout<<"cannot create log file"<<endl;
        return 0;
    }
	start=time(NULL);
	InitCROPParameters(CParameter);
	CParameter.Block[0]=TempBlockSize;
	CParameter.step=TempStep;
    switch(TempCtrl)
    {
        case 0:
             CParameter.Lower=1.5; 
             CParameter.Upper=2.5;
             break;
        case 1:
             CParameter.Lower=1;
             CParameter.Upper=1.5;
             break;
        case 2:
             CParameter.Lower=Lower;
             CParameter.Upper=Upper;
             break;
        default:
             cout<<"no such threhold options, -g will be used"<<endl;
             logfile<<"-g is used\n";
             CParameter.Lower=1.5; 
             CParameter.Upper=2.5;
             break;
    }
	CParameter.Maxiter=TempMaxiter;
	CParameter.RareCut=tRareCut;
	CParameter.MaxiterBlockSize=TempMaxiterBlockSize;
	CParameter.BlockResult=new BayesianClusteringResult*[CParameter.Maxiter];
	CParameter.BlockResult[0]=new BayesianClusteringResult[CParameter.Block[0]];
	CParameter.SeqNum=LoadingList(fname.c_str(),CParameter);
	logfile<<"-i "<<fname<<" -b "<<TempBlockSize<<" -e "<<TempStep<<" -m "<<TempMaxiter<<" -z "<<TempMaxiterBlockSize<<" Threshold: "<<CParameter.Lower<<" "<<CParameter.Upper<<endl;
	logfile<<CParameter.SeqNum<<endl;
    if(!ChompFile(fname.c_str(), CParameter.Block[0], CParameter))
    {
        cout<<"chop file fail"<<endl;
        logfile<<"chop file fail"<<endl;
        return 0;
    }
    #pragma omp parallel for
    for(int i=0;i<CParameter.Block[0];i++)
    {
        char Buf[ReadSize];
        memset(Buf, 0, sizeof(char)*ReadSize);
        Parameters Parameter;
        itoa(i+1,Buf,10);
        string BlockFnameNum(Buf);
        string BlockFnameString=fname;
        BlockFnameString+=".BlockSeq"+BlockFnameNum;
        #pragma omp critical
        {
            cout<<"Round: 0 Block: "<<i+1<<endl;
        }//omp critical
        BayesianClustering(BlockFnameString.c_str(), CParameter.step, CParameter.Lower, CParameter.Upper, CParameter.BlockResult[0][i],0,Parameter);
        #pragma omp critical
        {
            logfile<<"Round: 0 Block: "<<i+1<<": "<<CParameter.BlockResult[0][i].ClusterNumber<<endl;
        }
    } 
    for(int i=1;i<CParameter.Maxiter;i++)
    {
        InterResultRead(fname.c_str(), i-1, CParameter);
        if((CParameter.Finalk>=(double)0.9*CParameter.PrevIterk)&&(i>1))
        {
            CParameter.Maxiter=i-1;                 
            break;
        }//if
        if(CParameter.Finalk>CParameter.MaxiterBlockSize)
        {
            //2010-2-7 if merging is not efficient, try to increase the block size
            //********************************************************************
            if(CParameter.Finalk>=(double)0.8*CParameter.PrevIterk)
            {
                CParameter.PrevIterk=CParameter.Finalk;
                CParameter.Block[i]=MAX(CParameter.Finalk/MIN(1.5*CParameter.MaxiterBlockSize, LIMITB)+1,omp_get_num_procs());
            }
            else
            {
                CParameter.PrevIterk=CParameter.Finalk;
                CParameter.Block[i]=MAX(CParameter.Finalk/CParameter.MaxiterBlockSize+1,omp_get_num_procs());
            }
            //********************************************************************
            //2010-2-7 done
            if(!ChompInterFile(fname.c_str(), CParameter.Block[i], CParameter))
            {
                cout<<"Chomp Fail!"<<endl;
                logfile<<"Chomp Fail!"<<endl;
                return 0;
            }
            CParameter.BlockResult[i]=new BayesianClusteringResult[CParameter.Block[i]];
            #pragma omp parallel for
            for(int j=0;j<CParameter.Block[i];j++)
            {
                char Buf[ReadSize];
                Parameters Parameter;
                itoa(j+1,Buf,10);
                string BlockFnameNum(Buf);
                string BlockFnameString=fname;
                BlockFnameString+=".BlockSeq"+BlockFnameNum;
                #pragma omp critical
                {
                    cout<<"Round: "<<i<<" Block: "<<j+1<<endl;
                }//omp critical
                BayesianClustering(BlockFnameString.c_str(), CParameter.step, CParameter.Lower, CParameter.Upper, CParameter.BlockResult[i][j],i,Parameter);
                #pragma omp critical
                {
                    logfile<<"Round: "<<i<<" Block: "<<j+1<<": "<<CParameter.BlockResult[i][j].ClusterNumber<<endl;
                }
            }//j
            if(i==CParameter.Maxiter-1)
            {
                InterResultRead(fname.c_str(), i, CParameter);
                CParameter.Maxiter=i;
            }
        }//if
        else
        {
            CParameter.Maxiter=i-1;
            break;
        }//else
    }//i
    //finalout<<endl;
    //finalout<<"Merging Using BC"<<endl;
    logfile<<"Merging Using BC"<<endl;
    InterResultReadForBC(fname, CParameter.Maxiter, CParameter);
    MergeUsingBC(fname.c_str(), oname, CParameter);
    logfile<<"Writing Output"<<endl;
    logfile.close();
    /*
    for(int i=0;i<Block;i++)
        MergeResult(i);
    */       
    //***************output result***************************************
    string ofilename;
    ofilename=oname+".cluster";
    finalout.open(ofilename.c_str(), ios_base::out);
    if(finalout.fail())
    { 
        cout<<"cannot open output file 1"<<endl;
        return 0;
    }
    finalout<<CParameter.Finalk<<endl;
    for(int i=0;i<CParameter.Finalk;i++)
    {
        finalout<<CParameter.FinalCenterID[i]<<"\t"<<CParameter.FinalSize[i]<<"\t"<<CParameter.FinalSD[i]<<"\n";
    }
    finalout.close();
    ofilename=oname+".cluster.fasta";
    finalout.open(ofilename.c_str(), ios_base::out);
    if(finalout.fail())
    { 
        cout<<"cannot open output file 2"<<endl;
        return 0;
    }
    for(int i=0;i<CParameter.Finalk;i++)
    {
        finalout<<">"<<CParameter.FinalCenterID[i]<<"\n";
        finalout<<CParameter.FinalCenterSeq[i]<<"\n";
    }
    finalout.close();
    //cout<<"output done"<<endl;
    //******************************************************************
    //**********Free Memory*********************************************
    DeleteTmpFiles(fname, CParameter);
    for(int i=0;i<CParameter.Maxiter+1;i++)
    {
        for(int j=0;j<CParameter.Block[i];j++)
        {
            for(int l=0; l<CParameter.BlockResult[i][j].ClusterNumber; l++)
            {
                delete[] CParameter.BlockResult[i][j].ClusterMember[l];
                delete[] CParameter.BlockResult[i][j].ClusterMemberFlag[l];
                delete[] CParameter.BlockResult[i][j].ClusterCenterID[l];
                delete[] CParameter.BlockResult[i][j].ClusterCenterSeq[l];
                delete[] CParameter.BlockResult[i][j].ClusterMemberSeq[l];
            }//l
            delete[] CParameter.BlockResult[i][j].ClusterMemberSeq;
            delete[] CParameter.BlockResult[i][j].ClusterMember;
            delete[] CParameter.BlockResult[i][j].ClusterMemberFlag;
            delete[] CParameter.BlockResult[i][j].ClusterCenterID;
            delete[] CParameter.BlockResult[i][j].ClusterCenterSeq;
            delete[] CParameter.BlockResult[i][j].ClusterStandardDeviation;
            delete[] CParameter.BlockResult[i][j].ClusterSize;
            delete[] CParameter.BlockResult[i][j].ClusterIterationSize;
        }//j
        delete[] CParameter.BlockResult[i];
    }
    delete[] CParameter.BlockResult;
    delete[] CParameter.FinalCenterID;
    delete[] CParameter.FinalCenterSeq;
    delete[] CParameter.FinalSD;
    delete[] CParameter.FinalSize;
    CParameter.List.empty();
    //**********************************************
    stop=time(NULL);
    cout<<"Time Cost Is "<<stop-start<<" Sec"<<endl;
    return 1;    
}

void DeleteTmpFiles(string fname, CROPParameters &CParameter)
{
    char Buffer[10];
    for(int i=1;i<=CParameter.Block[0];i++)
    {
        string tempNum(itoa(i, Buffer, 10));
        string tempFile=fname+".BlockSeq"+tempNum;
        if(remove(tempFile.c_str())!=0)
            cout<<"Error deleting "<<tempFile<<endl;
    }
    string tempFile=fname+".TempCenters";
    if(remove(tempFile.c_str())!=0)
        cout<<"Error deleting "<<tempFile<<endl;
    return;
} 
