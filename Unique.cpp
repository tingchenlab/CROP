#include "Unique.h"

using namespace std;

int ExtractUnique(const char* fname)
{
    map<string,UniqueSeqs> Data;
	string Sequence, SequenceID, temp;
	unsigned int i = 0, len = 0, lentotal = 0, num_seq=0, unique_num_seq=0;
	fstream InputSeq;
	InputSeq.open(fname);
	if(InputSeq.fail())
	{
	   cout<<"Cannot open "<<fname<<endl;
	   return 0;
    }
    getline(InputSeq,temp);	
	while(!InputSeq.eof())    
    {
        if(temp[0]=='>')
        {
            num_seq++;
            SequenceID=temp;
            size_t found=SequenceID.find_first_of(" ,\t\r");
            SequenceID=SequenceID.substr(1,found-1);
            temp.clear();
            getline(InputSeq,temp);       
			while(temp[0]!='>')
            {
                if(temp[0]!='\n'&&temp[0]!='\0'&&temp[0]!='\r')
                {
                    size_t found=temp.find_first_of("\r\n");
                    Sequence+=temp.substr(0,found-1);
                }
                if(InputSeq.eof())                          
                    break;             
                getline(InputSeq,temp);
            }
            StringUpper(Sequence); //transform to uppercase
            map<string,UniqueSeqs>::iterator it = Data.find(Sequence);
		    if (it == Data.end())
            { 	//it's unique.
		        Data[Sequence].Name = SequenceID;
                Data[Sequence].Sequences = Sequence;  
		        Data[Sequence].Members = SequenceID;
		        Data[Sequence].Size = 1;
		        unique_num_seq++;
            }
            else 
            {  // its a duplicate.
                Data[Sequence].Members += "," + SequenceID;
	            Data[Sequence].Size++;
		    }
		    Sequence.clear();
        }           
    }	
	cout<<"num_seq = "<<num_seq<<endl;	
	InputSeq.close();
	//output unique sequences
    fstream OutputSeq1,OutputSeq2;
    string temponame1(fname);
    string temponame2;
    temponame1+=".unique";
    temponame2=temponame1+".list";
    OutputSeq1.open(temponame1.c_str(), ios_base::out);
    OutputSeq2.open(temponame2.c_str(), ios_base::out);
    if(OutputSeq1.fail())
    {
        cout<<"Can not open output file "<<temponame1<<endl;
        return 0;
    }
    if(OutputSeq2.fail())
    {
        cout<<"Can not open output file "<<temponame2<<endl;
        return 0;
    }
    map<string,UniqueSeqs>::iterator it = Data.begin();
    OutputSeq2<<unique_num_seq<<endl;
    while(it!=Data.end())
    {
        OutputSeq1<<">"<<(*it).second.Name<<"\n";
        OutputSeq1<<(*it).second.Sequences<<"\n";
        OutputSeq1<<(*it).second.Size<<"\n";
        OutputSeq2<<(*it).second.Name<<"\n";
        OutputSeq2<<(*it).second.Members<<"\n";
        it++;
    }
    OutputSeq1.close();
    OutputSeq2.close();
    num_seq=Data.size();
    Data.empty();
    cout<<"extracting unique sequences done"<<endl;
	return num_seq;    
}

void StringUpper(std::string& s)
{
     std::string::iterator i = s.begin();
     std::string::iterator end = s.end();	 
     while (i != end) 
     {
        *i = std::toupper((unsigned char)*i);
        ++i;
     }
     return;
}
