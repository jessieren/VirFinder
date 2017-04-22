#include <Rcpp.h>
//#include <tr1/unordered_map>
#include <unordered_map>
#include <vector>
//using namespace std::tr1;
using namespace std;
//
using namespace Rcpp;
//
//
////unsigned long power;//parameter, power = 4^(k-1)
//
////receptacle of kmer counts
unordered_map<unsigned long,unsigned long> HashTable;
////
//int seqlength=999999999;
//// seq: save read/genome(A C G T) in each line from file
//char seq[999999999];
//// inverse of seq
//char seq_inverse[999999999];
int ZI = 4;



void printFour(vector<unsigned long> four)
{
	cout << "print ";
	for(int it=0; it<four.size(); it++)
	{
		cout << four[it] << "," ;
	}
	cout << endl;
	
}


vector<int> ten2four(unsigned long ten, int k)
{
	vector<int> four (k,0);
	unsigned long tmp = ten;
	for(int currentPos = k-1; currentPos >=0; --currentPos)
	{
		four[currentPos]=tmp%ZI;
		tmp/=ZI;
	}
	//while(tmp>=(ZI-1)) {four[currentPos]=tmp%ZI;tmp/=ZI;currentPos--; }
	//four[currentPos] = tmp;
	return four;
}


int four2ten(vector<int> four, int k)
{
	unsigned long ten = 0;
	for(int currentPos=(k-1); currentPos >= 0; --currentPos)
	{
		int tmp = four[currentPos] * pow(ZI,(k-1 - currentPos));
		ten = ten + tmp;
		//cout << currentPos << " " << ten << endl;
	}
	return ten;
	
}


// find the reverse compilmentary of a word
vector<int> reverseFour(vector<int> Four)
{
	vector<int> reverseFour(Four.size(), 4);
	for(int revPos = 0; revPos < Four.size(); revPos++)
	{
		reverseFour[revPos] = 3 - Four[Four.size()- 1 - revPos];
	}
	return reverseFour;
	
}



unsigned long SeqKmerCountSingle(char* seqDNA, int k, unsigned long power)
{
	//..........................
	// int count: The length of char seq
	int count = 0;
	//..........................
	int i=0,j=0;
	unsigned long index = 0;
	unsigned long total = 0;//total number of the kmers counted in seq
	while(seqDNA[i])
	{
		//kmer in seq[i,i+1,i+2,...i+k-1] transfered to an index
		//current index = floor(previous index/4^(k-1))*4 + 0 or 1 or 2 or 3
		if(seqDNA[i]=='A'|| seqDNA[i] == 'a') {j++; }
		else if(seqDNA[i]=='C'|| seqDNA[i] == 'c') { j++; index++; }
		else if(seqDNA[i]=='G'|| seqDNA[i] == 'g') { j++; index+=2; }
		else if(seqDNA[i]=='T'|| seqDNA[i] == 't') { j++; index+=3; }
		else { j=0; index=0; }//If seq[i] is ambiguous, reset j and index
		
		if( j == k )
		{
			HashTable[index]++;
			total++;
			index %= power;// current index = floor(previous index/4^(k-1))
			j--;//the lengh of seq[i+1,i+2,...i+k-1]
		}
		index*=ZI;//current index = floor(previous index/4^(k-1))*4
		i++;
		count++;
	}
	
	
	return total;
}





// compute EuFeature using double strand
void loadToVector(int k, unsigned long total, vector<unsigned long>& kmerTen,  vector<double>& kmerCount)
{
	
	//int countWord = 0;
	for(unsigned long currentKmerTen = 0; currentKmerTen < pow(ZI, k); currentKmerTen++)
	{
		//cout << "test2.5" << endl;
		vector<int> currentKmerFour = ten2four(currentKmerTen, k);
		vector<int> currentKmerRevFour = reverseFour(currentKmerFour);
		
		unsigned long currentKmerRevTen = four2ten(currentKmerRevFour, currentKmerFour.size());
		//cout << currentKmerRevTen << endl;
		//printFour(currentKmerFour);
		//printFour(currentKmerFour);
		//printFour(currentKmerRevFour);
		//cout << "test2.6" << endl;
		
		if( currentKmerRevTen >= currentKmerTen )
		{
			kmerTen.push_back(currentKmerTen);
			kmerCount.push_back((HashTable[currentKmerTen] + HashTable[currentKmerRevTen])/double(2 * total));
			//countWord ++;
			
			//cout << currentKmerTen << "," << (HashTable[currentKmerTen] + HashTable[currentKmerRevTen])/double(2 * total) << endl;
			//featureOutput << currentKmerTen << "," << TransToReal(X_w) << endl;
		}
	}
	
	//printFour(kmerTen);
	//featureOutput.close();
	
	//return X_w;
}




// [[Rcpp::export]]
List countSeqFeatureCpp( CharacterVector RseqDNA,  int k) {
	
	// convert to C++ type
	//Rcout << "seq is " ;
	char seqDNAChar[RseqDNA.size()];
	for( int i=0; i < RseqDNA.size(); i++ ){
		//Rcout << RseqDNA(i);
		seqDNAChar[i] = Rcpp::as< char >(RseqDNA(i));
	}
	
	unsigned long power = 1; for( int i = 0; i < k-1; i++) power *= 4;
	HashTable.clear();

	// count kmer
	unsigned long total = SeqKmerCountSingle(seqDNAChar, k, power);
	
	// pair words and output count
	vector<unsigned long> kmerTen;
	vector<double> kmerCount;
	loadToVector(k, total, kmerTen, kmerCount);
	//Rcout << "\n total:" << total << endl;

	// convert to Rcpp type
	NumericVector RkmerTen(kmerTen.size());
	RkmerTen = kmerTen;
	NumericVector RkmerCount(kmerCount.size());
	RkmerCount = kmerCount;
	
	List ret;
	ret["kmerTen"] = kmerTen;
	ret["kmerCount"] = kmerCount;
	return ret;
	
}


