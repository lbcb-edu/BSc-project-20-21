#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <fstream>
#include "fasta_parser.hpp"
#include "fastq_parser.hpp"
#include "sequence.hpp"
#include "Pink_mapperConfig.h"
using namespace std;

static int help_flag;
static int version_flag;

void printHelpMessage(){
    	string line;
  	ifstream myfile ("help.txt");
  	if (myfile.is_open()){
    		while ( getline (myfile,line) ){
      			cout << line << '\n';
    		}
    		myfile.close();
  	}

  	else cout << "Unable to open file"; 

  	exit(0);
}

void swap(int *xp, int *yp){  
	int temp = *xp;  
	*xp = *yp;  
	*yp = temp;  
}  
  
int shellSort(int arr[], int n){   
	for (int gap = n/2; gap > 0; gap /= 2){        
		for (int i = gap; i < n; i += 1){            
            		int temp = arr[i]; 
            		int j;             
            		for (j = i; j >= gap && arr[j - gap] < temp; j -= gap) arr[j] = arr[j - gap];               
            		arr[j] = temp; 
        	} 
    	} 
    	return 0; 
}

//fragments sum length
int total_length(const vector<unique_ptr<biosoup::Sequence>> &list){
	int sum;
	for (int i = 0; i<list.size(); i++){
		int len = list[i] -> data.length();
		sum += len;
	}
	return sum;
}

//N50
void N50_length(const vector<unique_ptr<biosoup::Sequence>> &list){
	int sum;
	int i = 0;
	int total = total_length(list);
	
	int * len_array;
	len_array = new (nothrow) int [list.size()];
	if (len_array == nullptr) throw bad_alloc();
	
	for (int i = 0; i<list.size(); i++){
		len_array[i] = list[i] -> data.length();
	}
	
	shellSort(len_array, list.size());
	
	while(sum < total/2){
		sum += len_array[i];
		i += 1;
	}
	int N50 = len_array[i-1];
	cerr << "N50 fragments length: " << N50 << "\n";
	delete len_array;
	return;
}


//refrence genome name and length
void ref_genome(const vector<unique_ptr<biosoup::Sequence>> &list){
	for (int i = 0; i<list.size(); i++){
		cerr << i+1 << ". refrence genome name: " << list[i] -> name << "\n"
		<< i+1 << ". refrence genome length: " << list[i] -> data.length() << "\n";
	}
}

//number of fragments
inline void frag_number(const vector<unique_ptr<biosoup::Sequence>> &list){
	cerr << "Number of refrences in the fragments file: " << list.size() << "\n";
}

//avg, min, max fragment length
void frag_stats(const vector<unique_ptr<biosoup::Sequence>> &list){
	float avg;
	float sum = 0;
	int min = list[1] -> data.length();
	int max = list[1] -> data.length();
	for (int i = 0; i<list.size(); i++){
		int len = list[i] -> data.length();
		sum += len;
		if (len < min) min = len;
		if (len > max) max = len;
	}
	
	avg = sum/((float) list.size());
	cerr << "Average length of a fragment: " << avg << "\n";
	cerr << "Maximal length of a fragment: " << max << "\n";
	cerr << "Minimal length of a fragment: " << min << "\n";
}

int main(int argc, char **argv){
    int c;
    while(1){
        static struct option long_options[] =
            {
                {"help", no_argument, &help_flag, 1},
                {"version", no_argument, &version_flag, 1},
                {0, 0, 0, 0}
            };
        int option_index = 0;
        c = getopt_long(argc, argv, "hv", long_options, &option_index);
        if(c == -1){
            break;
        }
        switch (c)
        {
        case 0:
            break;
        case 'v':
            version_flag = 1;
            break;
        case 'h':
            help_flag = 1;
            break;
        case '?':
            break;
        default:
            abort();
        }
    }
    if(help_flag) printHelpMessage();
    if(version_flag){
        cout << "v" << pink_mapper_VERSION_MAJOR << "." << pink_mapper_VERSION_MINOR << "." << pink_mapper_VERSION_PATCH << endl;
        exit(0);
    }
    if(optind < argc){
        auto p = bioparser::Parser<biosoup::Sequence>::Create<bioparser::FastaParser>(argv[1]);

        // parse whole file
        auto s = p->Parse(-1);

        ref_genome(s);
        
        p = bioparser::Parser<biosoup::Sequence>::Create<bioparser::FastaParser>(argv[2]);

        // parse in chunks
        vector<unique_ptr<biosoup::Sequence>> ss;
        uint32_t chunk_size = 500 * 1024 * 1024;  // 500 MB
        for (auto t = p->Parse(chunk_size); !t.empty(); t = p->Parse(chunk_size)) {
            ss.insert(
                ss.end(),
                make_move_iterator(t.begin()),
                make_move_iterator(t.end()));
        }
        
        frag_number(ss);
        frag_stats(ss);
        N50_length(ss);
    }
    return 0;
}
