#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string>
#include <fstream>
#include <unordered_map>
#include <set>
#include <cmath>
#include <tuple>
#include "fasta_parser.hpp"
#include "fastq_parser.hpp"
#include "sequence.hpp"
#include "pink_alignment.hpp"
#include "pink_minimizers.hpp"
#include "Pink_mapperConfig.h"

using namespace std;
using namespace pink;

static int help_flag;
static int version_flag;
static int align_flag;
static int match_flag;
static int mismatch_flag;
static int gap_flag;
static int frequency_flag;
static int klength_flag;
static int wlength_flag;
static int cigar_flag;
static int thread_flag;

inline bool cmp(pair<unsigned int, unsigned int>& a,
		pair<unsigned int, unsigned int>& b){
	return a.second > b.second;
}

void cutoff(unordered_map<unsigned int, unsigned int>& m, float f){
	vector<pair<unsigned int,unsigned int>> A;
	for (auto& it: m){
		A.push_back(it);
	}
	sort(A.begin(),A.end(),cmp);
	int compare = -1;
	float param = 1.0f/m.size();
	float count = 0;
	vector<pair<unsigned int,unsigned int>>::iterator it;
	for(auto& it: A){
		cerr << it.first << ' ' << it.second << endl;
	}
	for(it=A.begin();it!=A.end();it++){
		count += param;
		if(count>=f){
			do{
				it++;
			}while(it->second==prev(it)->second);
			cerr << "Number of occurences of the most frequent minimizer without top " << f
				<< " occuring minimizers: "  << (it) -> second << "\n";
			break;
		}
	}
}

void call_minimizers(float frequency, unsigned int klength, unsigned int wlength,
					const vector<unique_ptr<biosoup::Sequence>> &ref_list,
					const vector<unique_ptr<biosoup::Sequence>> &frag_list){
	vector<tuple<unsigned int,unsigned int, bool>> mins;
	unordered_map<unsigned int,unsigned int> dist_mins;
	set<unsigned int> cache_set;

	cerr << "Minimizing refrence genome file...\n";
	for(int i=0;i<ref_list.size();i++){
		char seq[ref_list[i] -> data.length() + 1];
		strcpy(seq, (ref_list[i] -> data).c_str());
		mins = pink::Minimize(seq, strlen(seq), klength, wlength);
		for(auto min: mins){
			if(dist_mins.find(get<0>(min))==dist_mins.end()) dist_mins.insert({get<0>(min),1});
			else dist_mins[get<0>(min)]++;
		}
	}
	unordered_map<unsigned int, unsigned int>::iterator it;
	set<unsigned int>::iterator it_set;
	cerr << "Number of distinctive k-mers: " << dist_mins.size() << "\n\n";
	cerr << "Fraction of singletons:\n";
	for(it=dist_mins.begin();it!=dist_mins.end();it++){
		if(it->second == 1) cerr << it->first <<"\n";
	}
	cutoff(dist_mins,frequency);
	dist_mins.clear();
	cache_set.clear();
	cerr << "\nMinimizing fragments genome file...\n";
	for(int i=0;i<frag_list.size();i++){
		//cerr << "Minimizing sequence: " << i << "\n"; 
		char seq[frag_list[i] -> data.length() + 1];
		strcpy(seq, (frag_list[i] -> data).c_str());
		mins = pink::Minimize(seq, sizeof(seq)-1, klength, wlength);
		for(auto min: mins){
			//cout << get<0>(min) << " " << get<1>(min) << " " << get<2>(min) << endl;
			if(cache_set.find(get<0>(min))==cache_set.end()){
				cache_set.insert(get<0>(min));
			}
		}
		for(it_set=cache_set.begin();it_set!=cache_set.end();it_set++){
			if(dist_mins.find(*it_set)==dist_mins.end()) dist_mins.insert({*it_set,1});
			else dist_mins[*it_set]++;
		}
	}
	cerr << "\nNumber of distinctive k-mers: " << dist_mins.size() << "\n\n";

	cerr << "Fraction of singletons:\n";
	for(it=dist_mins.begin();it!=dist_mins.end();it++){
		if(it->second == 1) cerr << it->first << "\n";	
	}
	cerr << "\n";
	cutoff(dist_mins,frequency);
	
}

int call_align(string type, int match, int mismatch, int gap, int size, const vector<unique_ptr<biosoup::Sequence>> &list){
	time_t t;
	srand((unsigned) time(&t));
	//assign alignment type 
	AlignmentType a_type;
	if(type == "global") a_type = global;
	else if(type == "local") a_type = local;
	else if(type == "semiglobal") a_type = semiglobal;
	else { cerr << "Error: Alignment type " << type << " does not exist.\n"
		<< "Try typing 'pink_mapper --help for additional help\n";
		exit(1);
	}

	//pick two random sequences to calculate alignment score
	int q_index;
	int t_index;
	do{
		q_index = (rand() % (size - 0 + 1));
	} while(list[q_index] -> data.length() > 5000);

	do{
		t_index = (rand() % (size - 0 + 1));
	} while(list[t_index] -> data.length() > 5000);

	char query[list[q_index] -> data.length() + 1];
	strcpy(query, (list[q_index] -> data).c_str());

	char target[list[t_index] -> data.length() + 1];
	strcpy(target, (list[t_index] -> data).c_str());

	cerr << "Aligning sequences:\n" << list[q_index] -> name << "\n" << list[t_index] -> name << "\n\n";
	cerr << "Alignmet type: " << type << "\n" << "Match: " << match << "\n"
		<< "Mismatch: " << mismatch << "\n" << "Gap: " << gap << "\n\n";
	int result = Align(query, list[q_index] -> data.length(),
			target, list[t_index] -> data.length(),
			a_type,
			match, mismatch, gap);

	cerr << "alignment score: " << result << "\n\n"; 
	return result;
}

void printHelpMessage(){
	char help[] = {
	#include "help.txt"
	};
	for(int i=0;true;i++){
		if(help[i] == '\0') break;
		cerr << help[i];
	}
  	exit(0);
}

inline bool operator> (const tuple<unsigned int, unsigned int, bool> &t1, const tuple<unsigned int, unsigned int, bool> &t2){
	return get<1>(t1) > get<1>(t2);
}

inline bool operator< (const tuple<unsigned int, unsigned int, bool> &t1, const tuple<unsigned int, unsigned int, bool> &t2){
	return get<1>(t1) < get<1>(t2);
}

template <class T>  
void shellSort(vector<T> &arr, bool asc){   
	for (int gap = arr.size()/2; gap > 0; gap /= 2){        
		for (int i = gap; i < arr.size(); i += 1){            
            		T temp = arr[i]; 
            		int j;             
            		for (j = i; j >= gap && (asc ? arr[j - gap] > temp : arr[j - gap] < temp); j -= gap) arr[j] = arr[j - gap];               
            		arr[j] = temp; 
        	} 
    	}
	return; 
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
	int sum = 0;
	int i = 0;
	int total = total_length(list);
	
	vector<int> len_array;
	//len_array = new (nothrow) int [list.size()];
	//if (len_array == nullptr) throw bad_alloc();
	
	for (int i = 0; i<list.size(); i++){
		len_array.push_back(list[i] -> data.length());
	}
	
	shellSort<int>(len_array,0);
	
	while(sum < total/2){
		sum += len_array[i];
		i += 1;
	}
	int N50 = len_array[i-1];
	cerr << "N50 fragments length: " << N50 << "\n\n";
	//delete len_array;
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

//parse program options and their arguments
void parse_opt(int argc, char **argv,
				int &match, int &mismatch, int &gap, float &frequency,
				unsigned int &klength, unsigned int &wlength, string &type, int &tnum){
	int c;
    while(1){
        static struct option long_options[] =
            {
                {"help", no_argument, &help_flag, 1},
                {"version", no_argument, &version_flag, 1},
                {"align", required_argument, &align_flag, 1},
                {"match", required_argument, &match_flag, 1},
                {"mismatch", required_argument, &mismatch_flag, 1},
                {"gap", required_argument, &gap_flag, 1},
                {"frequency", required_argument, &frequency_flag, 1},
                {"klength", required_argument, &klength_flag, 1},
                {"wlength", required_argument, &wlength_flag, 1},
				{"cigar", no_argument, &cigar_flag, 1},
				{"threads", required_argument, &thread_flag, 1},
                {0, 0, 0, 0}
            };
        int option_index = 0;
        c = getopt_long(argc, argv, "hvct:a:m:n:g:f:k:w:", long_options, &option_index);
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
        case 'a':
            align_flag = 1;
            type = optarg;
            break;
        case 'm':
            match = atoi(optarg);
            break;
        case 'n':
            mismatch = atoi(optarg);
            break;
        case 'g':
            gap = atoi(optarg);
            break;
        case 'f':
            frequency = atof(optarg);
            break;
        case 'k':
            klength = atoi(optarg);
            break;
        case 'w':
            wlength = atoi(optarg);
            break;
		case 'c':
			cigar_flag = 1;
			break;
		case 't':
			tnum = atoi(optarg);
			break;
        case '?':
            break;
        default:
            abort();
        }
    }
    if(help_flag) printHelpMessage();
    if(version_flag){
        cerr << "v" << pink_mapper_VERSION_MAJOR << "." << pink_mapper_VERSION_MINOR << "." << pink_mapper_VERSION_PATCH << endl;
        exit(0);
    }
	
	
    if(align_flag == 1){
		string ndfile = argv[argc-1];
     	auto pp = bioparser::Parser<biosoup::Sequence>::Create<bioparser::FastaParser>(ndfile);
		auto ss = pp->Parse(-1);
		int interval = ss.size();
		call_align(type, match, mismatch, gap, interval-1, ss);
    }
}