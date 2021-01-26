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
#include "pink_demo.hpp"
#include <omp.h>
#define FRAG 4

using namespace std;
using namespace pink;

//longest increasing subsequence
pair<unsigned int,unsigned int> LIS(const vector<tuple<unsigned int, unsigned int, bool>> &matches, int length, int flength){
	vector<int> P(matches.size(),0);
	vector<int> M(matches.size()+1,0);
	int L=0;

	for(int i=0;i<matches.size();i++){
		int lo = 1;
		int hi = L;
		while(lo<=hi){
			int mid = ceil((lo+hi)/2);
			if(get<0>(matches[M[mid]]) < get<0>(matches[i])) lo = mid+1;
			else hi = mid-1;
		}
		int newL = lo;

		P[i] = M[newL-1];
		M[newL]=i;

		if(newL>L) L=newL;

	}
	vector<vector<tuple<unsigned int, unsigned int, bool>>> matrix;
	vector<tuple<unsigned int, unsigned int, bool>> S;
	pair<unsigned int,unsigned int> region;
	int k = M[L];
	for(int i=L-1;i>-1;i--){
		if(i<(L-1)){
			if((get<1>(S[0])-get<1>(matches[k]))>(flength)){
				matrix.push_back(S);
				S.clear();
				S.shrink_to_fit();
			}
		}
		S.insert(S.begin(),matches[k]);
		k=P[k];
	}
	if(matrix.size() ==0) matrix.push_back(S);;
	int max=0;
	for(int i=1;i<matrix.size();i++){
		if(matrix[i].size()>matrix[max].size()) max = i;
	}
	S = matrix[max];
	region.first = get<1>(S[0]);
	region.second = get<1>(S.back())+length+1;
	return region;
}

int main(int argc, char **argv){
    
    int match = 0, mismatch = 1, gap = 1;
    float frequency = 0.001;
	unsigned int klength = 15, wlength = 5;  
    string type = "global";
	int tnum = 1;

	//parse program options and their arguments
    parse_opt(argc, argv,
				match, mismatch, gap, frequency,
				klength, wlength, type, tnum);
	
	AlignmentType a_type;
	if(type == "global") a_type = global;
	else if(type == "local") a_type = local;
	else if(type == "semiglobal") a_type = semiglobal;
	else { cerr << "Error: Alignment type " << type << " does not exist.\n"
		<< "Try typing 'pink_mapper --help for additional help\n";
		exit(1);
	}

	//main function
    if(optind < argc){
		//read input file paths
		string stfile = argv[argc-2];
    	string ndfile = argv[argc-1];

		// parse files
        auto p = bioparser::Parser<biosoup::Sequence>::Create<bioparser::FastaParser>(stfile);
        auto s = p->Parse(-1);
		p = bioparser::Parser<biosoup::Sequence>::Create<bioparser::FastaParser>(ndfile);
		auto ss = p->Parse(-1);
        
		//print some basic input file statistics
		ref_genome(s);
		frag_number(ss);
        frag_stats(ss);
        N50_length(ss);

		//find refrence k-mers
		vector<tuple<unsigned int,unsigned int, bool>> ref_kmers;
		char seq[s[0] -> data.length() + 1];
		strcpy(seq, (s[0] -> data).c_str());
		ref_kmers = pink::Minimize(seq, strlen(seq), klength, wlength);

		//calculate number of occurences for each kmer found in refrence sequence
		unordered_map<unsigned int,unsigned int> dist_mins;
		unordered_map<unsigned int, unsigned int>::iterator it;
		for(auto min: ref_kmers){
			if(dist_mins.find(get<0>(min))==dist_mins.end()) dist_mins.insert({get<0>(min),1});
			else dist_mins[get<0>(min)]++;
		}

		//cutoff frequent minimizers
		double one = 1.0/(s[0] -> data.length());
		double thresh = 0;
		set<int> junk;
		vector<pair<unsigned int,unsigned int>> A;
		for (auto& it: dist_mins){
			A.push_back(it);
		}
		vector<pair<unsigned int,unsigned int>>::iterator it_pair = A.begin();
		sort(A.begin(),A.end(),cmp);
		while(thresh<frequency){
			thresh += (it_pair -> second)*one;
			junk.insert(it_pair -> first);
			it_pair++;
		}
		vector<tuple<unsigned int,unsigned int, bool>>::iterator iter;
		vector<tuple<unsigned int,unsigned int, bool>> new_kmers;

		for(auto m: ref_kmers){
			if(junk.find(get<0>(m))==junk.end()){
				new_kmers.push_back(m);
			}
		}

		//find selected read k-mers
		#pragma omp parallel for num_threads(tnum)
		for(int glava=0;glava<(ss.size());glava++){
			//if(ss[glava] -> data.length() >= 25000) continue;
			vector<tuple<unsigned int,unsigned int, bool>> frag_kmers;
			set<tuple<unsigned int, bool>>::iterator ite;
			set<tuple<unsigned int, bool>> frag_set;
			vector<tuple<unsigned int,unsigned int, bool>> ref_match;
			vector<tuple<unsigned int,unsigned int, bool>> frag_match;
			vector<tuple<unsigned int,unsigned int, bool>> ref_match_comp;
			vector<tuple<unsigned int,unsigned int, bool>> frag_match_comp;
			char* target = (char*) malloc(1*sizeof(char));
			char* query = (char*) malloc(1*sizeof(char));
			char seqq[ss[glava] -> data.length() + 1];
			strcpy(seqq, (ss[glava] -> data).c_str());
			frag_kmers = pink::Minimize(seqq, strlen(seqq), klength, wlength);
			//print read length

			//find matching kmers
			for(auto kmer: frag_kmers){
				frag_set.insert(make_tuple(get<0>(kmer), get<2>(kmer)));
			}
			set<tuple<unsigned int, bool>> ref_set;
			for(auto kmer: ref_kmers){
				if(get<2>(kmer)==true){
					ref_set.insert(make_tuple(get<0>(kmer), get<2>(kmer)));
					ite = frag_set.find(make_tuple(get<0>(kmer), get<2>(kmer)));
					if(ite != frag_set.end()){
						ref_match.push_back(make_tuple(get<0>(kmer), get<1>(kmer), get<2>(kmer)));
					}
					ite = frag_set.find(make_tuple(get<0>(kmer), false));
					if(ite != frag_set.end()){
						ref_match_comp.push_back(make_tuple(get<0>(kmer), get<1>(kmer), get<2>(kmer)));
					}
				}
			}
			for(auto kmer: frag_kmers){
				if(ref_set.find(make_tuple(get<0>(kmer), true)) != ref_set.end() && get<2>(kmer)==true){
					frag_match.push_back(make_tuple(get<0>(kmer), get<1>(kmer), get<2>(kmer)));
				}
				if(ref_set.find(make_tuple(get<0>(kmer), true)) != ref_set.end() && get<2>(kmer)==false){
					frag_match_comp.push_back(make_tuple(get<0>(kmer), get<1>(kmer), get<2>(kmer)));
				}
			}

			//sort matching kmers by their position in sequences in
			bool comp = ref_match.size() < ref_match_comp.size();
			shellSort<tuple<unsigned int, unsigned int, bool>>((comp ? ref_match_comp : ref_match),1);
			shellSort<tuple<unsigned int, unsigned int, bool>>((comp ? frag_match_comp : frag_match),1);

			//minimize all reads from fragment file (uncomment at your own risk)
			//call_minimizers(frequency,klength,wlength,s,ss);

			//longest increasing subsequence on sorted matches
			pair<unsigned int,unsigned int> target_interval = LIS((comp ? ref_match_comp : ref_match), klength, ss[glava] -> data.length());
			pair<unsigned int,unsigned int> query_interval = LIS((comp ? frag_match_comp : frag_match), klength, ss[glava] -> data.length());

			//prepare alignment candidates for alignment (convert strings to character arrays)

			target = (char*) realloc(target,((s[0] -> data.substr(target_interval.first-1,
									target_interval.second-target_interval.first)).length()+1)*sizeof(char));
			
			strcpy(target, s[0] -> data.substr(target_interval.first-1,
											target_interval.second-target_interval.first).c_str());
			
			query = (char*) realloc(query,((ss[glava] -> data.substr(query_interval.first-1,
									query_interval.second-query_interval.first)).length()+1)*sizeof(char));

			
			strcpy(query, ss[glava] -> data.substr(query_interval.first-1,
											query_interval.second-query_interval.first).c_str());

			//align found regions
			string cigar;
			cigar = "";
			int result = Align(query, (query_interval.second-query_interval.first),
				target, (target_interval.second-target_interval.first),
				a_type,
				match, mismatch, gap,
				&cigar);

			//count number of matches and total number of matches mismatches insertions and deletions
			int sum_cnt = 0;
			int match_cnt = 0;
			string tmp;
			for(int i=0;i<cigar.length();i++){
				if(isdigit(cigar[i])){
					tmp = tmp + cigar[i];
				}
				else{
					sum_cnt += stoi(tmp);
					if(cigar[i]=='='){
						match_cnt += stoi(tmp);
					}
					tmp = "";
				}
			}
			
			//PAF format
			cout << ss[glava] -> name << "    "
					<< (query_interval.second-query_interval.first) << "    "
					<< query_interval.first << "    "
					<< query_interval.second<< "    "
					<< (comp ? "-" : "+") << "    "
					<< s[0] -> name << "    "
					<< (target_interval.second-target_interval.first) << "    "
					<< target_interval.first << "    "
					<< target_interval.second<< "    "
					<< match_cnt << "    "
					<< sum_cnt << "    "
					<< "255" << "    "
					<< (cigar_flag ? ("cg:Z:" + cigar) : "" )<< endl;
			free(target);
			free(query);
		}
		
    }
    return 0;
}
