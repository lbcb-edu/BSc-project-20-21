#include <iostream>
#include <getopt.h>

#define VERSION "v0.1.0"

static int help_flag;
static int version_flag;

void printHelp(){
    std::cout << "\norange_mapper usage:\n"
    "./orange_mapper [OPTION]\n"
    "./orange_mapper [FILE1][FILE2]\n\n"
    "Options :\n"
    "-h or --help\t Prints this help message\n"
    "-v or --version\t Prints mapper version\n\n"
    "Mapper accepts two files arguments.\n"
    "The first file contains a reference genome in FASTA format.\n"
    "The second file contains a set of fragments in either FASTA or FASTQ format.\n\n"
    "Mapper parses two files and displays statistics for them which are:\n"
    "names of sequences in the reference file and their lengths,\n"
    "number of sequences in the fragments file, their average length,\n"
    "N50 length, minimal and maximal length.\n\n";
}

void printVersion(){
    std::cout << VERSION << std::endl;
}

int main(int argc, char *argv[]){
    int i;
    
    static struct option long_options[] = {
            {"help", no_argument, &help_flag, 1},
            {"version", no_argument, &version_flag, 1},
            {0, 0, 0, 0}    
    };
    
    while((i = getopt_long(argc, argv, "hv", long_options, nullptr)) != -1){
        
        switch (i)
        {
        case 0:
            break;

        case 'h':
            help_flag = 1;
            break;

        case 'v':
            version_flag = 1;
            break;

        case '?':
            break;    
        
        default:
            abort();
        }
    }

    if(help_flag){
        printHelp();
    } else if(version_flag){
        printVersion();
    }

    return 0;
}
