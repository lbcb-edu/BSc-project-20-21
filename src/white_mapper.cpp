#include <iostream>

#include <getopt.h>

#define VERSION "v1.0"

static int version_req;
static int help_req;

void version_print(){
    std::cout << VERSION << std::endl;
}

void help_print(){
    std::cout << "\nwhite_mapper usage:\n"
    "Two options :\n"
    "-h or --help\t Prints help message\n"
    "-v or --version\t Prints version\n\n"
    "Mapper accepts two files arguments.\n"
    "The first file  FASTA format.\n"
    "The second file  FASTA or FASTQ format.\n\n"
    "Mapper shows these statistics :\n"
    "names of sequences in the reference file and their lengths,\n"
    "number of sequences in the fragments file,\n"
    "their average length,\n"
    "N50 length,\n"
    "minimal and maximal length.\n\n";
}

int main(int argc, char *argv[]){

    int i;
    
    static struct option long_options[] = {
            {"help", no_argument, &help_req, 1},
            {"version", no_argument, &version_req, 1},
            {0, 0, 0, 0}    
    };
    
    while((i = getopt_long(argc, argv, "hv", long_options, nullptr)) != -1){
        
        switch (i)
        {
        case 0:
            break;

        case 'h':
            help_req = true;
            break;

        case 'v':
            version_req = true;
            break;

        case '?':
            break;    
        
        default:
            abort();
        }
    }

    if(version_req){
        version_print();
    } else if(help_req){
        help_print();
    }

    return 0;
}