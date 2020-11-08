#include <iostream>
#include <getopt.h>

#define VERSION "v0.1";

static option long_options[] = {
    {"help", no_argument, nullptr, 'h'},
    {"version", no_argument, nullptr, 'v'},
    {0, 0, 0, 0}
};

static std::string help = "brown_mapper \n"
                            "-h or --help for help\n"
                            "-v or --version for programs version\n"
                            "Program accepts two files as floating arguments.\n"
                            "The first file needs to contain a reference genome in FASTA format.\n"
                            "The second file needs to contain a set of fragments\n"
                            "in either FASTA or FASTQ format.";

int main(int argc, char* argv[]) {

    int c = getopt_long(argc, argv, "hv", long_options, 0);
    if (c != 1) {
        switch (c){
            case 'h' :
                std::cout << help << std::endl;
                break;
            case 'v' :
                std::cout << VERSION;
                break;
            default:
                break;
        }
    }

    return 0;
}