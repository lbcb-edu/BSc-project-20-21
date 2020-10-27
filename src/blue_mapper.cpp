#include <iostream>
#include <getopt.h>

const char *mapper_version = "v0.1";

static constexpr option options[] = {
    {"version", no_argument, nullptr, 'v'},
    {"help", no_argument, nullptr, 'h'},
    {nullptr, 0, nullptr, 0}};

void printHelp()
{
    std::cout << "usage: blue_mapper [options] <reference-genome> <fragments>\n"
                 "\n"
                 "  <reference-genome>\treference genome in FASTA format\n"
                 "  <fragments>       \tset of fragments in FASTA/FASTQ format\n"
                 "\n"
                 "  options\n"
                 "    -v, --version\tPrint version number\n"
                 "    -h, --help   \tPrint usage information\n";
}

int main(int argc, char *argv[])
{
    int opt;
    while ((opt = getopt_long(argc, argv, "hv", options, nullptr)) != -1)
    {
        switch (opt)
        {
        case 'v':
            std::cout << mapper_version << std::endl;
            return 0;
        case 'h':
            printHelp();
            return 0;
        default:
            return -1;
        }
    }

    return 0;
}
