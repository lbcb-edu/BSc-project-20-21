#include <iostream>
#include <string>
#include <getopt.h>
#include <vector>

#include "bioparser/fasta_parser.hpp"

class Sequence {
public: 
    Sequence(
        const char* name, std::uint32_t name_len,
        const char* data, std::uint32_t data_len) : name_(name, name_len), data_(data, data_len) {
    }

public:
    std::string name_, data_;
};

#define VERSION "v0.1.1"

/* Modificiran primjer https://www.gnu.org/software/libc/manual/html_node/Getopt-Long-Option-Example.html */
/* Pojasnjenje primjera https://www.gnu.org/software/libc/manual/html_node/Getopt-Long-Options.html */
 
static int help_flag = 0;     /* Flag set by ‘--help’.    */
static int version_flag = 0;  /* Flag set by ‘--version’. */

const std::string HELP_MESSAGE = "blonde_mapper usage: \n\n"
                                 "flags: \n"
                                 "-h or --help     prints help message \n"
                                 "-v or --version  prints version      \n\n"
                                 "blonde_mapper takes two filenames as command line arguments: \n"
                                 "The first file will contain a reference genome in FASTA format,\n" 
                                 "while the second file will contain a set of fragments in either\n"
                                 "FASTA or FASTQ format.\n\n";

int main (int argc, char **argv) {
    int c;
    
    while (1) {
        static struct option long_options[] =
            {
            /* These options set a flag. */
            {"help",    no_argument, &help_flag,    1},
            {"version", no_argument, &version_flag, 1},
            {0, 0, 0, 0}
            };
        /* getopt_long stores the option index here. */
        int option_index = 0;

        c = getopt_long (argc, argv, "hv",
                        long_options, &option_index);

        /* Detect the end of the options. */
        if (c == -1)
            break;

        switch (c) {
        case 0:
            break;

        case 'h':
            help_flag = 1;
            break;

        case 'v':
            version_flag = 1;
            break;

        case '?':
            /* getopt_long already printed an error message. */
            break;

        default:
            abort ();
        }
    }

    if (help_flag) {
        std::cout << HELP_MESSAGE;
    } else if (version_flag) {
        std::cout << VERSION << std::endl;
    } if (optind < argc) {
        auto genome_parser = bioparser::Parser<Sequence>::Create<bioparser::FastaParser>(argv[optind]);
        auto genomes = genome_parser->Parse(-1);

        auto fragment_parser = bioparser::Parser<Sequence>::Create<bioparser::FastaParser>(argv[optind + 1]);
        auto fragments = fragment_parser->Parse(-1);

        uint64_t length_sum = 0;
        std::vector<size_t> lengths(fragments.size());
        for (int i = 0; i < int(fragments.size()); i++) {
            lengths[i] = fragments[i]->data_.size();
            length_sum += lengths[i];
        }
        sort(lengths.begin(), lengths.end(), std::greater<size_t>());
        
        uint64_t N50 = -1, tmp_sum = 0;
        for (int i = 0; i < int(fragments.size()); i++) {
            tmp_sum += lengths[i];
            if (tmp_sum * 2 >= length_sum) {
                N50 = lengths[i];
                break;
            }
        }

        std::cerr << "Reference genome sequences:\n";
        for (int i = 0; i < int(genomes.size()); i++) {
            std::cerr << "\t" << genomes[i]->name_ << " , length = " << genomes[i]->data_.size() << '\n';
        }
        std::cerr << "Number of fragments: " << fragments.size() << '\n';
        std::cerr << "Average length: " << length_sum * 1.0 / fragments.size() << '\n';
        std::cerr << "N50 length: " << N50 << '\n';
        std::cerr << "Minimal length: " << lengths.back() << '\n';
        std::cerr << "Maximal length: " << lengths.front() << '\n';
    }

    return 0;
}