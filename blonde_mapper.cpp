#include <iostream>
#include <string>
#include <getopt.h>
#define VERSION "v0.1.0"

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

int main (int argc, char **argv)
{
  int c;

  while (1)
    {
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

      switch (c)
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
          /* getopt_long already printed an error message. */
          break;

        default:
          abort ();
        }
    }

  if (help_flag) {
      std::cout << HELP_MESSAGE;
  }

  if(version_flag) {
      std::cout << VERSION << std::endl;
  }

  /* Print any remaining command line arguments (not options). ()*/
  if (optind < argc)
    {
      std::cout << "non-option ARGV-elements: \n";
      while (optind < argc)
        std::cout << argv[optind++] << std::endl;
    }

  return 0;
}