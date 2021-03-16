/*  Pairwise-alignment-identity-checker (version 0.9.3)
 *  Copyright 2020 by Christoph Mayer
 *
 *  This source file is part of the Pairwise-alignment-identity-checker
 * 
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 *  For any enquiries send an Email to Christoph Mayer
 *  c.mayer.zfmk@uni-bonn.de
 *
 *  When publishing work that is based on the results please cite:
 *  ....
 *  
 */
#include "global-types-and-parameters.h"
#include "tclap/CmdLine.h"
#include <string>
#include "faststring2.h"
#include <climits>

using namespace TCLAP;
using namespace std;

faststring                       global_input_fasta_file;
faststring                       global_data_type_string;

unsigned                         global_min_overlap_len;
unsigned                         global_window_size;
unsigned                         global_window_offset;
unsigned                         global_verbosity;



bool file_exists(const char* fn)
{
  ifstream is;
  is.open(fn);
  if (is.fail())
  {
    return false;
  }
  is.close();
  return true;
}

void good_bye_and_exit(FILE *of, int error)
{
  if (error == 0)
    fprintf(of, "#\n#\n## Finished successfully. ##\n");
  fprintf(of, "\n" PROGNAME " says goodbye.\n\n");
  exit(error);
}

void good_bye_and_exit(ostream &os, int error)
{
  if (error == 0)
    os << "#\n#\n## Finished successfully. ##\n";
  os << "\n" PROGNAME " says goodbye.\n\n";
  exit(error);
}

void init_param()
{
  global_verbosity          = 1;
  global_min_overlap_len    = 15;
  global_window_size        = 30;
  global_window_offset      = 1;
}

void read_and_init_parameters(int argc, char** argv, ostream &logerr)
{
  init_param();

  unsigned tmp_min_overlap_len = 0;
  
  try
  {
    CmdLine cmd("**********************************************************\n"
		"The " PROGNAME  " program uses a sliding window approach to measure the degree of sequence identity/blosum "
		"score in multi sequence alignments. High levels of non-identity/low scores are indications of systematic "
		"problems in the alignment, mostly caused by none homologous sequences. "
		"Accepted input file format: Fasta files of aligned nucleotide or amino acid sequences.\n"
		"The algorithm works as follows: Pairwise identities are computed for all pairs of sequences "
		"in all windows of specified size (parameter window-size) and in steps of an offset (parameter window-offset). "
		PROGNAME " will report the number of sequence pair windows with an identity below "
		"50%, 25% and 10%. From our experience, the number of sequence pair windows below 50% identity "
		"are indicative of problematic alignments, since a 50% identity should normally be found for most "
		"properly aligned pairs of homologous sequences. Exceptions might exist and higher thresholds might be chosen.\n"

		"Recommendations regarding parameter: "
		"In most cases a window size in the range 30-60 bp is recommended, since it results in the detection of partial "
		"sequences that have a small identity."

		"In order to minimise edge effects, we recommend an offset much smaller than the window-size. "
		"In most cases, a window-offset of 1 is a good choice. Of course this will scale up the number of reported pairs, "
		"an effect which should be taken into account when looking at absolute values of pair counts."
		,
		' ', VERSION);

    ValueArg<unsigned> verbosity_Arg("", "verbosity",
				     "The verbosity option controls the amount of information " PROGNAME 
				     " writes to the console while running. 0: Prints only results. "
				     "1: Print also welcome message and essential "
				     "error messages that lead to exiting the program. "
				     "2: Prints progress. Default: 1.",
				     //	"1: report also warnings, 2: report also progress, "
				     // "3: report more detailed progress, >10: debug output. "
				     // "100: debug output."
				     // "Maximum 10000: write all possible diagnostic output."
				     //	" A value of 2 is required if startup parameters should be reported.\n",
				     false, global_verbosity, "unsigned integer");
    cmd.add( verbosity_Arg );
     
    ValueArg<unsigned> min_overlap_Arg("m", "minimum-pairs-in-window",
				       "Minimum number of non gap, non ambiguous nucleotide pairs in two sequences "
				       "in a window such that a meaningful comparison can be made. A good value seems to be "
				       "in the range 10-20. Default: 15, or halfs the window size, whichever is smaller.",
				       false, tmp_min_overlap_len, "unsigned integer");
    cmd.add( min_overlap_Arg );

    ValueArg<unsigned> window_size_Arg("w", "window-size",
				       "Window size in the sliding window analysis. Default: 30.",
				       false, global_window_size, "unsigned integer");
    cmd.add( window_size_Arg );

    ValueArg<unsigned> window_offset_Arg("o", "window-offset",
					 "Window offset in the sliding window analysis. In order to avoid edge effects, we "
					 "recommend a value of 1. Default: 1.",
					 false, global_window_offset, "unsigned integer");
    cmd.add( window_offset_Arg );

    ValueArg<string> data_type_Arg("d", "data-type",
					"Specifies the data type of the input data. For DNA sequences, pairwise identities are "
				        "computed, for protein sequences, blosumn62 scores are computed. "
					"Allowed values: DNA, PROTEIN. Default: DNA.",
					false, "DNA", "string");
    cmd.add( data_type_Arg );

    ValueArg<string> fasta_input_filename_Arg("i", "fasta-input-file-name",
					 "Name of nucleotide input file in fasta format. Sequences have to be aligned.",
					 true, "", "string");
    cmd.add( fasta_input_filename_Arg );

    

//************************************************
    cmd.parse( argc, argv );
//************************************************

    // Assigning parameters to variables:
    global_input_fasta_file                                             = fasta_input_filename_Arg.getValue().c_str();
    global_verbosity                                                    = verbosity_Arg.getValue();
    tmp_min_overlap_len                                                 = min_overlap_Arg.getValue();
    global_window_size                                                  = window_size_Arg.getValue();
    global_window_offset                                                = window_offset_Arg.getValue();
    global_data_type_string                                             = data_type_Arg.getValue().c_str();
  } // END void read_and_init_parameters(int argc, char** argv, ostream &logerr)


  catch (ArgException &e)
  {
    logerr << "ERROR: " << e.error() << " for arg " << e.argId() << "\n";
    good_bye_and_exit(stdout, -1);
  }

  global_data_type_string.ToUpper();
  
  if (global_data_type_string != "DNA" && global_data_type_string != "PROTEIN")
  {
    cerr << "The data type you specified: " << global_data_type_string << " is invalid. "
            "Allowed values are DNA, PROTEIN. Exiting." << endl;
    exit(-3);
  }

  if (global_data_type_string == "PROTEIN")
  {
    cerr << "NOTE: For the data type PROTEIN, the window size is automatically set to " << UINT_MAX << endl;
    global_window_size   = UINT_MAX;
    global_window_offset = 1;
  }
  
  if (tmp_min_overlap_len > 0)
    global_min_overlap_len = tmp_min_overlap_len;
  else
    if (global_window_size/2 < 15)
      global_min_overlap_len = global_window_size/2;
}



//void print_parameters(FILE *of, const char *s)
//{
//       fprintf(of, "%sYYY:                             <not specified>\n", s);
//}


void print_parameters(ostream &os, const char *s)
{
  os << s <<   PROGNAME " was started with the following command line parameters/default values:\n";
  os << s <<   "Input fasta file name:        " <<  global_input_fasta_file  << '\n';
  os << s <<   "Data type:                    " <<  global_data_type_string  << '\n';
  if (global_data_type_string=="DNA")
    os << s <<   "Mode of comparison:           " <<  "Pairwise DNA identity."          << '\n';
  if (global_data_type_string == "PROTEIN")
    os << s <<   "Mode of comparison:           " <<  "Pairwise protein blosum scores." << '\n';
  os << s <<   "Minimum number of nuc pairs:  " <<  global_min_overlap_len   << '\n';
  os << s <<   "Window size:                  " <<  global_window_size       << '\n';
  os << s <<   "Window offset:                " <<  global_window_offset     << '\n';
  os << s <<   "Verbosity:                    " <<  global_verbosity         << "\n\n";
}



void print_calling_command_line(FILE *of, unsigned argc, char ** argv)
{
  faststring cmd;
  
  for (unsigned i=0; i<argc; ++i)
  {
    cmd += argv[i];
    cmd += " ";
  }  

  fprintf(of, "\nThe "PROGNAME " program was called with the following command line:\n");
  fputs(cmd.c_str(), of);
  fputc('\n',of);  
}


void print_calling_command_line(ostream &os, unsigned argc, char ** argv)
{
  faststring cmd;
  
  for (unsigned i=0; i<argc; ++i)
  {
    cmd += argv[i];
    cmd += " ";
  }  

  os << "\nThe " PROGNAME " program was called with the following command line:\n";
  os << cmd.c_str() << '\n';  
}



// void print_output_creation_message(FILE *of, char *s)
// {
//   fprintf(of, "%sOutput created by " PROGNAME ", version " VERSION "\n", s);
//   fprintf(of, "%s\n", s);
//   print_parameters(of, s);
//   fprintf(of, "%s\n", s);
// }

// void print_output_creation_message(ostream &os, char *s)
// {
//   os << s << "Output created by " PROGNAME ", version " VERSION "\n";
//   os << s << '\n';
//   print_parameters(os, s);
//   os << s << '\n';
// }
