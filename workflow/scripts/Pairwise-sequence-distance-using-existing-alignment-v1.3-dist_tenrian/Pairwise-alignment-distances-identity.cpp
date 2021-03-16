//    Program: Pairwise-alignment-indentity
//
//    Christoph Mayer
//    Zoological museum Alexander Koenig
//    Adenauerallee 160
//    53113 Bonn
//    Germany
//

#include <iostream>
#include <fstream>
#include <vector>
#include <cstring>
#include <cstdio>
#include <climits>
#include "faststring2.h"
#include "distances_scores.h"
//#include "scoring-matrices/CScoreMatrix.h"
#include "CSequences2.h"
#include "CFile/CFile2_1.h"
#include "global-types-and-parameters.h"

//#include "fast-realloc-vector.h"

using namespace std;

//unsigned global_verbosity;
int      mode;

// Cleans only J characters.
// Reports an error and exits on other problems.
void check_clean_protein_sequence_for_blosum_matrix_compatibility(const faststring &filename, CSequences2 *pseqs)
{
  int i, N = pseqs->GetTaxaNum();
  CSequence_Mol *ps;
  unsigned      offending_char_pos;

  for (i=0; i<N; ++i)
  {
    ps = pseqs->get_seq_by_index(i);
    int res = ps->check_protein_sequence_for_blosum_compatibility(offending_char_pos);
    if (res < 0)
    {
      cerr << "ERROR: While checking the amino acid sequences for their compatibility with the Blosum matrix.\n"
	   << "A non valid character has been found in file: " << filename << endl
	   << "Sequence name: " << ps->getFullName() << endl
	   << "Position:      " << offending_char_pos+1 << endl
	   << "Symbol:        " << (*ps)[offending_char_pos] << endl;
      cerr << "Please resolve this problem before running this program again." << endl;
      exit(1);
    }
    if (res > 0)
    {
      cerr << "# Notice: The symbol J has been replaced by an X " << res << " times in " << endl
	   << "# File:     " << filename << endl
	   << "# Sequence: " << ps->getFullName() << endl << endl;
    }
  }
}


// void pairwise_alignment_scores(CSequences2 *pseqs1, CScoreMatrix *pmat, char outputformat, int special_mode, unsigned min_overlap, unsigned window_size, unsigned window_offset)

void pairwise_alignment_scores(CSequences2 *pseqs1, char outputformat, int special_mode, unsigned min_overlap, unsigned window_size, unsigned window_offset)
{
  // TODO: Extract values from pmat instead of using these values:
  float gap_open                = 3;
  float gap_ext                 = 1;
  float gap_match_score         = 1;
  float factor_front_back_gaps  = 0.5;

  // TODO: Implement use of min_overlap!!!! Not yet used.
  CDists_collection dc(*pseqs1, gap_open, gap_ext, gap_match_score, factor_front_back_gaps, special_mode, min_overlap, window_size, window_offset);

  const std::vector<float> &dists = dc.get_dists();

  cout << "Distances: " << endl;

  unsigned below050=0, below025=0, below010=0, below100=0, below200=0;
  unsigned numDist_good, numDist_undefined, numDist_total;
  
  for (unsigned i=0; i< dists.size(); ++i)
  {
    if (dists[i] != -FLT_MAX)
    {
      if (dists[i] < 0.1)
	++below010;
      if (dists[i] < 0.25)
	++below025;
      if (dists[i] < 0.5)
	++below050;
      if (dists[i] < 1.0)
	++below100;
      if (dists[i] < 2.0)
	++below200;
    }
    //    cout << dists[i] << endl;
  }

  numDist_good       = dc.get_number_of_good_distance_pairs();
  numDist_undefined  = dc.get_number_of_undefined_distance_pairs();
  numDist_total      = numDist_good+numDist_undefined;

  if (special_mode == 128) // DNA no gaps identiy
  {
    cout << "Total number of window pairs:     " << numDist_total     << endl;
    cout << "Number of good window pairs :     " << numDist_good      << endl;
    cout << "Number of undefined window pairs: " << numDist_undefined << endl;

    cout << "Proportion undefined:             " << (float) numDist_undefined/(numDist_total) << endl;
    cout << "Reported are the pairs with an identity below the following values " << endl;

    cout << "Proportion defined below 0.1:     " << (float) below010 /(numDist_good) << endl;
    cout << "Proportion defined below 0.25:    " << (float) below025 /(numDist_good) << endl;
    cout << "Proportion defined below 0.50:    " << (float) below050 /(numDist_good) << endl;
    cout << "Number defined below 0.1:         " << (float) below010  << endl;
    cout << "Number defined below 0.25:        " << (float) below025  << endl;
    cout << "Number defined below 0.50:        " << (float) below050  << endl;

    cout << "Per window pair mean identity     " << dc.get_mean_per_pair_score()  << endl;
    cout << "Per nuc pair mean identity        " << dc.get_mean_per_pos_score()   << endl;
  }
  else if (special_mode == 0) //Protein no gaps Blosum score
  {
    cout << "Total number of sequence pairs:   " << numDist_total     << endl;
    cout << "Number of good sequence pairs :   " << numDist_good      << endl;
    cout << "Number of undefined seq. pairs:   " << numDist_undefined << endl;

    cout << "Proportion undefined:             " << (float) numDist_undefined/(numDist_total) << endl;
    cout << "Reported are the pairs with a Blosum score below the following values " << endl;
    cout << "Proportion defined below 0.1:     " << (float) below010 /(numDist_good) << endl;
    cout << "Proportion defined below 0.25:    " << (float) below025 /(numDist_good) << endl;
    cout << "Proportion defined below 0.50:    " << (float) below050 /(numDist_good) << endl;
    cout << "Number defined below 0.1:         " << (float) below010  << endl;
    cout << "Number defined below 0.25:        " << (float) below025  << endl;
    cout << "Number defined below 0.50:        " << (float) below050  << endl;
    cout << "Number defined below 1.00:        " << (float) below100  << endl;
    cout << "Number defined below 2.00:        " << (float) below200  << endl;
    cout << "Per sequence pair mean score      " << dc.get_mean_per_pair_score()  << endl;
    cout << "Per AA pair mean score            " << dc.get_mean_per_pos_score()   << endl;
  }
}



void Usage_exit()
{
  cerr << "Usage:   Pairwise-alignment-distances-between-two-files (DNA|PROTEIN) fasta-file min-overlap window-size window-offset." << endl;
  cerr << endl;
  exit(-2);
}

void welcome(std::ostream &os)
{
  os << endl << endl;
  os << "      Welcome to " << PROGNAME << ", version " << VERSION << "."
     << endl;
}





int main(int argc, char **argv)
{


  // Read command line options
  read_and_init_parameters(argc, argv, cerr);

  if (global_verbosity >= 1)
  {
    welcome(cerr);
    cerr << endl << endl;
  }
 
  if (global_verbosity >= 1)
    print_parameters(cerr, "");

  if (global_min_overlap_len > global_window_size)
  {
    cerr << "The minimum overlap of sequence pairs in a window is larger than the window size, "
            "excluding all windows for all pairs from the analysis. Please adjust the parameters to "
            "reasonable values. Exiting.\n"; 
    exit(-13);
  }

  // //  cerr << "argc: " << argc << endl;
  // if (argc != 6)
  // {
  //   Usage_exit();
  // }

 
  char       datatype = 0;

  if (global_data_type_string == "PROTEIN")
    datatype = 1;
  
  /*
  faststring datatypestring  = argv[1];
  faststring filename1       = argv[2];
  faststring str_min_overlap = argv[3];
  faststring str_win_size    = argv[4];
  faststring str_win_offset  = argv[5];
  unsigned   min_overlap;
  unsigned   window_size;
  unsigned   window_offset;
  */

  /*
  if (globalstr_min_overlap.isAnUnsigned())
  {
    min_overlap    = str_min_overlap.ToUnsigned();
  }
  else
  {
    cerr << "The third parameter (minimum overlap) must be a positive integer. Exiting.\n";
    exit(-3);
  }

  if (str_win_size.isAnUnsigned())
  {
    window_size   = str_win_size.ToUnsigned();
    if (window_size == 0)
    {
      window_size = UINT_MAX/2;
    }
  }
  else
  {
    cerr << "The forth parameter (window-size) must be a positive integer. Exiting.\n";
    exit(-3);
  }

  if (str_win_offset.isAnUnsigned())
  {
    window_offset = str_win_offset.ToUnsigned();
    if (window_offset == 0)
    {
      window_offset = UINT_MAX/2;
    }
  }
  else
  {
    cerr << "The fifth parameter (window-offset) must be a positive integer. Exiting.\n";
    exit(-3);
  }
  */

  /*  
  if (global_da == "dna")
    datatype = 0;
  else if (datatypestring == "protein")
    datatype = 1;
  else
  {
    cerr << "ERROR: Unknown data type: " << datatypestring << endl;
    cerr << "Allowed data types are: \'dna\', \'protein\'." << endl << endl;
    Usage_exit();
  }
  */

  CSequences2 *seqs1;
  //  CScoreMatrix *pmat;

  if (datatype == 0)
  {
    seqs1 = new CSequences2(CSequence_Mol::dna);
    //    cerr << "Reading score matrix: standard-dna.mat" << endl;
    //    pmat = new CScoreMatrix("standard-dna.mat");
  }
  else
  {
    seqs1 = new CSequences2(CSequence_Mol::protein);
    //    cerr << "Reading score matrix: blosum50.mat" << endl;
    //    pmat = new CScoreMatrix("blosum50.mat");
  }

  CFile file1;
  CSequence_Mol::processing_flag pflag = CSequence_Mol::processing_flag(CSequence_Mol::convert_toupper);
  if (global_verbosity > 100)
    cerr << "Processing flag: " << pflag << endl;

  if (global_verbosity >= 2)
    cerr << "Reading fasta file: " << global_input_fasta_file << endl;

  file1.ffopen(global_input_fasta_file.c_str());
  if (file1.fail())
  {
    cerr << "Could not open specified file: " << global_input_fasta_file << endl;
    exit(-1);
  }
  if (global_verbosity > 100)
    seqs1->read_from_Fasta_File(file1, pflag, 0, -1, true);
  else
    seqs1->read_from_Fasta_File(file1, pflag, 0, -1, false);
  file1.ffclose();

  if (global_verbosity >= 2)
    cerr << "Finished reading sequences." << endl;

  if (datatype == 0) // No sequence alphabet check implemented yet.
  {
  }
  else
  {
    //    cerr << "Checking sequences that they are compatible with the blosum matrix.\n";
    check_clean_protein_sequence_for_blosum_matrix_compatibility(global_input_fasta_file, seqs1);
    //    cerr << "Finished checking sequences.\n";
  }

  //  cerr << "Printing sequences: " << endl;
  //  seqs1->ExportSequences_no_fill_in_ext_phylip_range(cerr, 0u, 200000000u);
  //  seqs2->ExportSequences_no_fill_in_ext_phylip_range(cerr, 0u, 200000000u);  
  //  cerr << "Finished printing sequences." << endl;

           char outputformat = '\0';
  unsigned char special_mode;
  
  // TODO:
  if (datatype == 0)
    special_mode = 128;
  else
    special_mode = 0;
  
  //  pairwise_alignment_scores(seqs1, pmat, outputformat, special_mode, global_min_overlap_len, global_window_size, global_window_offset);
  pairwise_alignment_scores(seqs1, outputformat, special_mode, global_min_overlap_len, global_window_size, global_window_offset);

  if (global_verbosity >= 1)
    cerr << "Finished without errors." << endl;
}


