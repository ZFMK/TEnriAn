#ifndef DISTANCES_SCORES_H
#define DISTANCES_SCORES_H

#include "faststring2.h"
#include "aa-matrices.h"
#include <cfloat>
#include <climits>

#include <iostream>
#include "CSequences2.h"
#include "CSequence_Mol2_1.h"
#include <vector>
#include "statistic_functions.h"
#include "basic-DNA-RNA-AA-routines.h"

// File in Klassen

/*
inline int get_index_for_index_pair(int i, int j)
{
  // In all_AA_distances, all_AA_distances_gaps we push the distances in the order:
  //    (1,0), (2,0), (2,1), (3,0), (3,1), (3,2), (4,0), ...
  // They 0,    1,     2,     3,  
  return i*(i-1)/2+j;
}
*/


// Fully implemented
inline void aaScore_BLOSUM62(const char* s1,
			     const char* s2,
			     unsigned N,
			     int      &blosum62_sum,
			     unsigned &len,
			     unsigned &unknown
			     )
{
  unsigned i;
  
  int      blosum62_sum_local  = 0;
  unsigned len_local           = 0;
  unsigned unknown_local       = 0;
  
  char c1, c2;
  int  res;
  
  for (i=0; i<N; ++i)
  {
    c1 = s1[i];
    c2 = s2[i];

    if ( !(c1 == '-' || c2 == '-') )
    {
      res = get_BLOSUM62_score(c1, c2);
      //      std::cout << "c1 c2 dist " << c1 << " " << c2 << " " << res << '\n';
 
      if (res == 255)
      {
	std::cerr << "Found unexpected symbol(s). No score defined for these symbols. " << c1 << " " << c2 << '\n';
	++unknown_local;
      }
      else
      {
	//	std::cerr << c1 << " " << c2 << " " << (int)res << '\n';
	blosum62_sum_local += res;
	++len_local;
      }
    }
  }
  blosum62_sum  = blosum62_sum_local;
  len           = len_local;
  unknown       = unknown_local;
}

inline void aaScore_BLOSUM62(const char* s1,
			     const char* s2,
			     unsigned range_beg,
			     unsigned range_end,
			     int      &blosum62_sum,
			     unsigned &len,
			     unsigned &unknown
			     )
{
  unsigned i;

  int      blosum62_sum_local  = 0;
  unsigned len_local           = 0;
  unsigned unknown_local       = 0;

  char c1, c2;
  int  res;

  for (i=range_beg; i<range_end; ++i)
  {
    c1 = s1[i];
    c2 = s2[i];

    if ( !(c1 == '-' || c2 == '-') )
    {
      res = get_BLOSUM62_score(c1, c2);

      //      std::cout << "c1 c2 dist " << c1 << " " << c2 << " " << res << '\n';
 
      if (res == 255)
      {
	std::cerr << "Found unexpected symbol(s). No score defined for these symbols. "
		  << c1 << " " << c2 << '\n';
	++unknown_local;
      }
      else
      {
	//	std::cerr << c1 << " " << c2 << " " << (int)res << '\n';
	blosum62_sum_local += res;
	++len_local;
      }
    }
  }
  blosum62_sum  = blosum62_sum_local;
  len           = len_local;
  unknown       = unknown_local;
}


// Fully implemented
inline void aaScore_BLOSUM62(faststring &s1,
			     faststring &s2, 
			     int        &blosum62_sum,
			     unsigned   &len,
			     unsigned   &unknown
			     )
{
  if (s1.length() != s2.length())
  {
    blosum62_sum  = 0;  // These values are indicative of this error or empty sequences.
    len           = 0;
    unknown       = 0;
    return;
  }
  aaScore_BLOSUM62(s1.c_str(), s2.c_str(), s1.length(), blosum62_sum, len, unknown);
}

// Fully implemented
inline void aaScore_BLOSUM62(faststring &s1,
			     faststring &s2,
			     unsigned   range_beg,
			     unsigned   range_end,
			     int        &blosum62_sum,
			     unsigned   &len,
			     unsigned   &unknown
			     )
{
  if (s1.length() < range_beg || s1.length() < range_end ||
      s2.length() < range_beg || s2.length() < range_end    )
  {
    blosum62_sum  = 0;  // These values are indicative of this error or empty sequences.
    len           = 0;
    unknown       = 0;
    return;
  }
  aaScore_BLOSUM62(s1.c_str(), s2.c_str(), range_beg, range_end, blosum62_sum, len, unknown);
}

// not verified !!!!!!!!!!!!

// factor_front_back_gaps should be >0. Should not be 0.
// Scoring gaps:
// AA--AA:
// AA--AA: 0 gap-openings, 2 gap-matches
//
// ----AA--:
// --AAAA--: 1 gap-opeing, 4 gap-matches
//
// AA----AA:
// AAAA--AA: 1 gap-opening, 4 gap-matches

inline void aaScore_BLOSUM62_gaps(const char* s1,
				   const char* s2,
				   unsigned N,
				   float    gap_open,
				   float    gap_ext,
		    		   float    gap_match_score,
				   float    factor_front_back_gaps,
				   float    &blosum62_sum,
				   unsigned &len,
				   unsigned &unknown
				   )
{
  int i;

  float    blosum62_sum_local = 0;
  unsigned len_local           = 0;
  unsigned unknown_local       = 0;

  int      gap_s1 = 0;
  int      gap_s2 = 0;

  char c1, c2;
  int  res;

  int start, end;

  bool gap_at_start_s1 = false;
  bool gap_at_start_s2 = false;

  bool in_gap;

  end   = N;

/*   // Determine start for comparison: */
/*   for (i=0; i< N; ++i) */
/*   { */
/*     if ( s1[i] != '-' || s2[i] != '-' ) */
/*       break; */
/*     else */
/*     { */
/*       blosum62_sum_local += gap_gap_match_score * factor_front_back_gaps; */
/*       ++len_local; */
/*     } */

/*   } */
/*   start = i; */
 
/*   for (i=N-1; i>=start; --i) */
/*   { */
/*     if ( s1[i] != '-' || s2[i] != '-' ) */
/*       break; */
/*     else */
/*     { */
/*       blosum62_sum_local += gap_gap_match_score * factor_front_back_gaps; */
/*       ++len_local; */
/*     } */
/*   } */
/*   end = i+1; // i is the position with */

  // In the range (start .. end-1) we start and end with a position that does
  // not have a gap in both sequences.

  start = 0;
  end   = N;

  if (start >= end)
  {
    blosum62_sum  = 0;
    len           = 0;
    unknown       = 0;
    return;
  }

  // Only one of the two ifs can be true:
  if (s1[start] == '-')
    gap_at_start_s1 = true;
  else if (s2[start] == '-')
    gap_at_start_s2 = true;

  for (i=start; i<end; ++i)
  {
    c1 = s1[i];
    c2 = s2[i];

    if ( c1 == '-' || c2 == '-' ) // => either gap-gap or gap-nongap or nongap-gap
    {
      if (c1 != c2)
      {
	++len_local; // We count this position
	if (c1 == '-')
	{
	  // This could be the end of a gap in s2:
	  if (gap_s2 > 0)
	  {
	    if (gap_at_start_s2)
	    {
	      blosum62_sum_local += (gap_open + gap_s2*gap_ext)*factor_front_back_gaps;
	      gap_s2 = 0;
	      gap_at_start_s2 = false;
	    }
	    else
	    {
	      blosum62_sum_local += (gap_open + gap_s2*gap_ext);
	      gap_s2 = 0;
	    }
	  }
	  if ( is_aa_blosum_code_index(c2) )
	  {
	    ++gap_s1;
	  }
	  else
	  {
	    ++unknown_local;
	    --len_local;
	  }
	}
	else // must be (c2 == '-')
	{
	  // This could be the end of a gap in s1:
	  if (gap_s1 > 0)
	  {
	    if (gap_at_start_s1)
	    {
	      blosum62_sum_local += (gap_open + gap_s1*gap_ext)*factor_front_back_gaps;
	      gap_s1 = 0;
	      gap_at_start_s1 = false;
	    }
	    else
	    {
	      blosum62_sum_local += (gap_open + gap_s1*gap_ext);
	      gap_s1 = 0;
	    }
	  }
	  if ( is_aa_blosum_code_index(c1) )
	  {
	    ++gap_s2;
	  }
	 else
	 {
	   ++unknown_local;
	   --len_local;
	 }
	}
      } // END if (c1 != c2)
      else // We have a gap in both sequences.
      {
	
      }
    }
    else // ! ( c1 == '-' || c2 == '-' ) => nongap-nongap
    {
      if (gap_s1 > 0)
      {
	if (gap_at_start_s1)
	{
	  blosum62_sum_local += (gap_open + gap_s1*gap_ext)*factor_front_back_gaps;
	  gap_s1 = 0;
	  gap_at_start_s1 = false;
	}
	else
	{
	  blosum62_sum_local += (gap_open + gap_s1*gap_ext);
	  gap_s1 = 0;
	}
      }
      else if (gap_s2 > 0) // DOUBLECHECK
      {
	if (gap_at_start_s2)
	{
	  blosum62_sum_local += (gap_open + gap_s2*gap_ext)*factor_front_back_gaps;
	  gap_s2 = 0;
	  gap_at_start_s2 = false;
	}
	else
	{
	  blosum62_sum_local += (gap_open + gap_s2*gap_ext);
	  gap_s2 = 0;
	}
      }

      ++len_local;
      res = get_BLOSUM62_score(c1, c2);
      if (res == 255)
      {
	std::cerr << "Found unexpected symbol(s). No score defined for these symbols. " << c1 << " " << c2 << '\n';
	++unknown_local;
	--len_local;
      }
      else
      {
	//	std::cerr << c1 << " " << c2 << " " << (int)res << '\n';
	blosum62_sum_local += res;
      }
    }
  } // END for (i=0; i<N; ++i)

  // Treat end gaps
  if (gap_s1 > 0)
  {
    blosum62_sum_local += (gap_open + gap_s1*gap_ext)*factor_front_back_gaps;
    gap_s1 = 0;
  }
  if (gap_s2 > 0)
  {
    blosum62_sum_local += (gap_open + gap_s2*gap_ext)*factor_front_back_gaps;
    gap_s2 = 0;
  }

  blosum62_sum  = blosum62_sum_local;
  len           = len_local;
  unknown       = unknown_local;
}

// Fully implemented
inline void all_AA_distances(CSequences2 &seqs, std::vector<float> *dists,
			     unsigned min_overlap, unsigned &sum_length,
			     float &sum_scores_of_pairs, float &sum_scores_of_positions,
			     unsigned &num_good_pairs, unsigned &num_undefined_pairs)
{
  unsigned num_good_pairs_local      = 0;
  unsigned num_undefined_pairs_local = 0;
  
  sum_length             = 0;
  sum_scores_of_pairs     = 0;
  sum_scores_of_positions = 0;
  
  unsigned N      = seqs.GetTaxaNum();
  unsigned posNum = seqs.GetPosNum();

  // CSequences2 is a class for aligned sequences, but we check whether all sequences
  // have equal length before we proceed.
  //  if (!seqs.equal_length_of_all_sequences() )
  //    return false;

  if (posNum == 0)
  {
    std::cerr << "Error: sequences passed to the CDists_collection object do not seem to have equal "
                 "lengths in all_AA_distances. It could also be that they have been read with "
                 "the remove gaps option.\n";
    std::exit(-97);
  }
  
  unsigned  i,j;
  int       tmp;
  unsigned  l;
  unsigned  u;

  float d;

  if (dists != NULL)
  {
    dists->reserve( (N-1)*N/2 + 1 );

    for (i=0; i<N; ++i)
    {
      for (j=0; j<i; ++j)
      {
	aaScore_BLOSUM62(seqs.get_Seq_Data(i),
			     seqs.get_Seq_Data(j),
			     posNum,
			     tmp, l, u);
	if (l >= min_overlap)
	{
	  ++num_good_pairs_local;
	  d = (float) tmp/l;
	  sum_scores_of_pairs      += d;
	  sum_scores_of_positions  += tmp;
	  sum_length += l;
	}
	else
	{
	  ++num_undefined_pairs_local;
	  d = -FLT_MAX;
	}
	//	std::cout << "(" << i << "," << j << "): " << d << '\n'; 
	//	std::cout << "Adding(1): " << d << '\n';
	//	std::cout << seqs.get_Seq_Data(i) << '\n';
	//	std::cout << seqs.get_Seq_Data(j) << '\n';
	//	std::cout << "with tmp: " << tmp << " l: " << l << '\n';
	//	std::cout << "New sum: " << sum_scores_of_pairs << '\n';
	//	std::cout << sum_scores_of_pairs << '\n';	
	dists->push_back(d);
      }
    }
  }
  else
  {
    for (i=0; i<N; ++i)
    {
      for (j=0; j<i; ++j)
      {
	aaScore_BLOSUM62(seqs.get_Seq_Data(i),
			     seqs.get_Seq_Data(j),
			     posNum,
			     tmp, l, u);
	if (l>= min_overlap)
	{
	  ++num_good_pairs_local;
	  d = (float) tmp/l;
	  sum_scores_of_pairs      += d;
	  sum_scores_of_positions  += tmp;
	  sum_length += l;
	}
	else
	{
	  ++num_undefined_pairs_local;
	}
	//	std::cout << "(" << i << "," << j << "): " << d << '\n'; 
	//	std::cout << "Adding(2): " << (float)tmp/l << '\n';
	//	std::cout << sum_length << '\n';

	//	std::cout << sum_length << '\n';
      }
    }
  }
  num_good_pairs      = num_good_pairs_local;
  num_undefined_pairs = num_undefined_pairs_local;
}

// Not verified ???
inline void all_AA_distances_gaps(CSequences2 &seqs, std::vector<float> *dists,
				  float gap_open, float gap_ext, float gap_match_score,
				  float factor_front_back_gaps, unsigned min_overlap, unsigned &sum_length,
				  float &sum_scores_of_pairs, float &sum_scores_of_positions,
				  unsigned &num_good_pairs, unsigned &num_undefined_pairs)
{
  unsigned num_good_pairs_local      = 0;
  unsigned num_undefined_pairs_local = 0;

  sum_length          = 0;
  sum_scores_of_pairs = 0;
  num_good_pairs      = 0;
  sum_scores_of_positions = 0;

  unsigned N      = seqs.GetTaxaNum();
  unsigned posNum = seqs.GetPosNum();

  if (posNum == 0)
  {
    std::cerr << "Error: sequences passed to the CDists_collection object do not seem to have equal "
                 "lengths in all_AA_distances_gaps. It could also be that they have been read with "
                 "the remove gaps option.\n";

    std::exit(-97);
  }
  
  // CSequences2 is a class for aligned sequences, but we check whether all sequences
  // have equal length before we proceed.
  //  if (!seqs.equal_length_of_all_sequences() )
  //    return false;

  unsigned  i,j;
  float     tmp;
  unsigned  l;
  unsigned  u;

  float d;

  if (dists != NULL)
  {
    dists->reserve( (N-1)*N/2 + 1 );

    for (i=0; i<N; ++i)
    {
      for (j=0; j<i; ++j)
      {
	aaScore_BLOSUM62_gaps(seqs.get_Seq_Data(i),
			      seqs.get_Seq_Data(j),
			      posNum,
			      gap_open, gap_ext, gap_match_score, factor_front_back_gaps,
			      tmp, l, u);
	if (l>= min_overlap)
	{
	  ++num_good_pairs;
	  d = (float) tmp/l;
	  sum_scores_of_pairs  += d;
	  sum_scores_of_positions  += tmp;
	  sum_length += l;
	  //	std::cout << "Adding: " << d << '\n';
	}
	else
	{
	  ++num_undefined_pairs;
	  d = -FLT_MAX;
	}
	//	std::cout << "Adding: " << d << '\n';
	dists->push_back(d);
      }
    }
  }
  else // if (dists != NULL)
  {
    for (i=0; i<N; ++i)
    {
      for (j=0; j<i; ++j)
      {
	aaScore_BLOSUM62_gaps(seqs.get_Seq_Data(i),
			      seqs.get_Seq_Data(j),
			      posNum,
			      gap_open, gap_ext, gap_match_score, factor_front_back_gaps,
			      tmp, l, u);

	if (l>= min_overlap)
	{
	  ++num_good_pairs;
	  d = (float) tmp/l;
	  sum_scores_of_pairs  += d;
	  sum_scores_of_positions  += tmp;
	  sum_length += l;
	}
	else
	{
	  ++num_undefined_pairs;
	}
	//	std::cout << "Adding: " << d << '\n';
      }
    }
  }
  num_good_pairs      = num_good_pairs_local;
  num_undefined_pairs = num_undefined_pairs_local;
}

// Count differences between two nucleotide (DNA) sequences:
// AA: identity is increased by 1. len is increased by 1.*
// AC: len is increased by 1.*
// RA: ambig is increased by 1. *
// RR: ambis is increased by 1.
// -A: gap is increased by 1.
// -R: gap is increased by 1. ambig is not increased by 1!
// --: gap is increased by 1.*
// OA: unknown is increased by 1.
// OR: unknown is increased. Ambig is not increases.
// O-: unknown is increased by 1. gap in not increased by one.
//
// Important: Does not check whether start and start+N < sequence-length.
inline void nucStrictIdentity(const char *s1,
			      const char *s2,
			      unsigned N, 
			      int      &identity, 
			      unsigned &len,
			      unsigned &ambig,
			      unsigned &gap,
			      unsigned &unknown,
			      unsigned start=0 )
{
  unsigned i;

  unsigned identity_sum_local = 0;
  unsigned len_dna_local      = 0;
  unsigned ambig_local        = 0;
  unsigned gap_local          = 0;
  unsigned gap_or_ambig_local = 0;
  unsigned unknown_local      = 0;

  char c1, c2;

  unsigned last_N = start+N;
  
  for (i=start; i < last_N; ++i)
  {
    c1 = toupper(s1[i]);
    c2 = toupper(s2[i]);

    /*    
    if (c1 == 0 || c2 == 0)
    {
      std::cerr << "Assertion failed here. if (c1 == 0 || c2 == 0)\n";
      exit(-1);
      //	      std::
    }
    */     
    
    if ( is_DNA_base(c1) && is_DNA_base(c2) ) // This is the most probale case.
    {
      ++len_dna_local;
      if (c1 == c2)
	++identity_sum_local;
    }
    else if (is_DNA_or_DNA_ambig_or_GAP(c1) && is_DNA_or_DNA_ambig_or_GAP(c2)) // This is the next most probable case.
    {
      // c1 and c2 are valid symbols. They are not both DNA symbols, so one is a gap or an ambig symbol.
      ++gap_or_ambig_local;
      if (c1 == '-' || c2 == '-')
      {
	++gap_local;
	//	--len_dna_local;  // VERIFY
      }
      else
      {
	++ambig_local;
      }
    }
    else // c1 or c2 is unknown
    {
      ++unknown_local;
    }
  }
  identity      = identity_sum_local;
  len           = len_dna_local;
  unknown       = unknown_local;
}

// Code taken from distdna.c 

enum nt_changes{
    AA = 0,
    AC,
    AG,
    AT,
    CA,
    CC,
    CG,
    CT,
    GA,
    GC,
    GG,
    GT,
    TA,
    TC,
    TG,
    TT
};



/**
 * pairwise frequencies
 *
 * compute all pairwise nucleotide combination frequencies among two sequences.
 *
 * @param char *s1	The first sequence
 * @param char *s2	The second sequence
 * @param int slen	Sequences length
 * @param double f[16]	An array to store all 4x4 combination possibilities
 *
 */

inline void pairwise_frequencies(char *s1, char *s2, int slen, double *f)
{
    long int i, validpos;

    /* reset frequency counts */
    for (i = 0; i < 16; i++) f[i] = 0.;

    /* count pairwise frequencies and valid positions */
    validpos = slen;
    for (i = 0; i < slen; i++) {
        switch (s1[i]) {
	case 'A':
	    switch(s2[i]) {
	        case 'A': f[AA] += 1.0; break;
		case 'G': f[AG] += 1.0; break;
		case 'C': f[AC] += 1.0; break;
		case 'T': f[AT] += 1.0; break;
		default: validpos--; break;
	    }
	    break;
	case 'C':
	    switch(s2[i]) {
	        case 'A': f[CA] += 1.0; break;
		case 'G': f[CG] += 1.0; break;
		case 'C': f[CC] += 1.0; break;
		case 'T': f[CT] += 1.0; break;
		default: validpos--; break;
	    }
	    break;
	case 'G':
	    switch(s2[i]) {
	        case 'A': f[GA] += 1.0; break;
		case 'G': f[GG] += 1.0; break;
		case 'C': f[GC] += 1.0; break;
		case 'T': f[GT] += 1.0; break;
		default: validpos--; break;
	    }
	    break;
	case 'T':
	    switch(s2[i]) {
	        case 'A': f[TA] += 1.0; break;
		case 'G': f[TG] += 1.0; break;
		case 'C': f[TC] += 1.0; break;
		case 'T': f[TT] += 1.0; break;
		default: validpos--; break;
	    }
	    break;
	default: validpos--; break;
	}
    }
    
    /* convert counts to frequencies */
    if (validpos != 0) {
        for (i = 0; i < 16; i++) f[i] /= (double) validpos;
    }
    /* else all f[i] == 0 */
}

/* double distance_p(char *s1, char *s2, char *s1_end) */
/* { */
/*   while (s1 < s1_end) */
/*   { */
/*     ++s1; */
/*     ++s2; */


/*   } */
/* } */


// TODO: The parameters gapwt and abmig do not make sense at the moment. The only reasonable values are 0 for both parameters.
//       But if the only reasonable values are 0, the parameters can be removed.

// Returns -FLT_MAX if the effective_len is to short such that we have a division by 0 or meaningless small number. (effective_len < 0.5).
inline double dist_p_distdna(const char *s1, const char *s2, int slen, double gapwt=0, int ambig=0)
{
    double matches;
    int    gaps, k;
    double dist;

    /* compute unambiguous sum of matches and gaps */
    gaps = 0;
    matches = 0.;
    for (k = 0; k < slen; k++) {
	if ((is_gap(s1[k])) || (is_gap(s2[k])))
	    gaps++;
	else if (is_DNA_ambig(s1[k]) || is_DNA_ambig(s2[k])) /*** Added by me ***/
	  ++gaps;
	//	else if (ambig)
	//	    matches += nt_ambig_match(s1[k], s2[k]);
	else if (s1[k] == s2[k])
	    matches += 1.0;
    }

    double effective_len = ((double) (slen - gaps) + ((double)gaps * gapwt));

    /* avoid divide by zero: return maximum distance */
    if (effective_len < 0.5)
    {
/*       std::cerr << "WARNING: effective_len<0.5: " << effective_len << '\n'; */
      return -FLT_MAX;
    }

    /* compute uncorrected distance */
    dist = 1.0 - (matches / effective_len) ;

/*     if (dist < 0) */
/*     { */
/*       std::cerr << "WARNING: Computed distance is smaller than 0: " << dist << '\n'; */
/*     } */
    return dist;
}


// TODO: The parameters gapwt and abmig do not make sense at the moment. The only reasonable values are 0 for both parameters.
//       But if the only reasonable values are 0, the parameters can be removed.

// Returns -FLT_MAX if the effective_len is to short such that we have a division by 0 or meaningless small number. (effective_len < 0.5).
inline double dist_p_distdna(const char *s1, const char *s2, int slen, int &len_no_gaps, double gapwt=0, int ambig=0)
{
    double matches;
    int    gaps, k;
    double dist;

    /* compute unambiguous sum of matches and gaps */
    gaps = 0;
    matches = 0.;
    for (k = 0; k < slen; k++) {
	if ((is_gap(s1[k])) || (is_gap(s2[k])))
	    gaps++;
	else if (is_DNA_ambig(s1[k]) || is_DNA_ambig(s2[k])) /*** Added by me ***/
	  ++gaps;
	//	else if (ambig)
	//	    matches += nt_ambig_match(s1[k], s2[k]);
	else if (s1[k] == s2[k])
	    matches += 1.0;
    }

    double effective_len = ((double) (slen - gaps) + ((double)gaps * gapwt));

    len_no_gaps = slen-gaps;

    /* avoid divide by zero: return maximum distance */
    if (effective_len < 0.5)
    {
      std::cerr << "WARNING: effective_len<0.5: " << effective_len << '\n';
      return -FLT_MAX;
    }

    /* compute uncorrected distance */
    dist = 1.0 - (matches / effective_len) ;

    if (dist < 0)
    {
      std::cerr << "WARNING: Computed distance is smaller than 0: " << dist << '\n';
    }
    return dist;
}



// TODO: The parameters gapwt and abmig do not make sense at the moment. The only reasonable values are 0 for both parameters.
//       But if the only reasonable values are 0, the parameters can be removed.

inline double dist_p_distPROTEIN(const char *s1, const char *s2, int slen, double gapwt=0, int ambig=0)
{
    double matches;
    int    gaps, k;
    double dist;

    /* compute unambiguous sum of matches and gaps */
    gaps = 0;
    matches = 0.;
    for (k = 0; k < slen; k++) {
	if ((is_gap(s1[k])) || (is_gap(s2[k])))
	    gaps++;
	else if (is_aa_ambig(s1[k]) || is_aa_ambig(s2[k])) /*** Added by me ***/
	  ++gaps;
	//	else if (ambig)
	//	    matches += nt_ambig_match(s1[k], s2[k]);
	else if (s1[k] == s2[k])
	    matches += 1.0;
    }

    double effective_len = ((double) (slen - gaps) + ((double)gaps * gapwt));

    /* avoid divide by zero: return maximum distance */
    if (effective_len < 0.5)
      return -FLT_MAX;

    /* compute uncorrected distance */
    dist = 1.0 - (matches / effective_len) ;

    return dist;
}

// Return the corrected JC distance.
// In case of too large a p-distance, such that a correction is impossible, it returns the negarive p distance.
// The calling program has to check for negative distances and deal with them.
// Negative distances occur here if the p-distance is above or equal to 75% or if the effective_len < 0.5.
// A distance above 0.75 should not occur in real data sets. If such a distance is found this hints at a problem,
// e.g. a sequence overlap which is too short.

// TODO: The parameters gapwt and abmig do not make sense at the moment. The only reasonable values are 0 for both parameters.
//       But if the only reasonable values are 0, the parameters can be removed.

inline double dist_jukes_cantor_dnadist(const char *s1, const char *s2, int slen, double gapwt=0, int ambig=0)
{
    double matches;
    int gaps, k;
    double dist;

    /* compute unambiguoug sum of matches and gaps */
    matches = gaps = 0;
    for (k = 0; k < slen; k++) {
	if ((is_gap(s1[k])) || (is_gap(s2[k])))
	    gaps++;
	else if (is_DNA_ambig(s1[k]) || is_DNA_ambig(s2[k])) /*** Added by me ***/
	  ++gaps;
	//	else if (ambig)
	//	    matches += nt_ambig_match(s1[k], s2[k]);
	else if (s1[k] == s2[k])
	    matches += 1.0;
    }

   double effective_len = ((double) (slen - gaps) + ((double)gaps * gapwt));

    /* avoid divide by zero: return maximum distance */
    if (effective_len < 0.5)
      return -FLT_MAX;

    /* compute uncorrected distance */
    /*
     * p = proportion of different nt
     * n = number of nt examined
     * k = number of nt pairs that differ
     * p ~= k/n = 1 - (matches/npos)
     */
    dist = 1.0 - (matches / effective_len) ;

    /* compute Jukes-Cantor distance */
    /*
     * d = -b log_e[1 - (p/b)]
     * assuming equal frequencies for all nt, f(nt) = 1/4 and b=3/4
     */
    if (dist >= 0.75)
      return -dist;  // This distance cannot be corrected.
    else
      return  -(3. / 4.) * log(1. - (dist / (3. / 4.)));
}


// Return the corrected JC distance.
// In case of too large a p-distance, such that a correction is impossible, it returns the negarive p distance.
// The calling program has to check for negative distances and deal with them.
// Negative distances occur here if the p-distance is above or equal to 75% or if the effective_len < 0.5.
// A distance above 0.75 should not occur in real data sets. If such a distance is found this hints at a problem,
// e.g. a sequence overlap which is too short.

// TODO: The parameters gapwt and abmig do not make sense at the moment. The only reasonable values are 0 for both parameters.
//       But if the only reasonable values are 0, the parameters can be removed.

inline double dist_jukes_cantor_dnadist(const char *s1, const char *s2, int slen, int &len_no_gaps, double gapwt=0, int ambig=0)
{
    double matches;
    int gaps, k;
    double dist;

    /* compute unambiguoug sum of matches and gaps */
    matches = gaps = 0;
    for (k = 0; k < slen; k++) {
	if ((is_gap(s1[k])) || (is_gap(s2[k])))
	    gaps++;
	else if (is_DNA_ambig(s1[k]) || is_DNA_ambig(s2[k])) /*** Added by me ***/
	  ++gaps;
	//	else if (ambig)
	//	    matches += nt_ambig_match(s1[k], s2[k]);
	else if (s1[k] == s2[k])
	    matches += 1.0;
    }

   double effective_len = ((double) (slen - gaps) + ((double)gaps * gapwt));

   len_no_gaps = slen - gaps;

    /* avoid divide by zero: return maximum distance */
    if (effective_len < 0.5)
      return -FLT_MAX;

    /* compute uncorrected distance */
    /*
     * p = proportion of different nt
     * n = number of nt examined
     * k = number of nt pairs that differ
     * p ~= k/n = 1 - (matches/npos)
     */
    dist = 1.0 - (matches / effective_len) ;

    /* compute Jukes-Cantor distance */
    /*
     * d = -b log_e[1 - (p/b)]
     * assuming equal frequencies for all nt, f(nt) = 1/4 and b=3/4
     */
    if (dist >= 0.75)
      return -dist;  // This distance cannot be corrected.
    else
      return  -(3. / 4.) * log(1. - (dist / (3. / 4.)));
}


/// TODO: HANDLE ambigs!!!!

// Return the corrected JC distance.
// In case of too large a p-distance >= 0.75, such that a correction is impossible, it returns the negarive p distance.
// If the effective_len is 0, (<0.5 to allow small floating point inacuracies), a distance of -FLT_MAX is returned.
// => A negative distance indicates a problem, the value reveals the reason.
// The calling program has to check for negative distances and deal with them.
// A distance >= 75% should not occur in real data sets. If such a distance is found this hints at a problem,
// e.g. a sequence overlap which is too short.
//
// TODO: The parameters gapwt and abmig do not make sense at the moment. The only reasonable values are 0 for both parameters.
//       But if the only reasonable values are 0, the parameters can be removed. gapwt (gapweight) is used).

inline double dist_jukes_cantor_distPROTEIN(const char *s1, const char *s2, int slen, double gapwt=0, int ambig=0)
{
    double matches;
    int gaps, k;
    double dist;

    /* compute unambiguoug sum of matches and gaps */
    matches = gaps = 0;
    for (k = 0; k < slen; k++) {
	if ((is_gap(s1[k])) || (is_gap(s2[k])))
	    gaps++;	
	/** I added this. Ambigs should not be treated as mismatches. **/
	else if (is_aa_ambig(s1[k]) || is_aa_ambig(s2[k]))
	  ++gaps;
	//	else if (ambig)
	//	    matches += nt_ambig_match(s1[k], s2[k]);
	else if (s1[k] == s2[k])
	    matches += 1.0;
    }

    /* avoid divide by zero */
    double effective_len = ((double) (slen - gaps) + ((double)gaps * gapwt));

    if (effective_len < 0.5)
      return -FLT_MAX;

    /* compute uncorrected distance */
    /*
     * p = proportion of different nt
     * n = number of nt examined
     * k = number of nt pairs that differ
     * p ~= k/n = 1 - (matches/npos)
     */
    dist = 1.0 - (matches / effective_len );

    /* compute Jukes-Cantor distance */
    /*
     * d = -b log_e[1 - (p/b)]
     *
     * assuming equal frequencies for all nt, f(nt) = 1/4 and b=3/4
     *
     */

    if (dist >= 0.75)
      return -dist;
    else
      return -(19. / 20.) * log(1. - (dist / (19. / 20.)));
}


// Return the corrected JC distance.
// In case of too large a p-distance >= 0.75, such that a correction is impossible, it returns the negarive p distance.
// If the effective_len is 0, (<0.5 to allow small floating point inacuracies), a distance of -FLT_MAX is returned.
// => A negative distance indicates a problem, the value reveals the reason.
// The calling program has to check for negative distances and deal with them.
// A distance >= 75% should not occur in real data sets. If such a distance is found this hints at a problem,
// e.g. a sequence overlap which is too short.
// TODO: The parameters gapwt and abmig do not make sense at the moment. The only reasonable values are 0 for both parameters.
//       But if the only reasonable values are 0, the parameters can be removed.

inline double dist_jukes_cantor_distPROTEIN(const char *s1, const char *s2, int slen, int &len_no_gaps, double gapwt=0, int ambig=0)
{
    double matches;
    int gaps, k;
    double dist;

    /* compute unambiguoug sum of matches and gaps */
    matches = gaps = 0;
    for (k = 0; k < slen; k++) {
	if ((is_gap(s1[k])) || (is_gap(s2[k])))
	    gaps++;
	/** I added this. Ambigs should not be treated as mismatches. **/
	else if (is_aa_ambig(s1[k]) || is_aa_ambig(s2[k]))
	  ++gaps;
	//	else if (ambig)
	//	    matches += nt_ambig_match(s1[k], s2[k]);
	else if (s1[k] == s2[k])
	    matches += 1.0;
    }

    /* avoid divide by zero */
    double effective_len = ((double) (slen - gaps) + ((double)gaps * gapwt));
    
    len_no_gaps = slen-gaps;

    if (effective_len < 0.5)
      return -FLT_MAX;

    /* compute uncorrected distance */
    /*
     * p = proportion of different nt
     * n = number of nt examined
     * k = number of nt pairs that differ
     * p ~= k/n = 1 - (matches/npos)
     */
    dist = 1.0 - (matches / effective_len );

    /* compute Jukes-Cantor distance */
    /*
     * d = -b log_e[1 - (p/b)]
     *
     * assuming equal frequencies for all nt, f(nt) = 1/4 and b=3/4
     *
     */

    if (dist >= 0.75)
      return -dist;
    else
      return -(19. / 20.) * log(1. - (dist / (19. / 20.)));
}


// Returns -1       if distance cannot be computed.
// Returns -FLT_MAX if number of residues is 0.

//**************************************
// http://emboss.sourceforge.net/apps/release/6.2/emboss/apps/distmat.html
//**************************************
// P = transitions/npos
// Q = transversions/npos

// npos - number of positions scored

// GC1 = GC fraction in sequence 1
// GC2 = GC fraction in sequence 2
// C = GC1 + GC2 - 2*GC1*GC2

// distance = -C ln(1-P/C-Q) - 0.5(1-C) ln(1-2Q)

// Reference:
// K. Tamura, Mol. Biol. Evol. 1992, 9, 678. 

inline double dist_tamura_dnadist(const char *s1, const char *s2, int slen)
{
    int k;
    int matches, transitions, transversions, npos;
    int gc1, at1, gc2, at2;
    double dist, P, Q, theta1, theta2, C;
    const char *nt = "AGCTU";
    const char *pur = "AGR";
    const char *pyr = "CTUY";
    const char *other_ambig = "NKHMVSBWD";

    /* compute transitions, transversions and GC content */
    matches = transitions = transversions = npos = 0;
    gc1 = at1 = gc2 = at2 = 0;
    for (k = 0; k < slen; k++) {
        /* this check is for compatibility with DISTMAT but 
	   shouldn't be done as it forces ignoring valid
	   nucleotides without a corresponding counterpart
	   in the other sequence *** JR *** */
        if (strchr("ATUGC", s1[k]) && strchr("ATUGC", s2[k])) { 

	/* check GC1 content */
	if (strchr("GCS", s1[k]))
	    gc1++;
	else if (strchr("ATUW", s1[k]))
	    at1++;
	/* check GC2 content */
	if (strchr("GCS", s2[k]))
	    gc2++;
	else if (strchr("ATUW", s2[k]))
	    at2++;

	}/* anything else cannot be identified as G+C and is ignored */

	/* gaps are ignored */
	if ((is_gap(s1[k])) || (is_gap(s2[k])))
	    continue;
	else if (is_DNA_ambig_non_RY(s1[k]) || is_DNA_ambig_non_RY(s2[k]))
	  continue;
	else if (strchr(nt, s1[k]) == strchr(nt, s2[k]))
	    /* both are the same unambiguous nucleotide */
	    matches++, npos++;
	else {
	    /* count transitions and transversions */
	    if (strchr(pur, s1[k]) && strchr(pur, s2[k]))
		/* R -> R */
		transitions++, npos++;
	    else if (strchr(pyr, s1[k]) && strchr(pyr, s2[k]))
		/* Y -> Y */
		transitions++, npos++;
	    else
		/* R -> Y */
		transversions++, npos++;
	}
    }
    
    /*
     * avoid divide by zero: if no matches return max. distance;
     * plus,if npos != o ==> there is at least one non-ambiguous nt
     * on each sequence and therefore (gc?+ac?) cannot be zero
     */
    if (npos == 0)
      return -FLT_MAX;

    /*
     * P=transitions/npos;
     * Q=transversions/npos;
     * theta1 = GC fraction in sequence 1
     * theta2 = GC fraction in sequence 2
     * C = theta1 + theta2 - 2*theta1*theta2
     * 
     * distance = -C ln(1-P/C-Q) - 0.5(1-C) ln(1-2Q)
     */
    P = (double) transitions / (double) npos;
    Q = (double) transversions / (double) npos;
    theta1 = (double) gc1 / (double) (gc1 + at1);
    theta2 = (double) gc2 / (double) (gc2 + at2);
    C = theta1 + theta2 - (2 * theta1 * theta2);

    if ((2. * Q) > 1. || P/C + Q > 1)
      return -1.;	/* cannot use Tamura's distance */

    /* compute Tamura distance */
    dist = -C * log(1 - P / C - Q) - (0.5 * (1 - C) * log(1 - (2 * Q)));

    return dist;
}

// Returns -1 if distance cannot be computed.
// Returns -FLT_MAX if number of residues is 0.
// Also provides len_no_gaps, that is the number of non-gap residues in interval.
inline double dist_tamura_dnadist(const char *s1, const char *s2, int slen, int &len_no_gaps)
{
    int k;
    int matches, transitions, transversions, npos;
    int gc1, at1, gc2, at2;
    double dist, P, Q, theta1, theta2, C;
    const char *nt = "AGCTU";
    const char *pur = "AGR";
    const char *pyr = "CTUY";

    /* compute transitions, transversions and GC content */
    matches = transitions = transversions = npos = 0;
    gc1 = at1 = gc2 = at2 = 0;
    for (k = 0; k < slen; k++) {
        /* this check is for compatibility with DISTMAT but 
	   shouldn't be done as it forces ignoring valid
	   nucleotides without a corresponding counterpart
	   in the other sequence *** JR *** */
        if (strchr("ATUGC", s1[k]) && strchr("ATUGC", s2[k])) { 

	/* check GC1 content */
	if (strchr("GCS", s1[k]))
	    gc1++;
	else if (strchr("ATUW", s1[k]))
	    at1++;
	/* check GC2 content */
	if (strchr("GCS", s2[k]))
	    gc2++;
	else if (strchr("ATUW", s2[k]))
	    at2++;

	}/* anything else cannot be identified as G+C and is ignored */

	/* gaps are ignored */
	if ((is_gap(s1[k])) || (is_gap(s2[k])))
	    continue;
	else if (is_DNA_ambig_non_RY(s1[k]) || is_DNA_ambig_non_RY(s2[k]))
	  continue;
	else if (strchr(nt, s1[k]) == strchr(nt, s2[k]))
	    /* both are the same unambiguous nucleotide */
	    matches++, npos++;
	else {
	    /* count transitions and transversions */
	    if (strchr(pur, s1[k]) && strchr(pur, s2[k]))
		/* R -> R */
		transitions++, npos++;
	    else if (strchr(pyr, s1[k]) && strchr(pyr, s2[k]))
		/* Y -> Y */
		transitions++, npos++;
	    else
		/* R -> Y */
		transversions++, npos++;
	}
    }
    
    /*
     * avoid divide by zero: if no matches return max. distance;
     * plus,if npos != o ==> there is at least one non-ambiguous nt
     * on each sequence and therefore (gc?+ac?) cannot be zero
     */

    len_no_gaps = npos;
    if (npos == 0)
      return -FLT_MAX;
    
    /*
     * P=transitions/npos;
     * Q=transversions/npos;
     * theta1 = GC fraction in sequence 1
     * theta2 = GC fraction in sequence 2
     * C = theta1 + theta2 - 2*theta1*theta2
     * 
     * distance = -C ln(1-P/C-Q) - 0.5(1-C) ln(1-2Q)
     */
    P = (double) transitions / (double) npos;
    Q = (double) transversions / (double) npos;
    theta1 = (double) gc1 / (double) (gc1 + at1);
    theta2 = (double) gc2 / (double) (gc2 + at2);
    C = theta1 + theta2 - (2 * theta1 * theta2);

    if ((2. * Q) > 1. || P/C + Q > 1)
      return -1.;	/* cannot use Tamura's distance */

    /* compute Tamura distance */
    dist = -C * log(1 - P / C - Q) - (0.5 * (1 - C) * log(1 - (2 * Q)));

    return dist;
}

inline void all_DNA_distances(CSequences2 &seqs, std::vector<float> *dists, int mode,
			      unsigned min_overlap, unsigned &sum_length,
			      float &sum_scores_of_pairs, float &sum_scores_of_positions,
			      unsigned &num_good_pairs, unsigned &num_undefined_pairs)
{
  unsigned num_good_pairs_local      = 0;
  unsigned num_undefined_pairs_local = 0;

  sum_length              = 0;
  sum_scores_of_pairs     = 0;
  num_good_pairs          = 0;
  sum_scores_of_positions = 0;

  if (min_overlap < 1)
    min_overlap = 1;

  unsigned seqNum = seqs.GetTaxaNum();
  unsigned posNum = seqs.GetPosNum();

  if (posNum == 0)
  {
    std::cerr << "Error: sequences passed to the CDists_collection object do not seem to have equal lengths in "
      "all_DNA_distances. It could also be that they have been read with the remove gaps option.\n";
    std::exit(-97);
  }
  
  // CSequences2 is a class for aligned sequences, but we check whether all sequences
  // have equal length before we proceed.
  //  if (!seqs.equal_length_of_all_sequences() )
  //    return false;

  unsigned  i,j;
  int       tmp;
  unsigned  l;
  unsigned  ambig;
  unsigned  gap;
  unsigned  u;

  float d;

  if (dists != NULL)
  {
    dists->reserve( (seqNum-1)*seqNum/2 + 1 );

    for (i=0; i<seqNum; ++i)
    {
      for (j=0; j<i; ++j)
      {
	if (mode == 128) // Identity, no gaps
	{
	  nucStrictIdentity(seqs.get_Seq_Data(i),seqs.get_Seq_Data(j), posNum, tmp, l, ambig, gap, u);

	  // if (tmp == 0 && l>= min_overlap)
	  // {
	  //   std::cerr << "FLAG RAISED (1).\n";

	  //   faststring s1 = seqs.get_Seq_Data(i);
	  //   faststring s2 = seqs.get_Seq_Data(j);
	  
	  //   std::cerr << s1 << '\n';
	  //   std::cerr << s2 << '\n';
	  // }
 
	  if (l >= min_overlap)
	  {
	    ++num_good_pairs_local;
	    d = (float) tmp/l;
	    sum_scores_of_pairs      += d;
	    sum_scores_of_positions  += tmp;
	    sum_length += l;
	  }
	  else
	  {
	    ++num_undefined_pairs_local;
	    d = -FLT_MAX;

	    faststring s1 = seqs.get_Seq_Data(i);
	    faststring s2 = seqs.get_Seq_Data(j);
	    
	    std::cerr << s1 << '\n';
	    std::cerr << s2 << '\n';
	  }
	  dists->push_back(d);
	  //	std::cout << "(" << i << "," << j << "): " << d << '\n'; 
	  //	std::cout << "Adding(1): " << d << '\n';
	  //	std::cout << seqs.get_Seq_Data(i) << '\n';
	  //	std::cout << seqs.get_Seq_Data(j) << '\n';
	  //	std::cout << "with tmp: " << tmp << " l: " << l << '\n';
	  //	std::cout << "New sum: " << sum << '\n';
	  //	std::cout << sum << '\n';
	}
	else if (mode == 130) // Distance, no gaps
	{
	  nucStrictIdentity(seqs.get_Seq_Data(i),seqs.get_Seq_Data(j), posNum, tmp, l, ambig, gap, u);
	  if (l >= min_overlap)
	  {
	    ++num_good_pairs_local;
	    d = 1.0f - (float) tmp/l;
	    sum_scores_of_pairs      += d;
	    sum_scores_of_positions  += tmp;	  
	    sum_length += l;
	  }
	  else
	  {
	    ++num_undefined_pairs_local;
	    d = -FLT_MAX;

	    // std::cerr << "INSPECT UNDEFINED (2):\n";

	    // faststring s1 = seqs.get_Seq_Data(i);
	    // faststring s2 = seqs.get_Seq_Data(j);
	    
	    // std::cerr << s1 << '\n';
	    // std::cerr << s2 << '\n';
	    
	  }
	  //	std::cout << "(" << i << "," << j << "): " << d << '\n'; 
	  //	std::cout << "Adding(1): " << d << '\n';
	  //	std::cout << seqs.get_Seq_Data(i) << '\n';
	  //	std::cout << seqs.get_Seq_Data(j) << '\n';
	  //	std::cout << "with tmp: " << tmp << " l: " << l << '\n';
	  //	std::cout << "New sum: " << sum << '\n';
	  //	std::cout << sum << '\n';
	  dists->push_back(d);
	}
      } // END for (j=0; j<i; ++j)
    }   // END for (i=0; i<seqNum; ++i)
  }
  else // Compute sum only.
  {
    for (i=0; i<seqNum; ++i)
    {
      for (j=0; j<i; ++j)
      {
	if (mode == 128) // Identity, no gaps
	{
	  nucStrictIdentity(seqs.get_Seq_Data(i),seqs.get_Seq_Data(j), posNum, tmp, l, ambig, gap, u);
	  if (l >= min_overlap)
	  {
	    ++num_good_pairs_local;
	    d                    =  (float) tmp/l;
	    sum_scores_of_pairs  += d;
	    sum_length           += l;
	  }
	  else
	  {
	    ++num_undefined_pairs_local;

	    // std::cerr << "INSPECT UNDEFINED (3):" << '\n';

	    // faststring s1 = seqs.get_Seq_Data(i);
	    // faststring s2 = seqs.get_Seq_Data(j);

	    // std::cerr << s1 << '\n';
	    // std::cerr << s2 << '\n';
	  }
	  //	std::cout << "(" << i << "," << j << "): " << d << '\n'; 
	  //	std::cout << "Adding(1): " << d << '\n';
	  //	std::cout << seqs.get_Seq_Data(i) << '\n';
	  //	std::cout << seqs.get_Seq_Data(j) << '\n';
	  //	std::cout << "with tmp: " << tmp << " l: " << l << '\n';
	  //	std::cout << "New sum: " << sum << '\n';
	  //	std::cout << sum << '\n';
	}
	else if (mode == 130) // Distance, no gaps
	{
	  nucStrictIdentity(seqs.get_Seq_Data(i),seqs.get_Seq_Data(j), posNum, tmp, l, ambig, gap, u);
	  if (l>= min_overlap)
	  {
	    ++num_good_pairs_local;
	    d = 1.0f - (float) tmp/l;
	    sum_scores_of_pairs  += d;
	    sum_length += l;
	  }
	  else
	  {
	    ++num_undefined_pairs_local;

	    // std::cerr << "INSPECT UNDEFINED (4):" << '\n';

	    // faststring s1 = seqs.get_Seq_Data(i);
	    // faststring s2 = seqs.get_Seq_Data(j);
	    
	    // std::cerr << s1 << '\n';
	    // std::cerr << s2 << '\n';
	  }
	  //	std::cout << "(" << i << "," << j << "): " << d << '\n'; 
	  //	std::cout << "Adding(1): " << d << '\n';
	  //	std::cout << seqs.get_Seq_Data(i) << '\n';
	  //	std::cout << seqs.get_Seq_Data(j) << '\n';
	  //	std::cout << "with tmp: " << tmp << " l: " << l << '\n';
	  //	std::cout << "New sum: " << sum << '\n';
	  //	std::cout << sum << '\n';
	}
      } // END for (j=0; j<i; ++j)
    }   // END for (i=0; i<N; ++i)
  }
  num_good_pairs      = num_good_pairs_local;
  num_undefined_pairs = num_undefined_pairs_local;
}


// Overload for sliding window identity analysis:

// Could late be unified with above (refactoring) by using default parameters:
// window_size = 0, window_offset = 0.

inline void all_DNA_distances(CSequences2 &seqs, std::vector<float> *dists, int mode,
			      unsigned min_overlap,       unsigned &sum_length,
			      float &sum_scores_of_pairs, float &sum_scores_of_positions,
			      unsigned &num_good_pairs,   unsigned &num_undefined_pairs,
			      unsigned window_size,       unsigned window_offset)
{
  unsigned num_good_pairs_local      = 0;
  unsigned num_undefined_pairs_local = 0;

  sum_length              = 0;
  sum_scores_of_pairs     = 0;
  num_good_pairs          = 0;
  sum_scores_of_positions = 0;

  if (min_overlap < 1)
    min_overlap = 1;

  unsigned seqNum = seqs.GetTaxaNum();
  unsigned posNum = seqs.GetPosNum();

  if (window_size == 0 || window_size > posNum)
   window_size = posNum;

  if (window_offset == 0 || window_offset > posNum)
   window_offset = posNum;
    
  if (posNum == 0)
  {
    std::cerr << "Error: sequences passed to the CDists_collection object do not seem to have equal lengths in all_DNA_distances. "
                 "It could also be that they have been read with the remove gaps option." << '\n';
    std::exit(-97);
  }

  // Values of window_size and window_offset larger than a quater of the maximum value are not allowed.
  // Note that these values need to be converted to ints below.

  // CSequences2 is a class for aligned sequences, but we check whether all sequences
  // have equal length before we proceed.
  //  if (!seqs.equal_length_of_all_sequences() )
  //    return false;

  unsigned  i,j;
  int       tmp;
  unsigned  l;
  unsigned  ambig;
  unsigned  gap;
  unsigned  u;

  float     d;
  unsigned  N = posNum - window_size + 1;

  if (dists != NULL)
  {
    dists->reserve( ((seqNum-1)*seqNum/2)*(1+(posNum-window_size)/window_offset)+1 );

    for (i=0; i<seqNum; ++i)
    {
      for (j=0; j<i; ++j)
      {
	if (mode == 128) // Identity, no gaps
	{
	  for (unsigned pos = 0; pos < N; pos += window_offset)
	  {
	    nucStrictIdentity(seqs.get_Seq_Data(i),seqs.get_Seq_Data(j), window_size, tmp,
			      l, ambig, gap, u, pos);

	    // Should be a rare case:
	    // if (tmp == 0 && l>= min_overlap)
	    // {
	    //   std::cerr << "FLAG RAISED (2)." << '\n';
	    //   std::cerr << "i,j,pos:     " << i+1 << " " << j+1 << " " << pos+1 << '\n';
	    //   std::cerr << "identity, l: " << tmp << " " << l   << '\n'; 
	    //   std::cerr << "window_size: " << window_size << '\n';
	      
	    //   faststring s1 = seqs.get_Seq_Data(i);
	    //   faststring s2 = seqs.get_Seq_Data(j);

	    //   std::cerr << s1.substr(pos, window_size) << '\n';
	    //   std::cerr << s2.substr(pos, window_size) << '\n';
	    // }
	    
	    if (l >= min_overlap)
	    {
	      ++num_good_pairs_local;
	      d = (float) tmp/l;
	      sum_scores_of_pairs      += d;
	      sum_scores_of_positions  += tmp;
	      sum_length += l;
	    }
	    else
	    {
	      ++num_undefined_pairs_local;
	      d = -FLT_MAX;

	      // std::cerr << "INSPECT UNDEFINED (5):" << '\n';

	      // faststring s1 = seqs.get_Seq_Data(i);
	      // faststring s2 = seqs.get_Seq_Data(j);
	    
	      // std::cerr << s1.substr(pos, window_size) << '\n';
	      // std::cerr << s2.substr(pos, window_size) << '\n';
	    }
	    dists->push_back(d);
	    //	std::cout << "(" << i << "," << j << "): " << d << '\n'; 
	    //	std::cout << "Adding(1): " << d << '\n';
	    //	std::cout << seqs.get_Seq_Data(i) << '\n';
	    //	std::cout << seqs.get_Seq_Data(j) << '\n';
	    //	std::cout << "with tmp: " << tmp << " l: " << l << '\n';
	    //	std::cout << "New sum: " << sum << '\n';
	    //	std::cout << sum << '\n';
	  } // END for (pos = 0; pos < N; pos += window_offset)0
	}
	else if (mode == 130) // Distance, no gaps
	{
	  for (unsigned pos = 0; pos < N; pos += window_offset)
	  {
	    nucStrictIdentity(seqs.get_Seq_Data(i),seqs.get_Seq_Data(j),
			      window_size, tmp, l, ambig, gap, u, pos);
	    if (l >= min_overlap)
	    {
	      ++num_good_pairs_local;
	      d = 1.0f - (float) tmp/l;
	      sum_scores_of_pairs      += d;
	      sum_scores_of_positions  += tmp;	  
	      sum_length += l;
	    }
	    else
	    {
	      ++num_undefined_pairs_local;
	      d = -FLT_MAX;

	      // std::cerr << "INSPECT UNDEFINED (6):" << '\n';

	      // faststring s1 = seqs.get_Seq_Data(i);
	      // faststring s2 = seqs.get_Seq_Data(j);

	      // std::cerr << s1.substr(pos, window_size) << '\n';
	      // std::cerr << s2.substr(pos, window_size) << '\n';
	    }
	    //	std::cout << "(" << i << "," << j << "): " << d << '\n'; 
	    //	std::cout << "Adding(1): " << d << '\n';
	    //	std::cout << seqs.get_Seq_Data(i) << '\n';
	    //	std::cout << seqs.get_Seq_Data(j) << '\n';
	    //	std::cout << "with tmp: " << tmp << " l: " << l << '\n';
	    //	std::cout << "New sum: " << sum << '\n';
	    //	std::cout << sum << '\n';
	    dists->push_back(d);

	  } // END for (pos = 0; pos < N; pos += window_offset)
	}
      } // END for (j=0; j<i; ++j)
    }   // END for (i=0; i<seqNum; ++i)
  } // if (dists != NULL)
  else // Compute sum only
  {
    for (i=0; i<seqNum; ++i)
    {
      for (j=0; j<i; ++j)
      {
        if (mode == 128) // Identity, no gaps
	{
	  for (unsigned pos = 0; pos < N; pos += window_offset)
	  {
	    nucStrictIdentity(seqs.get_Seq_Data(i),seqs.get_Seq_Data(j), window_size, tmp,
			      l, ambig, gap, u, pos);
	    if (l >= min_overlap)
	    {
	      ++num_good_pairs_local;
	      d                    =  (float) tmp/l;
	      sum_scores_of_pairs  += d;
	      sum_length           += l;
	    }
	    else
	    {
	      ++num_undefined_pairs_local;

	      // std::cerr << "INSPECT UNDEFINED (7):" << '\n';

	      // faststring s1 = seqs.get_Seq_Data(i);
	      // faststring s2 = seqs.get_Seq_Data(j);

	      // std::cerr << s1.substr(pos, window_size) << '\n';
	      // std::cerr << s2.substr(pos, window_size) << '\n';
	    }
	    //	std::cout << "(" << i << "," << j << "): " << d << '\n'; 
	    //	std::cout << "Adding(1): " << d << '\n';
	    //	std::cout << seqs.get_Seq_Data(i) << '\n';
	    //	std::cout << seqs.get_Seq_Data(j) << '\n';
	    //	std::cout << "with tmp: " << tmp << " l: " << l << '\n';
	    //	std::cout << "New sum: " << sum << '\n';
	    //	std::cout << sum << '\n';
	  } // END for (pos = 0; pos < N; pos += window_offset)
	}
	else if (mode == 130) // Distance, no gaps
	{
	  for (unsigned pos = 0; pos < N; pos += window_offset)
	  {
	    nucStrictIdentity(seqs.get_Seq_Data(i),seqs.get_Seq_Data(j), window_size, tmp,
			      l, ambig, gap, u, pos);
	    if (l>= min_overlap)
	    {
	      ++num_good_pairs_local;
	      d = 1.0f - (float) tmp/l;
	      sum_scores_of_pairs  += d;
	      sum_length += l;
	    }
	    else
	    {
	      ++num_undefined_pairs_local;

	      // std::cerr << "INSPECT UNDEFINED (8):" << '\n';

	      // faststring s1 = seqs.get_Seq_Data(i);
	      // faststring s2 = seqs.get_Seq_Data(j);

	      // std::cerr << s1 << '\n';
	      // std::cerr << s2 << '\n';
	    }
	    //	std::cout << "(" << i << "," << j << "): " << d << '\n'; 
	    //	std::cout << "Adding(1): " << d << '\n';
	    //	std::cout << seqs.get_Seq_Data(i) << '\n';
	    //	std::cout << seqs.get_Seq_Data(j) << '\n';
	    //	std::cout << "with tmp: " << tmp << " l: " << l << '\n';
	    //	std::cout << "New sum: " << sum << '\n';
	    //	std::cout << sum << '\n';
	  } // END for (pos = 0; pos < N; pos += window_offset)
	}
      } // END for (j=0; j<i; ++j)
    }   // END for (i=0; i<N; ++i)
  }
  num_good_pairs      = num_good_pairs_local;
  num_undefined_pairs = num_undefined_pairs_local;
}


class CDists_collection
{
 private:
  std::vector<float> dists;   // Entries in this vector have the order
                              // (1,0),(2,0),(2,1),(3,0),(3,1),(3,2),(4,0),(4,1),....

  std::vector<float> dists_without_undefined;   // Entries are not ordered.

  unsigned char      mode;    // 0:   AA Blosum score, no gaps in score,
                              // 1:   AA Blosum score, gaps in score
                              // 128: DNA identity, no gaps in score,
                              // 129: DNA identity, with gaps in score, // Not implemented
                              // 130: DNA distance, no gaps in score,   // Not implemented
                              // 131: DNA distance, with gaps in score, // Not implemented

  float              gap_open;
  float              gap_ext;
  float              gap_match_score;
  float              factor_front_back_gaps;

  unsigned           numTaxa;
  unsigned           numPos;

  unsigned           min_overlap;

  float              sum_scores_of_pairs;     // Each score is normalized for the length of the paire
  float              sum_scores_of_position;  // Unnormalized sum of scores.
  unsigned           sum_overlap_len;
  unsigned           num_good_pairs;
  unsigned           num_undefined_pairs;

  unsigned           window_size, window_offset;


  static int vector_index_of_pair(int i, int j)
  {
    // We want to compute pairs (1,0),(2,0),(2,1),(3,0),(3,1),(3,2),(4,0),(4,1),....
    // to index values: 0,1,2,3,4,5,......
    // The index can be computed using i*(i-1)/2+j if i>j
    // or the scores are symetric:     j*(j-1)/2+i if j>i

    if (i==j)
      return -1;

    if (i>j)
      return i*(i-1)/2+j;
    else
      return j*(j-1)/2+i;
  }


  public:
  CDists_collection(CSequences2 &seqs, float pgap_open, float pgap_ext, float pgap_match_score, float pfactor_front_back_gaps, int pmode, unsigned pmin_overlap):
    mode(pmode), gap_open(pgap_open), gap_ext(pgap_ext), gap_match_score(pgap_match_score),
    factor_front_back_gaps(pfactor_front_back_gaps),
    numTaxa(seqs.GetTaxaNum()), numPos(seqs.GetPosNum()),
    min_overlap(pmin_overlap),    
    sum_scores_of_pairs(0), sum_scores_of_position(0), sum_overlap_len(0), num_good_pairs(0), num_undefined_pairs(0),
    window_size(0), window_offset(0)
  {
    if (min_overlap == 0)
      min_overlap = 1;
    
    if (mode == 0) // Blosum scores, no gaps
    {
      all_AA_distances(seqs, &dists, min_overlap, sum_overlap_len, sum_scores_of_pairs, sum_scores_of_position, num_good_pairs, num_undefined_pairs);
    }
    else if (mode == 1) // Blosum scores, with gaps
    {
      all_AA_distances_gaps(seqs, &dists, gap_open, gap_ext, gap_match_score, factor_front_back_gaps, min_overlap, sum_overlap_len, sum_scores_of_pairs, sum_scores_of_position, num_good_pairs, num_undefined_pairs);
    }
    else if (mode == 128) // DNA identity, no gaps
    {
      all_DNA_distances(seqs, &dists, mode, min_overlap, sum_overlap_len, sum_scores_of_pairs, sum_scores_of_position, num_good_pairs, num_undefined_pairs);
    }
    else if (mode == 129) // DNA identity, with gaps
    {
      all_DNA_distances(seqs, &dists, mode, min_overlap, sum_overlap_len, sum_scores_of_pairs, sum_scores_of_position, num_good_pairs, num_undefined_pairs);
    }
    else if (mode == 130) // DNA distance, no gaps
    {
      all_DNA_distances(seqs, &dists, mode, min_overlap, sum_overlap_len, sum_scores_of_pairs, sum_scores_of_position, num_good_pairs, num_undefined_pairs);
    }
    else if (mode == 131) // DNA distance, with gaps
    {
      all_DNA_distances(seqs, &dists, mode, min_overlap, sum_overlap_len, sum_scores_of_pairs, sum_scores_of_position, num_good_pairs, num_undefined_pairs);
    }
    else
    {
      dists.clear();
    }
  }

  CDists_collection(CSequences2 &seqs, float pgap_open, float pgap_ext, float pgap_match_score, float pfactor_front_back_gaps, int pmode, unsigned pmin_overlap, unsigned pwindow_size, unsigned pwindow_offset):
    mode(pmode), gap_open(pgap_open), gap_ext(pgap_ext), gap_match_score(pgap_match_score),
    factor_front_back_gaps(pfactor_front_back_gaps),
    numTaxa(seqs.GetTaxaNum()), numPos(seqs.GetPosNum()),
    min_overlap(pmin_overlap),
    sum_scores_of_pairs(0), sum_scores_of_position(0), sum_overlap_len(0), num_good_pairs(0), num_undefined_pairs(0),
    window_size(pwindow_size), window_offset(pwindow_offset)
  {
    if (min_overlap == 0)
      min_overlap = 1;
    
    if (mode == 0) // Blosum scores, no gaps
    {
      all_AA_distances(seqs, &dists, min_overlap, sum_overlap_len, sum_scores_of_pairs, sum_scores_of_position, num_good_pairs, num_undefined_pairs);
    }
    else if (mode == 1) // Blosum scores, with gaps
    {
      all_AA_distances_gaps(seqs, &dists, gap_open, gap_ext, gap_match_score, factor_front_back_gaps, min_overlap, sum_overlap_len, sum_scores_of_pairs, sum_scores_of_position, num_good_pairs, num_undefined_pairs);
    }
    else if (mode == 128) // DNA identity, no gaps
    {
      all_DNA_distances(seqs, &dists, mode, min_overlap, sum_overlap_len, sum_scores_of_pairs, sum_scores_of_position, num_good_pairs, num_undefined_pairs, window_size, window_offset);
    }
    else if (mode == 129) // DNA identity, with gaps
    {
      all_DNA_distances(seqs, &dists, mode, min_overlap, sum_overlap_len, sum_scores_of_pairs, sum_scores_of_position, num_good_pairs, num_undefined_pairs, window_size, window_offset);
    }
    else if (mode == 130) // DNA distance, no gaps
    {
      all_DNA_distances(seqs, &dists, mode, min_overlap, sum_overlap_len, sum_scores_of_pairs, sum_scores_of_position, num_good_pairs, num_undefined_pairs, window_size, window_offset);
    }
    else if (mode == 131) // DNA distance, with gaps
    {
      all_DNA_distances(seqs, &dists, mode, min_overlap, sum_overlap_len, sum_scores_of_pairs, sum_scores_of_position, num_good_pairs, num_undefined_pairs, window_size, window_offset);
    }
    else
    {
      dists.clear();
    }
  }

  // Verify before using this function. 
  void print_debug(std::ostream &os)
  {
    os.setf(std::ios::fixed);
    os.precision(2);

    os << "Satus of CDists_collection object:\n";
    os << "Data type:       " << (mode < 128 ? "Data type: AA":"Data type: DNA") << '\n';
    os << "Mode (int):      " << (int)mode     << '\n';
    os << "Minimum overlap: " << min_overlap   << '\n';
    os << "Gap mode:   ";
    switch (mode)
    {
       case 0:   os << "AA without gaps"           << '\n'; break;
       case 1:   os << "AA with gaps"              << '\n'; break;
       case 128: os << "DNA identity without gaps" << '\n'; break;
       case 129: os << "DNA identity with gaps"    << '\n'; break;
       case 130: os << "DNA distance without gaps" << '\n'; break;
       case 131: os << "DNA distance with gaps"    << '\n'; break;
    }
    os << "gap_open:   " << gap_open  << '\n';
    os << "gap_ext:    " << gap_ext << '\n';
    os << "gap_match_score:        " << gap_match_score << '\n';
    os << "factor_front_back_gaps: " << factor_front_back_gaps << '\n';
    os << "numTaxa:                " << numTaxa  << '\n';
    os << "numPos:                 " << numPos   << '\n';
    
    if (window_size != 0)
    {
      os << "Sliding window: window size:   " << window_size << '\n';
      os << "Sliding window: window offset: " << window_offset << '\n';
    }

    os << "Number of dists in vector: " << dists.size() << '\n';
    os << "Raw dists: ";

    for (unsigned i=0; i<dists.size(); ++i)
      os << dists[i] << ", ";
    os << '\n';

    os << "Defined dists:\n";
    const std::vector<float> &mydists = get_dists_without_undefined();

    for (unsigned i=0; i<mydists.size(); ++i)
      os << mydists[i] << ", ";
    os << '\n';
    
    if (window_size == 0)
    {
      os << "Dists as lower matrix:\n";

      int k=0;
      for (unsigned i=0; i<numTaxa; ++i)
      {
	for (unsigned j=0; j<i; ++j)
	{
	  os << dists[k] << " ";
	  ++k;
	}
	os << '\n';
      }
      os << "\n\n";
    }

    os << "sum_scores_of_pairs:     " << sum_scores_of_pairs       << '\n';
    os << "mean_per_pair_score:     " << get_mean_per_pair_score() << '\n';
    os << "num_good_pairs:          " << num_good_pairs            << '\n';
    os << "num_undefined_pairs:     " << num_undefined_pairs       << '\n';
    os << "sum_scores_of_position:  " << sum_scores_of_position    << '\n';
    os << "sum_overlap_len:         " << sum_overlap_len           << '\n';
    os << "get_mean_per_pos_score:  " << get_mean_per_pos_score()  << '\n';
  }

  double get_dist_sum()
  {
    return sum_scores_of_pairs;
  }

  double get_mean_per_pair_score()
  {
    return sum_scores_of_pairs/num_good_pairs;
  }

  unsigned get_number_of_good_distance_pairs()
  {
    return num_good_pairs;
  }

  unsigned get_number_of_undefined_distance_pairs()
  {
    return num_undefined_pairs;
  }
  
  double get_mean_per_pos_score()
  {
    return sum_scores_of_position/(sum_overlap_len);
  }

  bool successfully_populated()
  {
    return dists.size()>0;
  }

  float get_dist(int i, int j)
  {
    if (i < j)
      return dists[vector_index_of_pair(i,j)];
    else if (i > j)
      return dists[vector_index_of_pair(j,i)];
    else
      return -FLT_MAX;
  }

  void get_min_max(float &min_score, float &max_score)
  {
    const std::vector<float> mydists = get_dists_without_undefined();
    vec_min_max(mydists, min_score, max_score);
  }

  void get_mean_sd(double &mean_score, double &sd_score)
  {
    const std::vector<float> mydists = get_dists_without_undefined();
    vec_mean_sd(mydists, mean_score, sd_score);
  }

  void get_mean_sd(double &mean_score, double &sd_score, float &sum, float &sum_squares)
  {
    const std::vector<float> mydists = get_dists_without_undefined();
    vec_mean_sd(mydists, mean_score, sd_score, sum, sum_squares);
  }



  //######################
  // These element functions are too specific. The user should do this.
  void get_median_quartile(double &Q1_score, double &median_score, double &Q3_score)
  {
    const std::vector<float> mydists = get_dists_without_undefined();
    vec_median_quartile_method3(mydists, Q1_score, median_score, Q3_score);
  }

  void get_median_quartile_outlier_bounds(double &Q1_score, double &median_score, double &Q3_score, double &O_lower, double &O_upper)
  {
    const std::vector<float> mydists = get_dists_without_undefined();
    vec_median_quartile_outlier_bounds_method3(mydists, Q1_score, median_score, Q3_score, O_lower, O_upper);
  }
  //######################


  const std::vector<float>& get_dists()
  {
    return dists;
  }

  const std::vector<float>& get_dists_without_undefined()
  {
    if (dists_without_undefined.empty())
    {
      dists_without_undefined.reserve(dists.size());
      unsigned i, N=dists.size();
      for (i=0; i<N; ++i)
	if (dists[i] != -FLT_MAX)
	  dists_without_undefined.push_back(dists[i]);
    }
    return dists_without_undefined;
  }

  // The vector returned can contain -FLT_MAX values, indicating undefined distances.
  void get_dists(unsigned i, std::vector<float> &v)
  {
    v.clear();
    
    unsigned j;
    for (j=0; j<numTaxa; ++j)
    {
      if (i!=j)
      {
	v.push_back(get_dist(i,j));
      }
    }
  }

  float get_mean_dist_to_others(unsigned i)
  {
    unsigned j;
    float    sum=0;
    float    d;
    unsigned count = 0;
    
    //    std::cout << "Averaging: ";

    for (j=0; j<numTaxa; ++j)
    {
      if (i==j)
      {
	//+0	
      }
      else
      {
	d = get_dist(i,j);
	if (d != -FLT_MAX)
	{
	  //	std::cout << "(" << i << "," << j << "):" << d << ",";
	  sum += d;
	  ++count;
	}
      }
    }
    return sum/count;
  }

  void get_vector_of_mean_dists_to_others(std::vector<float> &v)
  {
    unsigned j;

    v.clear();
    for (j=0; j<numTaxa; ++j)
    {
      v.push_back(get_mean_dist_to_others(j));
    }
  }


};




#endif

