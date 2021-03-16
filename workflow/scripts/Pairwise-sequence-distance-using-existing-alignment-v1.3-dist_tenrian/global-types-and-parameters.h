#ifndef GLOBAL_TYPES_AND_PARAMETERS_H
#define GLOBAL_TYPES_AND_PARAMETERS_H

#include <cstdio>
#include <iostream>
#include <vector>
#include <string>
#include <cstring>
#include "faststring2.h"
#include <climits>
#include <fstream>

#define DEBUG

#define PROGNAME "Pairwise-alignment-identity-checker"
#define VERSION  "1.3"

extern faststring                       global_input_fasta_file;
extern faststring                       global_data_type_string;

extern unsigned                         global_min_overlap_len;
extern unsigned                         global_window_size;
extern unsigned                         global_window_offset;
extern unsigned                         global_verbosity;

#define macromax(x,y) ((x)<(y) ? (y) : (x))
#define macromin(x,y) ((x)<(y) ? (x) : (y))

#ifdef  DEBUG
#define DEBUGOUT1(x)        std::cerr << x                << '\n';
#define DEBUGOUT2(x,y)      std::cerr << x << y           << '\n';
#define DEBUGOUT3(x,y,z)    std::cerr << x << y << z      << '\n';
#define DEBUGOUT4(x,y,z,w)  std::cerr << x << y << z << w << '\n';
#else
#define DEBUGOUT1(x)
#define DEBUGOUT2(x,y)
#define DEBUGOUT3(x,y,z)
#define DEBUGOUT4(x,y,z,w)
#endif


void good_bye_and_exit(FILE *of, int);
void good_bye_and_exit(std::ostream&, int);
void init_param();
void read_and_init_parameters(int argc, char** argv, std::ostream &);
//void print_parameters(FILE*, const char *s);
void print_parameters(std::ostream&, const char *s);
void print_calling_command_line(FILE*,unsigned, char**);
void print_calling_command_line(std::ostream&,unsigned, char**);

#endif
