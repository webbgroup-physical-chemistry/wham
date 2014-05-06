#ifndef BW_OPTION_PARSER_H_
#define BW_OPTION_PARSER_H_

#include <boost/program_options.hpp>
#include <iterator>
#include <iostream>
#include "bw_types.h"

void bw_option_parser(int argc, char *argv[], bw_options &options);

void bw_print_options(bw_options &options);

#endif