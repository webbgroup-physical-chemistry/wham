#ifndef PARSE_OPTIONS_H_
#define PARSE_OPTIONS_H_

#include <boost/program_options.hpp>
#include <iterator>
#include <iostream>
#include "wham_types.h"

void wham_options(int argc, char *argv[], t_options &options);
void multidimensional_option(std::vector<std::string> option,
                             t_options &options,
                             float defaultX,
                             std::string valuename,
                             std::vector<float> &X);
void print_options(t_options &options);

#endif