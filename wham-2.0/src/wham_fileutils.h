#ifndef WHAM_FILEUTILS_H_
#define WHAM_FILEUTILS_H_

#include <iostream>
#include <fstream>
#include <sstream>
#include "wham_types.h"

int get_n_lines( const std::string &filename );

bool fexists( const std::string &filename );

std::string backup( std::string &filename );


#endif