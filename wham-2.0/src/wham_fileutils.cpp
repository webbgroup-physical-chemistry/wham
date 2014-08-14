#include "wham_fileutils.h"

int get_n_lines(const std::string &filename)
{
    std::string line;
    int lines = 0;
    std::ifstream file(filename.c_str());
    if (file.is_open())
    {
        while (file.good())
        {
            getline(file, line);
            if (not line.empty())
            {
                lines ++;
            }
        }
    }
    else if (!file)
    {
        std::cerr << "\nError opening " << filename << "." << std::endl;
        std::exit(1);
    }
    return lines;
}

bool fexists( const std::string &filename )
{
	std::ifstream file (filename.c_str());
    if (file.good()) {
        file.close();
        return true;
    }
    else {
        file.close();
        return false;
    }
}

std::string backup( std::string &filename )
{
    int n = 0, nchars = 0;
    std::string newname = filename;
    if (fexists(filename))
    {
        while (fexists(newname))
        {
            char buffer[1054];
            newname.clear();
            nchars = sprintf(buffer,"#%s.%d#",filename.c_str(),n);
            newname = buffer;
            n++;
            if (n>100)
            {
                std::cout << "Will not make more than 100 copies of " << filename << std::endl;
                std::cout << "Enter a different file name to backup to:" << std::endl;
                getline(std::cin,filename);
                newname = filename;
            }
        }
        std::cout << "Backing up " << filename << " to " << newname << std::endl;
        rename(filename.c_str(),newname.c_str());
    }
    std::cout << "Writing to " << filename << std::endl;
    return filename;
}
