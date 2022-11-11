/*
 * ©2022 The Regents of the University of California.  All rights reserved.
 */

#include <iostream> // std

#ifndef AD_ERROR_HPP
#define AD_ERROR_HPP

class ADerror {
    ADerror(char *m, char *f, int l) : message(m), file(f), line(l) {}

public:
    std::string message;
    std::string file;
    int line;

    void ADcatch(ADerror *err) {
        // VERBOSE
//        std::cout << "ERROR:" << err->message << "\n in file " << err->file << " at line \n" << err->line;
        // SHORTER
        std::cout << "ERROR: \n in file " << err->file << " at line \n" << err->line;
        exit(1);
    }
};

#define throw(message) throw(ADerror(message,__FILE__,__LINE__));

#endif