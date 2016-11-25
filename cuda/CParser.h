/* 
 * File:   CParser.h
 * Author: Petr Klejch (klejcpet@fit.cvut.cz)
 *
 * Created on 14. března 2016, 9:30
 */

#ifndef CPARSER_H
#define	CPARSER_H

using namespace std;

#include <string>
#include <complex>


class CParser {
public:
    CParser();
    static void parse(float* &reals, float* &imags, int &N);
private:

};

#endif	/* CPARSER_H */

