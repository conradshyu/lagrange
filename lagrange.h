/*
 * lagrange.h
 *
 * Lagrange interpolating polynomials for free energy estimates
 * Copyright (C) 2014   Conrad Shyu (conradshyu@hotmail.com)
 * Department of Physics, University of Idaho
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * Author's comments
 * -----------------
 * note: the implementation has been verified to produce the correct result with a
 * simple integration on December 11, 2007
 *
 * written by Conrad Shyu (conradshyu@hotmail.com)
 *
 * first created on December 4, 2007
 * revised on December 11, 2007
 * revised on September 1, 2008 for production release
 * revised on March 6, 2014
*/
#ifndef _LAGRANGE_H
#define _LAGRANGE_H

#include <list>
#include <cmath>
#include <bitset>
#include <vector>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <cstring>
#include <iostream>

const unsigned int LAGRANGE_DEGREE = sizeof( unsigned int ) * 8;

typedef struct
{
    double x;   // positions on the x-axis
    double y;   // values on the y-axis, y=f(x)
} stLAGRANGE;

class Lagrange
{
public:
    Lagrange();
    Lagrange( const std::list<stLAGRANGE>& );
    Lagrange( const std::vector<double>&, const std::vector<double>& );
    ~Lagrange() {};

    bool GetEstimate( const std::string&, unsigned int ) const;

    double DoIntegral( bool = false ) const;
    double DoQuadrature( bool = false ) const;

    const std::vector<double>& GetPolynomial( bool = false ) const;
    const std::list<stLAGRANGE>& LoadData( const std::list<stLAGRANGE>& );
    const std::list<stLAGRANGE>& LoadData( const std::vector<double>&, const std::vector<double>& );

private:
    std::list<stLAGRANGE> sample;
    std::vector<double> factor;

    void ClearData();           // reset the data structure
    void DoPolynomial();        // construct the lagrange polynomial
    const std::vector<double>& GetPermute( std::vector<double>& ) const;
};  // class definition for lagragne interpolating polynomial

#endif  // _LAGRANGE_H
