/*
 * lagrange.cpp
 *
 * Lagrange interpolating polynomials for free energy estimates
 * Copyright (C) 2014   Conrad Shyu (conradshyu at hotmail.com)
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
 * note: the implementation has been verified to produce the correct
 * result with a simple integration on December 11, 2007
 *
 * written by Conrad Shyu (conradshyu at hotmail.com)
 *
 * first created on December 4, 2007
 * revised on December 11, 2007
 * revised on September 1, 2008 for production release
 * revised on March 6, 2014
*/

#include <lagrange.h>

/*
 * default class constructor
*/
Lagrange::Lagrange()
{
    ClearData();
}   // end of class constructor

/*
 * class constructor
*/
Lagrange::Lagrange(
    const std::list<stLAGRANGE>& _sample )
{
    LoadData( _sample );
}   // end of class constructor

/*
 * class constructor
*/
Lagrange::Lagrange(
    const std::vector<double>& _x,
    const std::vector<double>& _y )
{
    LoadData( _x, _y );
}   // end of class constructor

/*
 * reset and initialize essential variables
*/
const std::list<stLAGRANGE>& Lagrange::LoadData(
    const std::list<stLAGRANGE>& _sample )
{
    stLAGRANGE unit; ClearData();

    for ( std::list<stLAGRANGE>::const_iterator i = _sample.begin(); !( i == _sample.end() ); i++ )
    {
        unit.x = ( *i ).x; unit.y = ( *i ).y; sample.push_back( unit );
    }   // save a local copy of the data

    // perform interpolation with lagrange polynomials
    DoPolynomial(); return( sample );
}   // end of LoadData()

/*
 * reset and initialize essential variables
*/
const std::list<stLAGRANGE>& Lagrange::LoadData(
    const std::vector<double>& _x,
    const std::vector<double>& _y )
{
    stLAGRANGE unit; ClearData();

    for ( unsigned int i = 0; i < _x.size(); ++i )
    {
        unit.x = _x[ i ]; unit.y = _y[ i ]; sample.push_back( unit );
    }   // save a local copy of the data

    // perform interpolation with lagrange polynomials
    DoPolynomial(); return( sample );
}   // end of LoadData()

/*
 * clear all contents
*/
void Lagrange::ClearData()
{
    sample.clear(); factor.clear();
}   // end of ClearData()

/*
 * calculate the coefficients for the polynomials
 *
 * note: the calculations of coefficients for the polynomials have been verified
 * correctly on December 10, 2007
*/
const std::vector<double>& Lagrange::GetPermute(
    std::vector<double>& _x ) const
{
    unsigned int bit_value = ( 0x1 << _x.size() );
    std::vector<double> term( ( _x.size() + 1 ), 0.0 );
    std::bitset<LAGRANGE_DEGREE> permute;
    double unit;

    for ( unsigned int i = 0; i < bit_value; ++i )
    {
        unit = 1.0; permute = i;

        for ( unsigned int j = 0; j < _x.size(); ++j )
        {
            unit *= ( permute[ j ] ) ? ( -1.0 * _x[ j ] ) : 1.0;
        }   // calculate combinatoric terms; -a_1 * -a_2 * ...

        term[ permute.count() ] += unit;
    }   // iterate through all possible combinations of factors

    _x = term; return( _x );    // invoke copy constructor and overwrite contents
}   // end of GetPermute()

/*
 * set the constant terms for the lagrange interpolating polynomial
 * set the coefficients and factors for the lagrange interpolating polynomial
 * note: the calculuations of polynomial constants have been verified correctly on
 * December 10, 2007
*/
void Lagrange::DoPolynomial()
{
    unsigned int i = 0;
    std::vector<double> term;
    std::vector<double> constant( sample.size(), 1.0 ); // polynomial constant factors
    factor.resize( sample.size(), 0.0 );
    std::list<stLAGRANGE>::iterator p;
    std::list<stLAGRANGE>::iterator q;

    for ( i = 0, p = sample.begin(); !( p == sample.end() ); p++, ++i )
    {
        term.clear();

        for ( q = sample.begin(); !( q == sample.end() ); q++ )
        {
            if ( p == q )
            {
                continue;
            }   // calculate the factor only if i != j

            constant[ i ] *= ( ( *p ).x - ( *q ).x );
            term.push_back( ( *q ).x );     // construct list of coefficients
        }   // accumulate the constants ( x_i - x_j )

        constant[ i ] = ( ( *p ).y / constant[ i ] );
        GetPermute( term );     // calculate the coefficients for each variable

        for ( unsigned int j = 0; j < term.size(); ++j )
        {
            factor[ j ] += ( constant[ i ] * term[ j ] );
        }   // accumulate the factors with constants and coefficients
    }   // calculate the constants first

    term = factor; factor.clear();

    for ( unsigned int j = 0; j < term.size(); ++j )
    {
        factor.push_back( term[ ( term.size() - j - 1 ) ] );
    }   // need to reverse the order of coefficients
}   // end of DoPolynomial()

/*
 * reset and initialize the data points
 * perform integration on the lagrange polynomial
*/
double Lagrange::DoIntegral(
    bool _print ) const
{
    double lower = ( sample.front() ).x;
    double upper = ( sample.back() ).x;
    double area = 0.0;
    double power;

    for ( unsigned int i = 0; i < factor.size(); ++i )
    {
        power = static_cast<double>( i + 1.0 );
        area += ( ( pow( upper, power ) / power ) * factor[ i ] -
            ( pow( lower, power ) / power ) * factor[ i ] );
    }   // perform the integration with given upper and lower bounds

    if ( _print )
    {
        printf( "area under the curve: %.8f\n", area );
    }   // print out the integration result

    return( area );
}   // end of DoIntegral()

/*
 * calculate the area under the curve using quadrature
*/
double Lagrange::DoQuadrature(
    bool _print ) const
{
    std::list<stLAGRANGE>::const_iterator a = sample.begin();
    std::list<stLAGRANGE>::const_iterator b = sample.begin(); b++;
    double area = 0.0;

    while ( !( b == sample.end() ) )
    {
        area += ( ( *b ).y + ( *a ).y ) * 0.5 * ( ( *b ).x - ( *a ).x );
        a = b; b++;
    }   // iterate through the entire list

    if ( _print )
    {
        printf( "area under the curve: %.8f\n", area );
    }   // print out the integration result

    return( area );
}   // end of DoQuadrature()

/*
 * return the lagrange interpolating polynomial
*/
const std::vector<double>& Lagrange::GetPolynomial(
    bool _print ) const
{
    if ( _print )
    {
        printf( "Degree, Coefficients\n" );

        for ( unsigned int i = 0; i < factor.size(); ++i )
        {
            printf( "%6d, %.8f\n", i, factor[ i ] );
        }   // print out the lagrange polynomial
    }   // print out the coefficients of the polynomial, if necessary

    return( factor );
}   // end of GetPolynomial()

/*
 * get the estimate of f(x) with a given value
*/
bool Lagrange::GetEstimate(
    const std::string& _file,
    const unsigned int _step ) const
{
    std::ofstream ofs( _file.c_str(), std::ios::trunc );

    if ( ofs.bad() )
    {
        std::cout << "file " << _file << "cannot be opened" << std::endl; return( false );
    }   // make sure the file stream has been opened successfully

    char buffer[ 80 ]; double y = 0.0; double x = 0.0;
    double step = 1.0 / static_cast<double>( _step );

    for ( unsigned int s = 0; !( s > _step ); ++s )
    {
        for ( unsigned int i = 0; i < factor.size(); ++i )
        {
            y += ( pow( x, static_cast<double>( i ) ) * factor[ i ] );
        }   // substitute the given value and calculate the estimate

        sprintf( buffer, "%.4f, %.8f", x, y ); ofs << buffer << std::endl;
        y = 0.0; x += step;
    }   // iterate through the entire interval

    ofs.close(); return( true );
}   // end of GetEstimate()

/*
 * the main driver program for lagrange interpolation polynomial
*/
/*
int main( void )
{
    list<stLAGRANGE> sample;
    stLAGRANGE unit;

    sample.clear();

    unit.x = 0.0; unit.y =  51.49866347; sample.push_back( unit );
    unit.x = 0.1; unit.y =  23.92508775; sample.push_back( unit );
    unit.x = 0.2; unit.y =  10.35390700; sample.push_back( unit );
    unit.x = 0.3; unit.y =   2.58426990; sample.push_back( unit );
    unit.x = 0.4; unit.y =  -2.18351656; sample.push_back( unit );
    unit.x = 0.5; unit.y =  -5.41745387; sample.push_back( unit );
    unit.x = 0.6; unit.y =  -7.62452181; sample.push_back( unit );
    unit.x = 0.7; unit.y =  -9.25455804; sample.push_back( unit );
    unit.x = 0.8; unit.y = -10.45592989; sample.push_back( unit );
    unit.x = 0.9; unit.y = -11.39244138; sample.push_back( unit );
    unit.x = 1.0; unit.y = -12.12433704; sample.push_back( unit );

    Lagrange a( sample );
    a.DoIntegral( true ); a.GetPolynomial( true );

    return( 1 );
}   // end of main()
*/
