/*
 * driver.cpp
 *
 * Lagrange interpolating polynomials for free energy estimates
 * Copyright (C) 2014   Conrad Shyu (cornadshyu at hotmail.com)
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
 * the template for the main driver program
 *
 * written by Conrad Shyu (conradshyu at hotmail.com)
 *
 * revised on September 1, 2008
 * revised on March 6, 2014
*/

#include <lagrange.h>

const char* DELIMIT_TOKEN       = "\n\t,; ";
const unsigned int MAX_BUFFER   = 1024;

/*
 * open the file and load the data into memory
*/
unsigned int LoadData(
    std::ifstream& _file,
    std::list<stLAGRANGE>& _sample )
{
    char buffer[ MAX_BUFFER ];
    stLAGRANGE data;

    _sample.clear();

    while ( _file.getline( buffer, MAX_BUFFER ) )
    {
        if ( buffer[ 0 ] == '#' )
        {
            continue;
        }   // skip the comments

        /*
         * windows nags about the deprecated function
         * linux somehow doesn't complain about it
        */
        data.x = atof( strtok( buffer, DELIMIT_TOKEN ) );   // lambda
        data.y = atof( strtok( NULL, DELIMIT_TOKEN ) );     // dg/dl value
        _sample.push_back( data );
    }   // parse the file

    return( _sample.size() );
}   // end of LoadData()

/*
 * the main program
*/
int main( int argc, char** argv )
{
    // print out GLP licencse V3
    std::cout << "Copyright (C) 2014 Conrad Shyu (conrads@uidaho.edu)" << std::endl;
    std::cout << "This free software and comes with ABSOLUTELY NO WARRANTY." << std::endl << std::endl;

    if ( argc < 2 )
    {
        std::cout << argv[ 0 ] << " input_file [plot_file [data_points]]" << std::endl;
        std::cout << " input_file: file contains thermodynamic integration data" << std::endl;
        std::cout << "  plot_file: file for the plot data [optional]" << std::endl;
        std::cout << "data_points: number of data points for plot [optional]" << std::endl;
        std::cout << std::endl << "see readme.txt for more information" << std::endl;
        return( 0 );
    }   // check the number of parameters

    std::ifstream ifs( argv[ 1 ], std::ios::in );

    if ( ifs.bad() )
    {
        std::cout << "failed to open the file " << argv[ 1 ] << std::endl; return( 1 );
    }   // make sure the file can be opened successfully

    std::list<stLAGRANGE> sample; LoadData( ifs, sample ); ifs.close();

    Lagrange a( sample ); a.GetPolynomial( true );

    printf( "\nFree energy difference\n Lagrange: %.8f\nTrapezoid: %.8f\n",
        a.DoIntegral(), a.DoQuadrature() );

    if ( argc > 2 )
    {
        a.GetEstimate(
            argv[ 2 ],      // filename for plot data
            ( argc > 3 ) ? atoi( argv[ 3 ] ) : ( sample.size() - 1 ) );
    }   // need to write the estimates

    return( 0 );
}   // end of main
