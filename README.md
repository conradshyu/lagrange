# Lagrange Interpolating Polynomial for Free Energy Estimates

Written by [Conrad Shyu](mailto:conradshyu@hotmail.com)<br>
*release 2, September 8, 2009*<br>
*release 3, March 11, 2014*<br>
*last updated on October 30, 2026*

## Synopsis

This archive contains the implementation of Lagrange interpolating polynomial for estimates of free energy
differences using thermodynamic integration. It also contains an example data file to jumpstart use of the
software. The implementation makes heavy use of STL (standard template library) which requires a C++ compiler.
This program is self-contained and should not require any additional libraries.

To compile the program, simply type:

`make`

The make program will invoke the Makefile and compiles the source code. Alternatively, to compile the files
manually, issue the command:

`g++ -I. driver.cpp lagrange.cpp -o lagrange -lm`

The project contains the following files:

| Filename | Description |
| --- | --- |
| `driver.cpp` | the user interface program |
| `lagrange.cpp` | the main implementation of Lagrange interpolation polynomial |
| `lagrange.h` | the header file for Lagrange interpolation polynomial |
| `example.csv` | example file of simulations data from thermodynamic integration |
| `Makefile` | the makefile to compile the source code |
| `README.md` | this file |
| `LICENSE` | GNU General Public License, version 3 or later |
| `.gitignore` | patterns for generated build artifacts |

The program constructs a Lagrange interpolating polynomial that fits the thermodynamic integration data and
reports its integral alongside a trapezoidal estimate over the first and last input lambda values. It can also
evaluate the polynomial at regularly spaced points and write those estimates to a file.

The Lagrange interpolation program takes one required and two optional command-line parameters:
`input_file`, `plot_file`, and `data_points`. The following sections describe their formats and meanings.

## `input_file`

The input file must be plain text. Each data row contains a lambda value followed by its ensemble average from the
simulation. Separate the two values with a comma, semicolon, or whitespace. Lines whose first character is `#` are
ignored. Use distinct lambda values in ascending order so the interpolation and integration are well-defined.
For example, the following is the contents of the example input file:

```
0.0,  51.49866347
0.1,  23.92508775
0.2,  10.35390700
0.3,   2.58426990
0.4,  -2.18351656
0.5,  -5.41745387
0.6,  -7.62452181
0.7,  -9.25455804
0.8, -10.45592989
0.9, -11.39244138
1.0, -12.12433704
```

To run the analysis program with only the required parameter, type:

`lagrange example.csv`

The program will output the interpolating polynomial and the estimate of free energy using the Lagrange
interpolating polynomial and trapezoidal rule. For example, the example file should generate the following output:

```
Degree, Coefficients
     0,      51.49866347
     1,    -440.06771988
     2,    2724.49944200
     3,  -16355.52325204
     4,   74074.45611540
     5, -226406.51705359
     6,  454849.19406749
     7, -590543.95494381
     8,  476225.01322752
     9, -216633.75369268
    10,   42443.03080908

Free energy difference
 Lagrange: 0.65644834
Trapezoid: 1.02220063
```

## `plot_file`

This parameter specifies the name of the file to which estimates from the interpolating polynomial are written.
This parameter is optional; if omitted, no plot data is generated. The output contains two comma-separated
columns: lambda and the estimated ensemble average. Estimates are generated over the fixed range 0 to 1, inclusive.
For example, to write plot data using the default number of intervals, run:

`lagrange example.csv plot.csv`

The file, example.csv, is the input file that contains the thermodynamic integration data. The file, plot.csv, is
the output file that contains the estimates from the interpolating polynomial.

## `data_points`

This optional parameter specifies the number of intervals used to sample the polynomial over 0 to 1. The output
contains one more row than the interval count because both endpoints are included. If omitted, the program uses
one fewer interval than the number of input rows, producing the same number of output rows as input rows. For
example, to generate 101 rows (100 intervals), run:

`lagrange example.csv plot.csv 100`

## Notes

It is strongly recommended to limit the number of lambda values to less than a dozen. High degree of Lagrange
polynomials suffers greatly from numerical instability. The analysis from the program can fluctuate widely and
becomes completely unreliable (see the reference for more information).

Please report any problems or send comments to [me](mailto:conradshyu@hotmail.com).

> **Note**: Lagrange and Newton polynomials are mathematically equivalent. The numerical differences are only
notable after the eighth digit. The Newton polynomial is generally considered more computationally efficient than
Lagrange.

## Reference

C. Shyu and F.M. Ytreberg (2011). Accurate estimation of solvation free energy using polynomial fitting
techniques. *Journal of Computational Chemistry*, 32(1): 134-141, doi: 10.1002/jcc.21609.

## License

This software is distributed under the GNU General Public License, version 3 or (at your option) any later version.
See [`LICENSE`](LICENSE) for the full license text.
