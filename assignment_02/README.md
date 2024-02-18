## Assignment_2 README

The pdf report is located in this directory *assignment_02/* in the *report.pdf*
file. I'm not really sure how best to report on the coding portion. All I did
was display the *output.txt* with latex. Should I manually type the results
using latex?

In order to run the code, change directory
```bash
cd code/
```

and compile and run with
```bash
make
```

which outputs to the *output.txt* file. All subroutines used are in the
*LinAl.f90* file and the driver file is *Driver_LinAl.f90* which calls the
required subroutines. I use the subroutine print_mat(A) for printing (which
doesn't require dimensions) but I have a subroutine printMat(A, dims) which
does. The trace is calculated in subroutine traceMat(A, m, trace), column norms
are calculated and printed with subroutine printColumnNorm(A) which uses the
twoNorm(vec,n,norm) subroutine. Gaussian elimination with partial pivoting uses
the partGaussElim(A, B, dimsA, dimsB, SINGLR) subroutine, and reduces the
system to an upper triangular matrix and them uses the subroutine
backSubstitutionU(A,B,X) to find the solution.

I have two LU decomposition routines, and the one I use in the report stores
the matrices L and U inside the matrix A. This is in the subroutine
LU_partial(A, m, s, SINGLR), and it uses the subroutine
backSubstitutionLU(A,m,B,s,X) which unpacks the LU decomposition from A and
then calls the subroutines backSubstitutionL(L,B,Y) and backSubstitution(U,Y,X)
to find the solution.

### A very basic application
The code to plot the plane equation is located in *./code/plane.jl* file. Its
a julia file and uses the Plotly backend as well as the RowEchelon package (I
know I know, I initially just wanted to verify the solutions between julia and
fortran were the same, but then I ran out of time to read the output from the
fortran solution). I use the *P.dat* file to solve the system in fortran. I can
use python or matlab next time if necessary, but if you have julia you can add
the packages to run the code in the following way.

First start the julia repl and type this key to open the package manager
```bash
]
```

then type
```bash
add Plotly
add RowEchelon
```

and run the code
```bash
julia plane.jl
```

this will produce the *../figures/plane_plot.html* file and can be opened in
your web browser. I saved this as *../figures/newplot_.png* for the report.

The system I used to solve was 3x4 and had a free variable *d* for the constant
in the plane equation. I figured I could use any number for *d*, and take the
dot product with the normal vector <a,b,c> and any of the three points *A,B,C*
to obtain the correct value, but that didn't seem to work. In my report I use
the incorrect value of *d=1* for the plot so the points do not lie in the
plane. I'm not sure what I was doing wrong but taking that dot product produced
the value *d = -36* (approximately), but it seems its actually about *d = -12.8* which I manually determined after the deadline. 
