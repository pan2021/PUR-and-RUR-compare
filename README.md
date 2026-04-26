# PUR-and-RUR-compare
Implementation of the PUR algorithm for zero-dimensional polynomial ideals in Maple. Comparison with the classic RUR algorithm. Includes ideal deflation and linear algebra-based symbolic computation methods. All source codes are .mpl files. ALL results are .mw files.
##INFORMATION:

"PUR-and-RUR-compare" contains six files:"Algorithm3_PUR.mpl" and "Rouillier_RUR.mpl"  and "Example_results_PUR.mw" and "Example_results_RUR(Rouillier).mw" and "Example_results_RURIrr.mw" and "Example_results_RURVars.mw".
These algorithms are implemented by Jian Pan (panjian201806@163.com).
The Algorithm3_PUR are based on the paper: "Jian Pan, et al : Computing PUR of Zero-dimensional Ideals." 
The Rouillier_RUR is based on the paper: "Rouillier F, Solving zero-dimensional systems through the Rational Univariate Representation, {\it Appl. Algebra Engrg. Comm. Comput.}, 1999,  {\bf 9}(5): 433--461." 

##PROCEDURES:
All the algorithm contains the following sub-algorithms:

1.GetM(ns::list, rv::table, G::{list, set, PolynomialIdeal}, tord::(:-ShortMonomialOrder)) 
It is used to compute the MultiplicationMatrix of every variable.
(ns,rv : Normal set; G: Groebner basis; tord: term ordering (tdeg(x,y,z)))

2 .BM (s, N, x)
3 .LCMpoly (p, q, x)
4 . minimalpoly  (A::Matrix, startIndex::posint := RowDimension(A))
Sub-algorithms 2-4 are used to compute minimalpolynomial of matrix.  
(s: sequence; N:half of the length of sequence;p,q: polynomial)
5.  GetSquarefree (f::polynom)
It is used to compute the squarefree of the minimalpolynomial.

6. IdealDeflation  (n::posint, ns::list, rv::table, G::{list, set, PolynomialIdeal}, tord::(:-ShortMonomialOrder)) 
It is used to compute an ideal of breadth at most one and return minimal polynomial and transformation basis matrix.
(n:the number of variables; ns,rv : Normal set; G: Groebner basis; tord: term ordering (graded reverse lexicographic order)

6. GetPUR  (F::{list, set, PolynomialIdeal}, tord::(:-ShortMonomialOrder))
It is used to solve the linear symbolic systems and get PUR.
 
##Note on the normal set assumption

Our algorithm works correctly for ideals where the normal set (basis of the quotient ring) either contains all variables or is missing exactly one variable. In both cases, the corresponding univariate representations for all variables, including the one missing from the normal set, are fully computed and output.

Our algorithm works correctly for ideals where the normal set (basis of the quotient ring) either contains all variables or is missing exactly one variable. In both cases, the corresponding univariate representations for all variables, including the one missing from the normal set, are fully computed and output. This covers many well-known benchmark ideals, such as cyclic 3, cyclic 5, which often have one variable missing from the normal set under graded reverse lexicographic order.  If two or more variables are missing from the normal set, the algorithm will still output correct representations for all variables present in the normal set, but will omit the corresponding polynomials for the missing variables. The algorithm is not designed to handle such cases, and the representations for the missing variables are not guaranteed to be computed correctly. Users are encouraged to verify that their ideal meets the above assumptions before use.

##Comparison with Related Work
This implementation is compared with two representative methods for zero-dimensional polynomial system solving:
the RUR algorithm by Rouillier (1999) and the recent method proposed by Xiao et al. (2025). For the algorithm and experiments of Xiao et al. (2025), please refer to their original paper and official GitHub repository:

• Paper:  Xiao S, Zeng G, Irredundant Decomposition of the Radicals of Polynomial Ideals Based on Rational Univariate Representations. {\it Journal of Systems Science & Complexity}, 2025, {\bf 38}(6): 277-2778.

• GitHub: https://github.com/realalgebra/RURIrr

##Experimental Results
* "Example_results_PUR.mw"  — Output and runtime results of Algorithm3.
* "Example_results_RUR(Rouillier).mw" — Output and runtime results of Rouillier's algorithm (1999) .
*  "Example_results_RURIrr.mw" and "Example_results_RURVars.mw" —Output and runtime results of Xiao's algorithm (2025).
All experiments are conducted on the same of benchmark examples.
