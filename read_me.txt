NFORMATION:
"PUR-and-RUR-compare" contains five files:"myAlgorithm3_PUR.mpl" and "myAlgorithm3_PUR.mpl"  and "Exampleresults_PUR.mw" and "Exampleresults_RUR.mw" and "10Examples.mw".


These algorithms are implemented by Jian Pan (panjian201806@163.com).


The myAlgorithm3_PUR are based on the paper: "Jian Pan, et al : Computing PUR of Zero-dimensional Ideals." 


The Rouillier_RUR is based on the paper: "Rouillier F, Solving zero-dimensional systems through the Rational Univariate
  Representation, {\it Appl. Algebra Engrg. Comm. Comput.}, 1999,  {\bf 9}(5): 433--461." 


PROCEDURES:
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
(n:the number of variables; ns,rv : Normal set; G: Groebner basis; tord: term ordering (tdeg(x,y,z)))



6. GetPUR  (F::{list, set, PolynomialIdeal}, tord::(:-ShortMonomialOrder))
It is used to solve the linear symbolic systems and get PUR.
 

Note on the normal set assumption
Our algorithm assumes that the normal set (basis of the quotient ring) contains all variables ​. This condition is satisfied for many ideals arising from applications and is essential for extracting a univariate representation for each variable. If a variable is missing from the normal set, the algorithm will simply omit its corresponding polynomial from the output, while the representations for the other variables may still be correct. Users should verify that their ideal meets this assumption, or adapt the code accordingly.