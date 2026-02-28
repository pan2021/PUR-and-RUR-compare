
with(LinearAlgebra):
with(Groebner):
with(ListTools): with(PolynomialTools):
with(PolynomialIdeals):
##GetM :Compute Multiplicationtensor  ##
 GetM := proc (ns, rv::table, G::{list, set, PolynomialIdeal}, tord::(:-ShortMonomialOrder))
 local M, i; M := Vector(nops(ns), 0);
 for i to nops(ns) do M[i] := Groebner:-MultiplicationMatrix(ns[i], ns, rv, G, tord) end do; 
 return M; end proc: 
 ##GetVtr1 :compute Trace ##
 GetVtr1 := proc (M::(Vector(Matrix)))
 local i, vtr1; vtr1 := Vector(ArrayNumElems(M)); 
 for i to ArrayNumElems(M) do vtr1[i] := LinearAlgebra:-Trace(M[i]); end do;
 return vtr1; end proc:
##GetQ1 :compute Q1 for obtaining the number of zeros of ideal##
 GetQ1 := proc (M::(Vector(Matrix)), vtr1::Vector) 
 local Q1, i, j, d; d := ArrayNumElems(M);
 Q1 := Matrix(d, d, shape = symmetric); 
 for i to d do for j to i do Q1[i, j] := LinearAlgebra:-DotProduct(Vector(M[i][j, 1 .. -1]), vtr1); end do; end do; 
 return Q1; end proc: 
 ##GetChiT :compute Characteristic polynomial##
 GetChiT := proc (vtr1::Vector, p::polynom, ns::list, rv::table, G::{list, set, PolynomialIdeal}, tord::(:-ShortMonomialOrder)) 
 local mm, N, v, V, a, i; 
 mm := Groebner:-MultiplicationMatrix(p, ns, rv, G, tord);
 N := Vector[row](nops(ns)+1, 0); 
 V := Vector[row](nops(ns)+1, 0); N[1] := ArrayNumElems(vtr1); v := Vector[row](nops(ns), 0); v[1] := 1; V[1] := v; a := Vector[row](nops(ns)+1, 0); a[1] := 1;
 for i to nops(ns) do 
 v := Vector[row](v) . mm; V[i+1] := v; N[i+1] := LinearAlgebra:-DotProduct(Vector(v), Vector(vtr1)); end do; 
 for i to nops(ns) do 
 a[i+1] := -add(a[i-j+1]*N[j+1], j = 1 .. i)/i; 
 end do; 
 return PolynomialTools:-FromCoefficientList(convert(a, list), t, termorder = reverse), V; end proc:
##GetSquarefree :compute  the squarefree part of Characteristic polynomial
 GetSquarefree := proc (f::polynom) 
 local s, r, i; 
 s := sqrfree(f); 
 r := 1;
 for i to nops(s[2]) do r := r*s[2][i][1]; end do; 
 r := expand(r); return r/PolynomialTools:-CoefficientList(r, op(indets(r)))[-1]; end proc: 
 ##computing RUR  ##
 GetGt1T := proc (ChiT::polynom)
 local diffChiT; 
 diffChiT := diff(ChiT, op(indets(ChiT))); 
 return simplify(diffChiT/gcd(diffChiT, ChiT)); end proc:

 GetGtVT := proc (v::polynom, V::Vector, ChiTSquarefree::polynom, sepele::polynom, vtr1::Vector, ns::list, rv::table, G::{list, set, PolynomialIdeal}, tord::(:-ShortMonomialOrder))
 local mmv, GetHorner, Hi, vtrv, c, c1, d, gtvt, i, j, k; 
 GetHorner := proc (f::polynom, j::nonnegint) local c; if nops(indets(f)) <> 1 then error "The number of indeterminates must be one!"; return  end if; 
 c := PolynomialTools:-CoefficientList(f, op(indets(f))); 
 return PolynomialTools:-FromCoefficientList(c[-j-1 .. -1], op(indets(f))); end proc:
 mmv := Groebner:-MultiplicationMatrix(v, ns, rv, G, tord); vtrv := mmv . Vector[column](vtr1); d := degree(ChiTSquarefree, op(indets(ChiTSquarefree))); 
 gtvt := 0;
 for i from 0 to d-1 do 
 gtvt := gtvt+LinearAlgebra:-DotProduct(Vector(V[i+1]), Vector(vtrv))*GetHorner(ChiTSquarefree, d-1-i); end do;
 return gtvt; end proc:

##GetSeparatingElement :compute Separating Element ##
 GetSeparatingElement := proc (vtr1::Vector, d::posint, n::posint, ns::list, rv::table, G::{list, set, PolynomialIdeal}, tord::(:-ShortMonomialOrder)) 
 local i, j, c, t, D, poly, sfpoly, V, sepeles, vars;
 D := nops(ns);
 vars := PolynomialIdeals:-IdealInfo[Variables](PolynomialIdeals:-PolynomialIdeal(G)); 
 for i from 0 to (1/2)*(n-1)*D*(D-1) do
 c := [1, seq(i^j, j = 1 .. n-1)]; t := add(vars[j]*c[j], j = 1 .. n); 
 poly, V := GetChiT(vtr1, t, ns, rv, G, tord); sfpoly := GetSquarefree(poly); 
 if evalb(degree(sfpoly, op(indets(poly))) = d) then return t, poly, V; end if; end do; return 0, 0, 0; end proc: 
 
 #GetRUR : main function #
 GetRUR := proc (F::{list, set, PolynomialIdeal}, tord::(:-ShortMonomialOrder)) 
 
 local G, d, ns, rv, M, vtr1, n, sepele, ChiT, V, g0, g, vars, i; 
 G := Groebner:-Basis(F, tord); 
 ns, rv := Groebner:-NormalSet(G, tord); 
 M := GetM(ns, rv, G, tord); 
 vtr1 := GetVtr1(M); d := LinearAlgebra:-Rank(GetQ1(M, vtr1)); 
 vars := PolynomialIdeals:-IdealInfo[Variables](PolynomialIdeals:-PolynomialIdeal(F));
 n := nops(vars); sepele, ChiT, V := GetSeparatingElement(vtr1, d, n, ns, rv, G, tord); 
 g0 := simplify((diff(ChiT, t))/gcd(diff(ChiT, t), ChiT)); 
 g := [seq(0, i = 1 .. n)]; 
 for i to n do g[i] := GetGtVT(vars[i], V, GetSquarefree(ChiT), sepele, vtr1, ns, rv, G, tord);  end do; 
 return sepele, ChiT, g0, op(g); end proc:

