
with(Groebner): with(ListTools): with(LinearAlgebra): with(PolynomialIdeals): with(PolynomialTools):
##GetM :Compute MultiplicationMatrix ##
GetM := proc (ns::list, rv::table, G::{list, set, PolynomialIdeal}, tord::(:-ShortMonomialOrder)) 
local Mmatrix, i, vars; 
vars := PolynomialIdeals:-IdealInfo[Variables](PolynomialIdeals:-PolynomialIdeal(G)); 
Mmatrix := Vector(nops(vars), 1, []);
 for i to nops(vars) do 
 Mmatrix[i] := Groebner:-MultiplicationMatrix(vars[i], ns, rv, G, tord);
 end do; 
 return Mmatrix; 
 end proc:
 ##Compute minimalpolynomial of MultiplicationMatrix with the BerlekampMassey algorithm##
 BM := proc (s, N, x)
 local C, B, L, m, denom, n, d, i, T;
 if andmap(`=`, s, 0) then return 1 end if;
 C := 1; B := 1; L := 0; m := 1; denom := 1; 
 for n from 0 to 2*N-1 do 
 d := s[n+1]; 
 for i to L do 
 d := d+coeff(C, x, i)*s[n-i+1];
 end do;
 if d = 0 then m := m+1 else T := C; C := expand(C-d*x^m*B/denom); 
 if 2*L <= n then L := n+1-L; B := T; denom := d; m := 1 else m := m+1 end if end if end do; 
 expand(x^L*subs(x = 1/x, C)); end proc: 
 
 LCMpoly := proc (p, q, x) 
 local l; 
 if p = 1 then return q end if; if q = 1 then return p end if; 
 l := lcm(p, q); expand(primpart(l, x)); end proc:
 
 minimalpoly := proc (A::Matrix, startIndex::posint := RowDimension(A))
 local D, g, d, i, base, b, e, powers, s_list, f, local_g, krylov_vectors, result_mat, j, SFg; 
 D := RowDimension(A);
 if startIndex < 1 or D < startIndex then error "startIndex must be between 1 and %1", D end if; 
 g := 1; d := 0; krylov_vectors := []; e := [seq(Row(IdentityMatrix(D), k), k = 1 .. D)];
 for base from startIndex to D while d < D do
 b := e[base]; powers := Array(1 .. 2*D); powers[1] := b; 
 for j from 2 to 2*D do 
 powers[j] := VectorMatrixMultiply(powers[j-1], A);
 end do; 
 local_g := 1; 
 for i to D while degree(local_g, T) < D do
 s_list := [seq(powers[j][i], j = 1 .. 2*D)];
 f := BM(s_list, D, T); 
 if f <> 1 then local_g := LCMpoly(local_g, f, T) end if end do; 
 g := LCMpoly(g, local_g, T); d := degree(g, T); 
 if base = startIndex then krylov_vectors := [seq(powers[j], j = 1 .. D)] end if end do; 
 result_mat := Matrix(D, D, []); 
 if d = D then result_mat := Matrix(D, D, proc (row, col) options operator, arrow; krylov_vectors[D+1-row][col] end proc); return g, result_mat; 
 elif d < D then return g, result_mat; end if; end proc:



  ## GetSquarefree : Compute the squarefree part of the minimalpolynomial   ##
 GetSquarefree := proc (f::polynom) 
 local s, r, i;
 s := sqrfree(f); r := 1;
 for i to nops(s[2]) do r := r*s[2][i][1]; end do; 
 r := expand(r); 
 return r/PolynomialTools:-CoefficientList(r, op(indets(r)))[-1]; end proc:
 
 
 
 ##IdealDeflation : Compute an ideal of breadth at most one ##
 IdealDeflation := proc (n::posint, ns::list, rv::table, G::{list, set, PolynomialIdeal}, tord::(:-ShortMonomialOrder))
 local i, j, k, m, c, t, D, vars, Mlist, M_t, g, Atran, SFg, List_new, SFgX, SFgX1, M_sfg, M_row, A_g, Ns, G_new, ns_new, rv_new, D_new, g1, Atran1;
 vars := PolynomialIdeals:-IdealInfo[Variables](PolynomialIdeals:-PolynomialIdeal(G)); 
 ns_new := ns; rv_new := table(rv); G_new := G; D_new := nops(ns_new); Mlist := GetM(ns_new, rv_new, G_new, tord); 
 for i from 0 to n-1 do 
 c := [1, seq(i^j, j = 1 .. n-1)]; t := add(vars[j]*c[j], j = 1 .. n); M_t := add(Mlist[j]*c[j], j = 1 .. n);
 g, Atran := minimalpoly(M_t); 
 SFg := GetSquarefree(g); List_new := []; 
 if evalb(degree(g, op(indets(g))) = D_new) then return g, Atran, ns_new;
 elif evalb(degree(g, op(indets(g))) < D_new) and degree(SFg, op(indets(SFg))) < degree(g, op(indets(g))) then SFgX := expand(subs(T = t, SFg)); 
 SFgX1 := NormalForm(SFgX, G_new, tord); M_sfg := Groebner:-MultiplicationMatrix(SFgX1, ns_new, rv_new, G_new, tord); M_row := LinearAlgebra:-ReducedRowEchelonForm(M_sfg); Ns := convert(ns_new, Vector[column]); 
 A_g := MatrixVectorMultiply(M_row, Ns); List_new := select(proc (A_g) options operator, arrow; A_g <> 0 end proc,
 [seq(A_g(k, 1), k = 1 .. RowDimension(A_g))]); G_new := [op(G_new), op(List_new)]; G_new := InterReduce(G_new, tord); 
 ns_new, rv_new := NormalSet(G_new, tord); ns_new := Reverse(ns_new); rv_new := table([seq(ns_new[m] = m, m = 1 .. nops(ns_new))]); 
 D_new := nops(ns_new); Mlist := GetM(ns_new, rv_new, G_new, tord);
 elif evalb(degree(g, op(indets(g))) < D_new) and degree(g, op(indets(g))) = degree(SFg, op(indets(SFg))) then 
 ns_new := ns_new; rv_new := table(rv_new); G_new := G_new; D_new := D_new; Mlist := Mlist end if end do;
 for i from n to (1/2)*(3*n-3)*D_new*(D_new-1) do
 c := [1, seq(i^j, j = 1 .. n-1)]; t := add(vars[j]*c[j], j = 1 .. n); M_t := add(Mlist[j]*c[j], j = 1 .. n);
 g1, Atran1 := minimalpoly(M_t); if evalb(degree(g1, op(indets(g1))) = D_new) then return g1, Atran1, ns_new; end if; end do; 
 return 0, 0, 0; end proc: 
 
 
 
  ## GetPUR : main function##
 GetPUR := proc (F::{list, set, PolynomialIdeal}, tord::(:-ShortMonomialOrder))
 local i, j, k, D_0, S, S1, n, P, L, U, V, X, Y, T, t, P1, m, l, G, Ns, ns, rv, Atran, g, vars, g1, ns_new;
 G := Groebner:-Basis(F, tord); ns, rv := Groebner:-NormalSet(G, tord);
 ns := ListTools:-Reverse(ns); rv := table([seq(ns[i] = i, i = 1 .. nops(ns))]);
 vars := PolynomialIdeals:-IdealInfo[Variables](PolynomialIdeals:-PolynomialIdeal(G));
 n := nops(vars); 
 g, Atran, ns_new := IdealDeflation(n, ns, rv, G, tord);
 D_0 := nops(ns_new); Ns := convert(ns_new, Vector[column]); Y := Matrix(D_0, 1, []); l := 0; m := 0;
 for i to D_0 do if 1 < degree(Ns(i, 1), {op(vars)}) then l := l+1; elif degree(Ns(i, 1), {op(vars)}) < 2 then break; end if; end do;
 m := D_0-l; S := Vector(D_0, proc (i) options operator, arrow; T^(D_0-i) end proc); 
 P, L, U := LinearAlgebra:-LUDecomposition(Atran);
 P1 := LinearAlgebra:-Transpose(P); 
 S1 := LinearAlgebra:-MatrixVectorMultiply(P1, S);
 Y := LinearAlgebra:-LinearSolve(L, S1); X := Matrix(D_0, 1, []); 
 for k from D_0 by -1 to D_0-m+1 do for j from k+1 to D_0 do Y(k, 1) := Y(k, 1)-U(k, j)*X(j, 1) end do; 
 X(k, 1) := Y(k, 1)/U(k, k) end do;
 g1 := [seq(0, i = 1 .. m-1)]; for i to m-1 do g1[i] := X(D_0-m+i, 1) end do; return g1, g; 
 end proc:
