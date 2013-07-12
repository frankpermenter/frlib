frlib is matlab code for partial SDP facial reduction.  It cheaply converts an SDP with no strictly
feasible solution into a smaller, equivalent SDP using inner approximations of the PSD cone.  Better
inner approximations lead to more reductions at the cost of more computation.  Inner approximations
currently supported include non-negative diagonal matrices and PSD diagonally dominant matrices.

To reduce an SDP in dual form, pass sedumi formatted inputs A,b,c,K to the function

prg = frlibPrg(A,b,c,K);

To reduce the dual using diagonally dominant ('dd') or diagonal ('d') inner approximations of the PSD cone, issue the call:

prgR = prg.ReduceDual('dd');

To solve the reduced SDP, call:

[~,y] = prgR.Solve();

The solution y solves the original SDP.

To reduce the primal, call

prgR = prg.ReducePrimal('dd');

To recover solutions x0 to the unreduced primal problem call

[xFr]=prgR.Solve();
x0 = prgR.RecoverPrimal(xFr);


