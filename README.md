To reduce an SDP in dual form, pass sedumi formatted inputs A,b,c,K to the function

prg = frlibPrg(A,b,c,K);

To reduce the dual using diagonally dominant ('dd') or diagonal ('d') inner approximations of the PSD cone, issue call of the form:

prgR = prg.ReduceDual('dd');

To solve the reduced SDP, call:

[~,yFr]=prgR.Solve();

To reduce the primal, call

prgR = prg.ReducePrimal('dd');


To recover solutions x0 to the unreduced problem call

[xFr]=prgR.Solve();

x0 = prgR.RecoverPrimal(xFr);


