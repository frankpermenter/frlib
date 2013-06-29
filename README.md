To reduce an SDP in dual form, pass sedumi formatted inputs A,b,c,K to the function

prg = frlibPrg(A,b,c,K);

Then call to reduce the dual using diagonally dominant ('dd') or diagonal ('d') inner approximations of the PSD cone.

prgR = prg.ReduceDual('dd');

To solve the reduced SDP, call:

[~,yFr]=prgR.Solve();

To reduce the primal, call

prgR = prg.ReducePrimal('dd');

[xFr]=prgR.Solve();

To recover solutions x0 to the unreduced problem call

xO = prgR.RecoverPrimal(xFr);


