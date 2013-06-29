To reduce an SDP in dual form, pass sedumi formatted inputs A,b,c,K to the function

prg = frlibPrg(A,b,c,K);

Then call.
prgR = prg.ReduceDual('dd');


[~,yFr]=prgR.Solve();
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Reduction of Primal Problems
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%pass in sedumi formatted A,b,c,K

prg = frlibPrg(A,b,c,K);

%Reduce primal

prgR = prg.ReducePrimal('dd');

[x,y] = prgR.Solve();

%Recover solution to unreduced problem

xO = prgR.RecoverPrimal(x);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



Inner approximations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%The following commands reduced an SDP using the non-negative diagonal  matrices as an inner approximation of the PSD cone

prg = prgR.ReduceDual('d');
prg = prgR.ReducePrimal('d');


%The following commands reduced an SDP using the diagonally dominant matrices as an inner approximation of the PSD cone

prg = prgR.ReduceDual('dd');
prg = prgR.ReducePrimal('dd');



