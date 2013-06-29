Reduction of Dual Problems:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%pass in sedumi formatted A,b,c,K
prg = frlibPrg(A,b,c,K);

%reduce the dual
prgR = prg.ReduceDual('dd');

%solve. solutions to reduced dual are solutions to unreduced
[~,yFr]=prgR.Solve();
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Reduction of Primal Problems
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%pass in sedumi formatted A,b,c,K
Z = frlibPrg(A,b,c,K);

%Reduce primal
prg = Z.ReducePrimal('dd');
[x,y] = prg.Solve();

%Recover solution to unreduced problem
xO = prg.RecoverPrimal(x);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



Inner approximations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%The following commands reduced an SDP using the non-negative
%diagonal  matrices as an inner approximation of the PSD cone

prg = Z.ReduceDual('d');
prg = Z.ReducePrimal('d');


%The following commands reduced an SDP using the diagonally
%dominant matrices as an inner approximation of the PSD cone

prg = Z.ReduceDual('dd');
prg = Z.ReducePrimal('dd');



