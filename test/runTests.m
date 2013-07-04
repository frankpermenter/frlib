
clear all;
load hybridLyap.mat;
A = A';
Z = frlibPrg(A,b,c,K);

%diagonal
prg = Z.ReducePrimal('d');
[x,y] = prg.Solve();
xO = prg.RecoverPrimal(x);
pass  = prg.CheckPrimal(x);
pass = pass & all( prg.K.s ==[6 56 11 1 1 0 11 1 1 0 11 11]);
if ~(pass)
    error('Test case failed')
end

%diagonally dominant
prg = Z.ReducePrimal('dd');
[x,y] = prg.Solve();
xO = prg.RecoverPrimal(x);
pass  = prg.CheckPrimal(x);
pass = pass & all( prg.K.s == [6 34 8 1 1 0 8 1 1 0 9 7]);
if ~(pass)
    error('Test case failed')
end
prg
break;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load testDual.mat;


Z = frlibPrg(A,[],c,K);
prg = Z.ReduceDual('dd');
[x,y] = prg.Solve();
pass  = prg.CheckDual(y);
if ~(pass)
    error('Test case failed')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



display('Test cases passed!')



































