
clear all;
load hybridLyap.mat;
A = A';
Z = frlibPrg(A,b,c,K);
prg = Z.ReducePrimal('dd');
[x,y] = prg.Solve();
xO = prg.RecoverPrimal(x);
pass  = prg.CheckPrimal(x);
if ~(pass)
    error('Test case failed')
end
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



































