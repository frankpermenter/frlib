clear all;

testPass = [];

load hybridLyap.mat;
A = A';
prg = frlibPrg(A,b,c,K);

%diagonal
display('Checking diagonal fr')
prgD = prg.ReducePrimal('d');
[x,y] = prgD.Solve();
x0 = prgD.RecoverPrimal(x);
pass  = prg.CheckPrimal(x0);
testPass(end+1) = pass & all( prgD.K.s == [6 56 11 1 1 0 11 1 1 0 11 11]);
if ~(testPass(end))
    warning('Test case failed')
end
display('Checking diagonally dominant fr')
%diagonally dominant
prgDD = prgD.ReducePrimal('dd');
[x,y] = prgDD.Solve();
x0 = prgDD.RecoverPrimal(x);
x1 = prgD.RecoverPrimal(x0);
pass  = prg.CheckPrimal(x1);
testPass(end+1) = pass;

if ~(testPass(end))
    warning('Test case failed')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load testDual.mat;
prg = frlibPrg(A,[],c,K);
prgDD = prg.ReduceDual('dd');
[~,y] = prgDD.Solve();

testPass(end+1)  = prg.CheckDual(y) & prgDD.K.f == 3 & prgDD.K.s == 2;

if ~(testPass(end))
    warning('Test case failed')
end


prgD = prg.ReduceDual('d');
[~,y] = prgD.Solve();

testPass(end+1)  = prg.CheckDual(y) & prgDD.K.f == 3 & prgDD.K.s == 2;

if ~(testPass(end))
    warning('Test case failed')
end


pass = runHorn;
testPass = [testPass,pass];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('%d of %d tests passed \n',sum(testPass),length(testPass))



































