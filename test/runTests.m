clear all;
testPass = [];

load hybridLyap.mat;
A = A';
prg = frlibPrg(A,b,c,K);

%diagonal
display('Checking diagonal fr')
prgD = prg.ReducePrimal('d');
opts.useQR = 1;
[x,y] = prgD.Solve(opts);
x0 = prgD.RecoverPrimal(x);
pass  = prg.CheckPrimal(x0,10^-4);
testPass(end+1) = pass & all( prgD.K.s == [6 56 11 1 1 0 11 1 1 0 11 11]);
if ~(testPass(end))
    error('Test case failed')
end

%diagonally dominant
display('Checking diagonally dominant fr')
prgDD = prgD.ReducePrimal('dd');
opts.useQR = 1;
[x,y] = prgDD.Solve(opts);
x1 = prgDD.RecoverPrimal(x);
x2 = prgD.RecoverPrimal(x1);
pass  = prg.CheckPrimal(x2,10^-4);
testPass(end+1) = pass;

if ~(testPass(end))
    error('Test case failed')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load testDual.mat;
prg = frlibPrg(A,[],c,K);
prgDD = prg.ReduceDual('dd');
opts.useQR = 1;
[~,y] = prgDD.Solve(opts);

testPass(end+1)  = prg.CheckDual(y,10^-4) & prgDD.K.s == 2;

if ~(testPass(end))
    error('Test case failed')
end


prgD = prg.ReduceDual('d');
[~,y] = prgD.Solve(opts);

testPass(end+1)  = prg.CheckDual(y,10^-4) & prgDD.K.s == 2;

if ~(testPass(end))
    error('Test case failed')
end


pass = runHorn;
testPass = [testPass,pass];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('%d of %d tests passed \n',sum(testPass),length(testPass))

