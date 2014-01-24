clear all;

testPass = [];
pass = [];
K.f = 0;
K.l = 3;
K.q = 2;
A = [[0,1,-1,1,0];...
     [0,0,0,0,1]...
     ];

c = [1,0,0,0,0];
prg = frlibPrg(A,[],c,K);
prgRed = prg.ReduceDual('d');
[~,y] = prgRed.Solve();

% -y_1 \ge 0, y_1 \ge 0
% (-y1,-y2) \in SOC(2)
% implies y1 = y2 = 0;
if prgRed.K.l == prg.K.l-2 && prgRed.K.q == 0 && (norm(y) < 10^-12)
    pass(end+1) = 1;
else
    pass(end+1) = 0;
end


%%%%%%%%%
clear K;
K.f = 1;
K.l = 2;
K.q = 2;
K.s = 1;
A = [[0,0,-1,0,0,1];...
     [0,-1,0,1,0,0]...
     ];

c = [0,0,0,0,0,0];

prg = frlibPrg(A,[],c,K);
prgRed = prg.ReduceDual('d');
%[~,y] = prgRed.Solve();

% -y_1 \ge 0, y_1 \ge 0
% (-y1,-y2) \in SOC(2)
% implies y1 = y2 = 0;
if prgRed.K.l == prg.K.l-2 && prgRed.K.q == 0 && prgRed.K.s(1) == 0
    pass(end+1) = 1;
else
    pass(end+1) = 0;
end

%%%%
clear K;
K.f = 1;
K.l = 2;
K.q = 2;
K.s = 1;
A = [[0,0,0,-1,0,1];...
     [0,0,0,0,0,0]...
     ];

c = [0,0,0,0,0,0];

prg = frlibPrg(A,[],c,K);
prgRed = prg.ReduceDual('d');
%[~,y] = prgRed.Solve();

% -y_1 \ge 0, y_1 \ge 0
% (-y1,-y2) \in SOC(2)
% implies y1 = y2 = 0;
if prgRed.K.l == prg.K.l-2 && prgRed.K.q == 0 && prgRed.K.s(1) == 0
    pass(end+1) = 1;
else
    pass(end+1) = 0;
end

%%%%%%
clear K;
K.f = 1;
K.l = 1;
K.q = 3;
K.s = 1;
A = [[0,0,1,1,0,-1];...
     [0,1,0,0,1,0]...
     ];

c = [0,0,0,0,0,0];

prg = frlibPrg(A,[],c,K);
prgRed = prg.ReduceDual('d');

if prgRed.K.l == prg.K.l-1 && prgRed.K.q == 0 && prgRed.K.s(1) == 0
    pass(end+1) = 1;
else
    pass(end+1) = 0;
end

%%%
clear K;
K.f = 1;
K.l = 1;
K.q = 3;
K.s = 1;
A = [[0,0,1,1,0,0];...
     [0,1,0,0,1,-1]...
     ];

c = [0,0,0,0,0,0];

prg = frlibPrg(A,[],c,K);
prgRed = prg.ReduceDual('d');

if prgRed.K.l == prg.K.l-1 && prgRed.K.q == 3 && prgRed.K.s(1) == 0
    pass(end+1) = 1;
else
    pass(end+1) = 0;
end

%%%%%%%%%%
clear K;

K.q = 4;
A = [[1,1,0,0];...
     [0,0,0,0]...
     ];

c = [0,0,0,0];

prg = frlibPrg(A,[],c,K);
prgRed = prg.ReduceDual('dd');

if prgRed.K.l == 1 && prgRed.K.q == 0 
    pass(end+1) = 1;
else
    pass(end+1) = 0;
end

%%%%%%%%%%
clear K;
K.q = 4;
K.s = 2;
A = [[1,1,0,0,0,0,0,0];...
     [0,0,1,0,1,0,0,0]...
     ];

c = [0,0,0,0,0,0,0,0];

prg = frlibPrg(A,[],c,K);
prgRed = prg.ReduceDual('dd');

if prgRed.K.l == 1 && prgRed.K.q == 0 && prgRed.K.s == 0 
    pass(end+1) = 1;
else
    pass(end+1) = 0;
end

break;




%load hybridLyap.mat;
%A = A';
%prg = frlibPrg(A,b,c,K);

%diagonal
% display('Checking diagonal fr')
% prgD = prg.ReducePrimal('d');
% [x,y] = prgD.Solve();
% x0 = prgD.RecoverPrimal(x);
% pass  = prg.CheckPrimal(x0);
% testPass(end+1) = pass & all( prgD.K.s == [6 56 11 1 1 0 11 1 1 0 11 11]);
% if ~(testPass(end))
%     warning('Test case failed')
% end
% display('Checking diagonally dominant fr')
% %diagonally dominant
% prgDD = prgD.ReducePrimal('dd');
% [x,y] = prgDD.Solve();
% x0 = prgDD.RecoverPrimal(x);
% x1 = prgD.RecoverPrimal(x0);
% pass  = prg.CheckPrimal(x1);
% testPass(end+1) = pass;
% 
% if ~(testPass(end))
%     warning('Test case failed')
% end


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

