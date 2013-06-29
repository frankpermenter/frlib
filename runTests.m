
clear all;
load ./test/hybridLyap.mat;
A = A';
Z = frlibPrg(A,b,c,K);
prg = Z.ReducePrimal('dd');
[x,y] = prg.Solve();
xO = prg.RecoverPrimal(x);
success = prg.CheckPrimal(x);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load bench/files/hinf12.mat;


Z = frlibPrg(A,b,c,K);
prg = Z.ReducePrimal('dd');
[x,y] = prg.Solve();
xO = prg.RecoverPrimal(x);
success = prg.CheckPrimal(x);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%







































clear all;
load bench/files/hinf13.mat;

%load bench/files/minphase.mat;
%A = At';

Z = frlibPrg(A,b,c,K);
prg = Z.ReducePrimal('sdd');
[x,y] = prg.Solve();
xO = prg.RecoverPrimal(x);
success = prg.CheckPrimal(x);
break;

syms y1 y2 y3

M = [y1,0,0;0,-y1,0;0,0,y2];
[A,c,K]=LMI2Sedumi(M);



M = [y1,y1,0;y1,y1,0;0,0,y2];
[A,c,K]=LMI2Sedumi(M);
T = FacialRed(A,[],c,K);
[An,bn,cn,Kn] = T.ReduceDual('dd');
Kn
break;





clear all;
load ./tests/prob_data_2_for_frank.mat;
Z = FacialRed(A,b,c,K);
[Anew,bnew,cnew,Knew,T] = Z.ReducePrimal('d');;

[x,y]=sedumi(Anew,bnew,cnew,Knew);

CheckPrimal(T*x,A,b,K);
break;
%%%%%%%%%%%%%%%%%%%%%
% Lorentz and Linear
%%%%%%%%%%%%%%%%%%%%%
clear A b K c
K.l = 5;
K.q = [4];

A = sparse(1,9);

A(1,1) = 1;
A(1,6) = 1;
b = 0;

[Knew,Anew,T] = PolyFR(K,A,b);

if (Knew.l == 4 && Knew.q == 0)
   display('Test Pass')   
else
   display('Test Failed') 
end

%%%%%%%%%%%%%%%%%%%%%
% Lorentz and Lorentz
%%%%%%%%%%%%%%%%%%%%%
K.l = 0;
K.q = [5,4];
K.s = 1;

A = sparse(1,10);

A(1,1) = 1;
A(1,6) = 1;
b = 0;

[Knew,Anew,T] = PolyFR(K,A,b);
if (Knew.q(1) == 0 && Knew.q(2) == 0)
   display('Test Pass')   
else
   display('Test Failed') 
end

%%%%%%%%%%%%%%%%%%
% Lorentz and PSD  
%%%%%%%%%%%%%%%%%%
K.l = 0;
K.q = [5,4];
K.s = 3;

NumVar = K.l+sum(K.q) + sum(K.s^2);
A = sparse(1,NumVar);

%first two lorentz sum to zero
A(1,1) = 1;
A(1,6) = 1;
b(1,1) = 0;

%diagonal and a variable in lorentz sums to zero
A(2,3) = 1;
A(2,10) = 1;
b(2,1) = 0;

[Knew,Anew,T] = PolyFR(K,A,b);

if (Knew.s == K.s-1)
   display('Test Pass')   
else
   display('Test Failed') 
end



%%%%%%%%%%%%%%%%%%%%%
% Lorentz and Rotated
%%%%%%%%%%%%%%%%%%%%%
K.l = 0;
K.q = [2];
K.r = [4,4];
K.s = 0;

A = sparse(1,10);

A(1,1) = 1;
A(1,7) = 1;
b = 0;
[Knew,Anew,T] = PolyFR(K,A,b);

if (Knew.l == 1 && Knew.r(1) == 4)
   display('Test Pass')   
else
   display('Test Failed') 
end


%%%%%%%%%%%%%%%%%%%%%
% Lorentz and Rotated
%%%%%%%%%%%%%%%%%%%%%
K.l = 0;
K.q = [2];
K.r = [4,4];
K.s = 0;

%set q and and and second r vanish
A = sparse(1,10);
A(1,1) = 1;
A(1,7) = 1;
b(1,1) = 0;

%make first r vanish
A(2,3) = 1;
b(2,1) = 0;

[Knew,Anew,T] = PolyFR(K,A,b);

if (Knew.l == 2 )
   display('Test Pass')   
else
   display('Test Failed') 
end







