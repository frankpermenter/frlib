clear all;
syms y1 y2 y3 y4
y = [y1,y2,y3,y4];

T=[1,-y1,0,-y3; ...
    -y1,2*y2-1,y3,0;...
    0,y3,2*y1-1,-y2;...
    -y3,0,-y2,1];

[A,c,K] = LMI2Sedumi(T);


prg = frlibPrg(A,[],c,K);
prgR = prg.ReduceDual('dd');

if (prgR.K.s ~= 2)
%    error('Reduced K incorrect size') 
end


K = [];
K.l = 5; 
K.q = [4]; 
 
A = sparse(1,9); 
 
A(1,1) = 1; 
A(1,6) = 1; 
b = 0; 
 
prg = frlibPrg(A,b,[],K); 
prgR = prg.ReducePrimal('d');
 
if (prgR.K.l == 4 && prgR.K.q == 0) 
   display('Test Pass')    
else 
   display('Test Failed')  
end 

% if (1,-c): if norm(c) < 1 then set everything to zero
% if ||c|| is one then f =  (f(y) + b) = \lambda (1,c) for \lambda >=0
% write in terms of (f(y) + b) V = 0 ( f(y) + b)/c_i \ge 0.


% primal: replace with zeros or with new variable \lambda.  define
% linear operator that takes \lambda to \lambda(1,c) to recover sol.