function [A,b]= clean_linear(A,b)

[~,R,P]=qr([A,b]',0);
eqKeep = P( find(abs(diag(R)) > eps));

A = A(eqKeep,:);
b = b(eqKeep,:);


