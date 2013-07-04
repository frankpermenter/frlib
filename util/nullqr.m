function [ns,U] = nullqr(A)

n = size(A,1);
[q,r,e] = qr(A');
p = find( abs(diag(r)) < 10^-4);

if ~isempty(p) 
    indxNull = min(p):size(q,2);
    indxOrth = setdiff(1:size(q,2),indxNull);
else
    indxNull = n+1:size(q,2);
    indxOrth = setdiff(1:size(q,2),indxNull);
end

ns = q(:,indxNull);
[r,c]= find(e(:,indxOrth));
U  = A(:,r);
