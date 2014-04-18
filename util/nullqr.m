function [ns,U] = nullqr(A)

    n = size(A,1);
    [q,r,e] = qr(A');
   

    %find vectors not in span of A'
    p = find( abs(diag(r)) < 10^-8);
    if ~isempty(p)
        indxNull = min(p):size(q,2);
    else
        indxNull = n+1:size(q,2);
    end

    ns = q(:,indxNull);




