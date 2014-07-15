function [nullA,rangeAt] = NullQR(A,eps)

    if ~exist('eps','var')
        eps = 10^-8;
    end

    n = size(A,1);
    [q,r,e] = qr(A');
    [e,~] = find(e);

    %find vectors not in span of A'
    if size(r,2) == 1 || size(r,1) == 1
        p = find(abs(r(1,1)) < eps);
    else
        p = find(abs(diag(r)) < eps);
    end
    
    if ~isempty(p)
        indxNull = min(p):size(q,2);
    else
        indxNull = n+1:size(q,2);
    end

    nullA = q(:,indxNull);
    Atperm = A(e,:)';
    rangeAt = Atperm(:,setdiff(1:n,indxNull));
   




