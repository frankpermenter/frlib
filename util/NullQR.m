function [nullA,rangeAt,rangeAtOrth] = NullQR(A,eps)

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
        startNull = min(p);
        indxNull = startNull:size(q,2);
        indxRange = 1:startNull-1;
    else
        indxNull = n+1:size(q,2);
        indxRange = 1:n;
    end

    nullA = q(:,indxNull);
    Atperm = A(e,:)';
    rangeAt = Atperm(:,indxRange);
    rangeAtOrth = q(:,indxRange);




