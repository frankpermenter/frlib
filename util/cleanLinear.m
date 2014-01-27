function [A,b,T]= cleanLinear(A,b)
    T = speye(size(A,1));
    
    if (b == 0)
        b = sparse(size(A,1),1);
    end

    R = qr(sparse([A,b]'));
    [r,c] = find(R);

    %the first non-zero entry on a row implies
    %the column (e.g. equation) is linearly independent
    [~,indx] = unique(r,'first');
    eqKeep = c(indx);

    eqRmv = setdiff(1:size(A,1),eqKeep);

    %a map that maps dual variables for original system to updated
    T = speye(size(A,1));
    T(eqRmv,:) = 0;
    T = T(:,eqKeep);

    A = A(eqKeep,:);
    b = b(eqKeep,:);
