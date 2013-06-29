function [A,b,T]= clean_linear(A,b)
    
    [~,R,P]=qr([A,b]',0);
    eqKeep = P( find(abs(diag(R)) > eps));
    eqRmv = setdiff(1:size(A,1),eqKeep);
    eqKeep = sort(eqKeep);

    T = speye(size(A,1));
    T(eqRmv,:) = 0;
    T = T(:,eqKeep);

    A = A(eqKeep,:);
    b = b(eqKeep,:);


