function [A,b,T]= cleanLinear(A,b)
    
    R = qr(sparse([A,b]'));
    [r,c] = find(R);
    %the first non-zero entry on a row implies
    %the column (e.g. equation) is linearly independent
    [~,indx] = unique(r,'first');	
    eqKeep = c(indx);
	
    %must be non-zero and must be first col
    eqRmv = setdiff(1:size(A,1),eqKeep);
	
    T = speye(size(A,1));
	
    T(eqRmv,:) = 0;
    T = T(:,eqKeep);

    A = A(eqKeep,:);
    b = b(eqKeep,:);
	
	

