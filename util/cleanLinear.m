function [A,b,T]= cleanLinear(A,b,useQR)
  
    if (b == 0)
        b = sparse(size(A,1),1);
    end
   
    if (exist('useQR','var'))
        R = qr(sparse([A,b]'));
        [r,c] = find(R);
        %the first non-zero entry on a row implies
        %the column (e.g. equation) is linearly independent
        [~,indx] = unique(r,'first');
        eqKeep = c(indx);
    else
        eqKeep = find(any(A,2));
    end
    
  
    %a map that maps dual variables for original system to updated
    %Dual variable for equation we've removed gets mapped to 0
    %Dual variable for equation we keep gets mapped to itself
    numEqRed = length(eqKeep);
    numEqOrig = size(A,1);
    T = sparse(eqKeep,1:numEqRed,1,numEqOrig, numEqRed);

    A = A(eqKeep,:);
    b = b(eqKeep,:);
