function [A,C] = makeSymmetric(A,C,K)
    Z = coneHelp(A,[],C,K);
    for i=1:length(K.s)
        [startPos,endPos] = Z.GetIndx('s',i); 
        for j=1:size(A,1)
            Atemp = mat(A(j,startPos:endPos)*1/2);
            Atemp = Atemp+Atemp'; 
            A(j,startPos:endPos) = Atemp(:)';
        end

        C = C(:)';
        Ctemp = mat(C(1,startPos:endPos)*1/2);
        Ctemp = Ctemp+Ctemp'; 
        C(1,startPos:endPos) = Ctemp(:)';
    end
end
