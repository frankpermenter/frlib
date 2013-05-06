function [nneg] = sedumiNonNeg(K)

    nneg = [];
    offset = K.f;
    
    if K.l > 0
        nneg = [1:K.l]'+offset;
        offset = nneg(end);
    end
    
    %1st variable of each lorentz cone constraint
    if any(K.q > 0)  
        firstVar = [0,cumsum(K.q)]+1 + offset;
        firstVar = firstVar(K.q > 0);
        nneg = [nneg;firstVar'];
        %ending position of lorentz cone constraints
        offset = offset + sum(K.q);
    end
    
    %1st and 2nd variable of each rotated lorentz cone constraint
    if any(K.r > 0)
        firstVar = [0,cumsum(K.r)]+1+offset;
        firstVar = firstVar(K.r>0);
        secondVar = firstVar+1; 
        nneg = [nneg;firstVar';secondVar'];
        %ending position of socp constraints
        offset = offset + sum(K.r);
    end
    
    for i = 1:length(K.s)
        if K.s(i) == 0
            continue;
        end

        N = K.s(i);
        diagPos = [1:N+1:N*N]';
        nneg = [nneg;diagPos+offset]; 
        %ending position of last SDP constraint is a diagonal element
        offset = nneg(end);
    end

