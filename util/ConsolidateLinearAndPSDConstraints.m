function [Apsd,K] = ConsolidateLinearAndPSDConstraints(A,K)

    if (length(K.s)+K.l <= 1)
        return;
    end

    t = coneBase(K);
    K = t.cleanK(K);
    K.s = [ones(1,K.l),K.s];
    K.l = 0;

    sLin = t.Kstart.l;
    eLin = t.Kend.l;
    sPsd = t.Kstart.s(1);
    A = [A(:,1:t.Kend.f),A(:,eLin+1:sPsd-1),A(:,sLin:eLin),A(:,sPsd:end)];

    t = coneBase(K);

    if (t.NumVar ~= size(A,2))
       error('Num var mismatch') 
    end

    for i=1:size(A,1)
        m = [];
        for j=1:length(t.Kstart.s)   
           m = blkdiag(m,mat(A(i,t.Kstart.s(j):t.Kend.s(j))));
        end
        Apsd(i,:) = m(:)';
    end

    K.s = sum(K.s);