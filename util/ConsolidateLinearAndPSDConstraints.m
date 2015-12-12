function [Apsd,Kpsd] = ConsolidateLinearAndPSDConstraints(A,K)

    K = coneBase(K); K = K.K;
    
    if (length(K.s)+K.l <= 1) && (K.f == 0)
       Apsd = A; Kpsd = K;
       return;
    end

    %Throw away everything but lin/psd vars
    sLin = 1 + K.f;
    eLin = sLin + K.l-1;
    sPsd = K.f + K.l + sum(K.r) + sum(K.q) + 1;
    A = [A(:,sLin:eLin),A(:,sPsd:end)];

    NumVars = K.l+sum(K.s.^2);
    linearVars = K.l;
    K.s = [ones(1,K.l),K.s];
    K.r = 0; K.q = 0; K.l = 0;
       
    if (NumVars ~= size(A,2))
       error('Num var mismatch') 
    end

   
    indx = sparse(1:linearVars,1:linearVars,ones(linearVars,1),sum(K.s),sum(K.s));
    
    if isempty(linearVars)
        offset = 0;
    else
        offset = linearVars;
    end

    for i=(offset+1):length(K.s)   
       s = sum(K.s(1:i-1))+1;
       e = K.s(i)+s-1;
       indxTemp = ones(K.s(i));     
       indx(s:e,s:e) = indxTemp;
    end
    indx = find(indx(:));
   
    Kpsd.s = sum(K.s);
    Apsd(:,indx) = A;

   