function [Knew,Anew,T] = DiagFR(K,A,b,useDD)
    
    Knew = K;
    if (size(A,1) > size(A,2))
       Anew = A'; 
    else
       Anew = A; 
    end
    
    T = speye(size(Anew,2));
    
    if ~exist('useDD','var') 
        useDD = 0;
    end
    
    enableDD = 0;
    
    while (1)
        [varRmv,Knew,Anew,Tform]=doIter(Knew,Anew,b,enableDD);
        T = T*Tform;
        if isempty(varRmv)
            
           if (enableDD == 1) || useDD == 0
                return
           else
               enableDD = 1;
           end
        end
    end

end









function [varRmv,Knew,Anew,T] = doIter(K,A,b,DD)

    varRmv=[];Knew=[];Anew=[];T=[];col=[];
    K = cleanK(K);

    NumVars = K.f + K.l + sum(K.q) + sum(K.r) + sum(K.s.^2);
    
    %get variables that are non-negative if conic
    %constraints in K hold
    nneg = sedumiNonNeg(K);
    
    A_ineq = sparse(length(nneg),NumVars);
    for i=1:length(nneg)
        A_ineq(i,nneg(i)) = -1;
    end

    if (DD == 1)
        V = ones(2);

        A_ineq1 = -sedumiNbyN(K,V);
        V(2,1) = -1;V(1,2) = -1;

        A_ineq2 = -sedumiNbyN(K,V);
        A_ineq1 = [A_ineq1;A_ineq2];
        A_ineq = [A_ineq;A_ineq1];
    else
        A_ineq1 = [];
    end


    b_ineq = zeros(size(A_ineq,1),1);

    y = lpact(A_ineq,b_ineq,A,b);

    N = size(A_ineq,1)-size(A_ineq1,1);
        
    if (0)%DD == 0)
        Knew = K;
        save test.mat Knew A_ineq A_ineq1 A_ineq2 y A b
        zeroVar2 = removeDiag(K,A_ineq1,y(N+1:end));
    else
        zeroVar2 = [];
    end

    zeroVar2 = [];
    %active constraints

    activeCnst = find(y(1:N));
    [~,zeroVar1] = find(A_ineq(activeCnst,:));

    indxArray = [zeroVar1;zeroVar2];
    %other vars that must vanish
    [varRmv,Knew] = FindMustVanish(K,[zeroVar1;zeroVar2]);
    %build map: oldVar = T*newvar
    T=speye(NumVars);
    T(varRmv,:) = 0;
    
    indxKeep = setdiff(1:NumVars,varRmv);
    T = T(:,indxKeep);
    %move lorentz/rotated lorentz constraints
    %with one variable to K.l
    convertToLinear = find(Knew.q == 1);
    offsetQ = Knew.f + Knew.l;
    
    T2 = speye(length(indxKeep));
    for i = 1:length(convertToLinear)
        indxPaste = offsetQ + sum(Knew.q(1:convertToLinear(i)-1))+1;
        indxCut = Knew.l + 1;
        T2 = CutAndPaste(T2,indxCut,indxPaste);
        Knew.l = Knew.l + 1;
        Knew.q(convertToLinear(i)) = 0;
    end

    convertToLinear = find(Knew.r == 1);
    offsetR = Knew.f + Knew.l+sum(Knew.q);
    for i = 1:length(convertToLinear)
        indxPaste = offsetR + sum(Knew.r(1:convertToLinear(i)-1))+1;
        indxCut = Knew.l + 1;
        T2 = CutAndPaste(T2,indxCut,indxPaste);
        Knew.l = Knew.l + 1;
        Knew.r(convertToLinear(i)) = 0;
    end
    T = T*T2;
    Anew = A * T;
  
end


function T = CutAndPaste(T,indxCut,indxPaste)


if (indxCut == indxPaste)
    return
end
% cut < paste

if (indxPaste > indxCut)

T = [T(1:indxCut-1,:);T(indxCut+1:indxPaste,:);T(indxCut,:);T(indxPaste+1:end,:)];
else
    T = [T(1:indxPaste,:);T(indxCut,:);T(indxPaste+1:indxCut-1,:);T(indxCut+1:end,:)];
end

end

function offset = CnstToOffset(K,type,num)

    
    switch type

        case 's'
            offset = K.f + K.l + sum(K.q) + sum(K.r) + sum(K.s(1:num-1).^2);

        case 'r'
            offset = K.f + K.l + sum(K.q)+ sum(K.r(1:num-1));

        case 'q'
            offset = K.f + K.l + sum(K.q(1:num-1));

    end


end




function K = cleanK(K)

    if ~isfield(K,'f')
        K.f = 0;
    end

    if ~isfield(K,'l')
        K.l = 0;
    end

    if ~isfield(K,'q')
        K.q = 0;
    end

    if ~isfield(K,'s')
        K.s = 0;
    end
    
    if ~isfield(K,'r')
        K.r = 0;
    end
    
 
end

    

    


