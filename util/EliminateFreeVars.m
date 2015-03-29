function [A,b,K] = EliminateFreeVars(A,b,K)

    t = coneBase(K);
    Kfend = t.Kend.f;
    Af=A(:,1:Kfend);
    Aqrs = A(:,Kfend+1:end);
    if ~isempty(Af)       
        nullSpaceAt =  NullQR(Af');
        A = nullSpaceAt'*Aqrs;
        b =  nullSpaceAt'*b;
    end
    K.f = 0;