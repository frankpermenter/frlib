function [A,b] = EliminateFreeVars(A,b,Kf)

    Af=A(:,1:Kf);
    Aqrs = A(:,Kf+1:end);
    if ~isempty(Af)       
        nullSpaceAt =  NullQR(Af');
        A = nullSpaceAt'*Aqrs;
        b =  nullSpaceAt'*b;
    end
   