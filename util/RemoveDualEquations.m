function [A,b,c,K,T,y0] = RemoveDualEquations(A,b,c,K)

    Z = ConeBase(K);
    [s,e] = Z.GetIndx('f',1);
    
    if (isempty(s))
       T = eye(size(A,1)); 
       y0 = 0;
       return
    end
    
    Deq = A(:,s:e)';
    feq =  c(s:e);
     
    [y0,T] = lsol(Deq,feq(:));
    
    Anew = T'*A;
    Cnew = c(:) - A'*y0;
    
    A = Anew(:,e+1:end);
    b = T'*b;
    c = Cnew(e+1:end);
    K.f = [];
    
