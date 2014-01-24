
function T =  UtoT(U)
    
    [N,M]=size(U);

    %linear map of q(:) to (U*q)(:)
    T1 = U;
    for i=2:M
        T1=blkdiag(T1,U); 
    end

    %linear map of (U*q)(:) to (U*q*U')(:)
    %for symmetric q
    T2 = sparse([]);
    
    for i=1:N
        cindx = sub2ind([N,M],i*ones(1,M),1:M);
        rindx = [1:N]+(i-1)*N;
        T2(rindx, cindx) = U;  
    end
    
    T = T2*T1;

end
