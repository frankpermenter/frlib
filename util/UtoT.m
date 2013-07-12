
function T =  UtoT(U)
    
    [N,M]=size(U);

    %linear map of q(:) to (U*q)(:)
    T1 = U;
    for i=2:M
        T1=blkdiag(T1,U); 
    end

    %linear map of (U*q)(:) to (U*q*U')(:)
    %for symmetric q
    Ut = U';
    T2 = sparse([]);
    for i=1:N
        indx = sub2ind([N,M],i*ones(1,M),1:M);
        for j=1:N
   %         T2(j+(i-1)*N,indx) = Ut(:,j);
        end
    end
    
    for i=1:N
        cindx = sub2ind([N,M],i*ones(1,M),1:M);
        rindx = [1:N]+(i-1)*N;
        T2(rindx, cindx) = Ut';  
    end

   
    
    
    T = T2*T1;
    %Q1 = U*q*U';
    %Q2 = mat(T2*T1*q(:));

end
