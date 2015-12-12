function [rangeAtperp,rangeAt,rangeAtOrth] = NullQR(A,eps)

    if ~exist('eps','var')
        eps = 10^-8;
    end
    
    try
        [q,r,e] = qr(A','vector');
    catch
        [q,r,e] = qr(A');
        [e,~] = find(e);
    end

    if size(r,2) == 1 || size(r,1) == 1
       r_diag = abs( r(1) ); 
    else
       r_diag = abs( diag(r) );
    end
   
    rangeAtOrth = q(:,r_diag  > eps );
    q(:,r_diag  > eps) = [];
    rangeAtperp = q;
    
    Atperm = A(e,:)';
    rangeAt = Atperm(:,r_diag > eps);
    
   



