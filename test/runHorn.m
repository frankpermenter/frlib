function pass = runHorn(opts)

    eps = 10^-4;
    costeps = 10^-1;
    pass = [];
    
    if ~exist('opts','var')
       opts.verbose = 1; 
    end
    
    pars.verbose = opts.verbose;

    for i=2:5
        
        load(['horn',num2str(i),'.mat']);
        c = [1:K.s]'; c = c*c'; c = c(:);
        
        p = frlibPrg(A,b,c,K);
        xorig = sedumi(A,b,c,K,pars);
        costRef = c'*xorig;

        pred = p.ReducePrimal('dd');
        opts.useQR = 1;
        [xr,yr] = pred.Solve(opts);
        x = pred.Recover(xr,yr);
        
        costRed = c'*x;
        pass(end+1) =  norm(costRed-costRef)/norm(costRef) < costeps;
        if (pass(end) == 0)
            xerror('dfdf')
        end
  
        load(['hornD',num2str(i),'.mat']);
        b = [1:size(A,1)]';
        d = frlibPrg(A,b,c(:),K);
        [~,yorig] = sedumi(A,b,c(:),K,pars);
        costRef = b'*yorig;

        opts.removeDualEq = 0;
        dred = d.ReduceDual('dd',opts);
        [xr,yr] = dred.Solve();
        [~,y] = dred.Recover(xr,yr);

        costRed = b'*y;
        pass(end+1) =  norm(costRed-costRef)/norm(costRef) < costeps;
        if (pass(end) == 0)
            xerror('fail')
        end
    
        pass(end+1) = all(dred.K.s == pred.K.s) & d.CheckDual(y,eps) ...
            & p.CheckPrimal(x,eps) && all(dred.K.s < d.K.s);

        if (pass(end) == 0)
            xerror('fail')
        end

        opts.removeDualEq = 1;
        drednoeq = d.ReduceDual('dd',opts);
        [xr,yr] = drednoeq.Solve();
        [~,y] = drednoeq.Recover(xr,yr);

        costRed = b'*y;
        pass(end+1) =  norm(costRed-costRef)/norm(costRef) < costeps;
        if (pass(end) == 0)
            xerror('fail')
        end


        pass(end+1) = all(drednoeq.K.s == pred.K.s) & d.CheckDual(y,eps) ...
            & p.CheckPrimal(x,eps) && all(dred.K.s < d.K.s);   
        if (pass(end) == 0)
            xerror('fail')
        end



    end
end

function xerror(string)
    error(string)
end
    
    
    
    


