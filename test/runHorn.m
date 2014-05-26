function pass = runHorn

    pass = [];
    for i=2:5
        
        load(['horn',num2str(i),'.mat']);
        p = frlibPrg(A,b,[],K);
        pred = p.ReducePrimal('dd');
        opts.useQR = 1;
        [x,~] = pred.Solve(opts);
        x = pred.RecoverPrimal(x);
  
        load(['hornD',num2str(i),'.mat']);
        d = frlibPrg(A,[],c(:),K);
        dred = d.ReduceDual('dd');
        [~,y] = dred.Solve();

        eps = 10^-4;
        pass(end+1) = all(dred.K.s == pred.K.s) & d.CheckDual(y,eps) ...
            & p.CheckPrimal(x,eps) && all(dred.K.s < d.K.s);
      
    end
    
    
    
    
    
    


