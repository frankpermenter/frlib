function pass = runHorn

    pass = [];
    for i=2:5
        
        load(['horn',num2str(i),'.mat']);
        p = frlibPrg(A,b,[],K);
        pred = p.ReducePrimal('dd');

        load(['hornD',num2str(i),'.mat']);
        d = frlibPrg(A,[],c,K);
        dred = d.ReduceDual('dd');
        pass(end+1) = all(dred.K.s == pred.K.s);
        
    end
    
    
    


