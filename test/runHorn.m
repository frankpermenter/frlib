%function pass = runHorn
    pass = [];
    for i=2:5
        profile on;
        load(['horn',num2str(i),'.mat']);
        p = frlibPrg(A,b,[],K);
        pred = p.ReducePrimal('dd');


        display('unreduced')
        p.GetSubSpaceDim()
        p.K.s

        display('reduced')
        pred.GetSubSpaceDim()
        pred.K.s

        load(['hornD',num2str(i),'.mat']);
        d = frlibPrg(A,[],c,K);
        dred = d.ReduceDual('dmax');

        dred2 = dred.ReduceDual('dd');
        pass(end+1) = all(dred.K.s == pred.K.s);

        profile off;
        
    end


