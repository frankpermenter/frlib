function frlibTests
    testPass = [];
    opts.useQR = 1;
    opts.verbose = 0;
    
    load hybridLyap.mat;
    A = A';
    prg = frlibPrg(A,b,c,K);

    %diagonal
    TestDisplay('Checking reduction of primal (diagonal)');
    prgD = prg.ReducePrimal('d',opts);
    PrintStats(prgD);
    
    [x,y] = prgD.Solve(opts);   
    x0 = prgD.Recover(x,y);
    pass  = prg.CheckPrimal(x0,10^-4);
    testPass(end+1) = pass & all(prgD.K.s == [6 56 11 1 1 0 11 1 1 0 11 11]);
    if ~(testPass(end))
        error('Test case failed')
    else
       TestDisplay('Pass') 
    end

    %diagonally dominant
    TestDisplay('Checking reduction of primal (diagonally dominant)');
    prgDD = prgD.ReducePrimal('dd',opts);
    PrintStats(prgDD);
    
    [x,y] = prgDD.Solve(opts);
    [x1,y1] = prgDD.Recover(x,y);
    pass  = prgD.CheckPrimal(x1,10^-4);
    testPass(end+1) = pass;

    x2 = prgD.Recover(x1,y1);
    pass = prg.CheckPrimal(x2,10^-4);
    testPass(end+1) = pass;

    if ~(testPass(end))
        error('Test case failed')
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    TestDisplay('Checking reduction of dual');
    load testDual.mat;
    prg = frlibPrg(A,[],c,K);
    prgDD = prg.ReduceDual('dd');
    PrintStats(prgDD);
    
    [~,y] = prgDD.Solve(opts);

    testPass(end+1)  = prg.CheckDual(y,10^-4) & prgDD.K.s == 2;

    if ~(testPass(end))
        error('Test case failed')
    end


    prgD = prg.ReduceDual('d');
    PrintStats(prgD);
    [~,y] = prgD.Solve(opts);
    
    testPass(end+1)  = prg.CheckDual(y,10^-4) & prgDD.K.s == 2;

    if ~(testPass(end))
        error('Test case failed')
    end


    pass = runHorn(opts);
    testPass = [testPass,pass];


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('\n\n%d of %d tests passed \n',sum(testPass),length(testPass))

end


function pass = runHorn(opts)
    TestDisplay('Checking primal and dual reductions (horn form)');
    pass = [];
    for i=2:5
        
        load(['horn',num2str(i),'.mat']);
        p = frlibPrg(A,b,[],K);
        pred = p.ReducePrimal('dd',opts);
        PrintStats(pred);
        
        [x,y] = pred.Solve(opts);
        [x,y] = pred.Recover(x,y);
  
        load(['hornD',num2str(i),'.mat']);
        d = frlibPrg(A,[],c(:),K);
        dred = d.ReduceDual('dd',opts);
        PrintStats(dred);
        [xr,yr] = dred.Solve(opts);
        [~,y] = dred.Recover(xr,yr);

        eps = 10^-4;
        pass(end+1) = all(dred.K.s == pred.K.s) & d.CheckDual(y,eps) ...
            & p.CheckPrimal(x,eps) && all(dred.K.s < d.K.s);
      
    end
    
end

function PrintStats(prg)
    
    prg.PrintStats();

end

function TestDisplay(msg)
   fprintf(['\n']);
   display('**************************************************************');
   fprintf([msg,'\n']);
   display('**************************************************************');
end


