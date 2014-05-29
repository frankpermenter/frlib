function frlibTests
    testPass = [];
    opts.useQR = 1;
    opts.verbose = 0;
    
    load hybridLyap.mat;
    A = A';
    prg = frlibPrg(A,b,c,K);

    %diagonal
    testdisplay('Checking reduction of primal (diagonal)');
    prgD = prg.ReducePrimal('d');
    PrintStats(prgD);
    
    [x,~] = prgD.Solve(opts);   
    x0 = prgD.RecoverPrimal(x);
    pass  = prg.CheckPrimal(x0,10^-4);
    testPass(end+1) = pass & all(prgD.K.s == [6 56 11 1 1 0 11 1 1 0 11 11]);
    if ~(testPass(end))
        error('Test case failed')
    end

    %diagonally dominant
    testdisplay('Checking reduction of primal (diagonally dominant)');
    prgDD = prgD.ReducePrimal('dd');
    PrintStats(prgDD);
    
    [x,y] = prgDD.Solve(opts);
    x1 = prgDD.RecoverPrimal(x);
    x2 = prgD.RecoverPrimal(x1);
    pass  = prg.CheckPrimal(x2,10^-4);
    testPass(end+1) = pass;

    if ~(testPass(end))
        error('Test case failed')
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    testdisplay('Checking reduction of dual');
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
    testdisplay('Checking primal and dual reductions (horn form)');
    pass = [];
    for i=2:5
        
        load(['horn',num2str(i),'.mat']);
        p = frlibPrg(A,b,[],K);
        pred = p.ReducePrimal('dd');
        PrintStats(pred);
        
        [x,~] = pred.Solve(opts);
        x = pred.RecoverPrimal(x);
  
        load(['hornD',num2str(i),'.mat']);
        d = frlibPrg(A,[],c(:),K);
        dred = d.ReduceDual('dd');
        PrintStats(dred);
        [~,y] = dred.Solve(opts);
        
        eps = 10^-4;
        pass(end+1) = all(dred.K.s == pred.K.s) & d.CheckDual(y,eps) ...
            & p.CheckPrimal(x,eps) && all(dred.K.s < d.K.s);
      
    end
    
end

function PrintStats(prg)
    
    prg.PrintStats();

end

function testdisplay(msg)
   fprintf(['\n']);
   display('**************************************************************');
   fprintf([msg,'\n']);
   display('**************************************************************');
end


