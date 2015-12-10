function frlibTests
    testPass = [];
    opts.useQR = 1;

    eps = 10^-4;
    directory =  fileparts(which('frlibTests.m'));
    load([directory,'/hybrid/hybridLyap.mat']);
    A = A';
    prg = frlibPrg(A,b,c,K);

    %diagonal
    TestDisplay('Checking reduction of primal (diagonal)');
    
    prgD = prg.ReducePrimal('d',opts);
    PrintStats(prgD);
    

    pass = TestSolution(prg,prgD,eps);
    testPass(end+1) = pass & all(prgD.K.s == [6 56 11 1 1 0 11 1 1 0 11 11]);

    %diagonally dominant
    TestDisplay('Checking reduction of primal (diagonally dominant)');
    prgDD = prgD.ReducePrimal('dd',opts);
    PrintStats(prgDD);
    
    testPass(end+1) = TestSolution(prgD,prgDD,eps);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    TestDisplay('Checking reduction of dual');
    load testDual.mat;
    prg = frlibPrg(A,[],c,K);
    prgDD = prg.ReduceDual('dd');
    PrintStats(prgDD);
    
    testPass(end+1) = TestSolution(prg,prgDD,eps);
    testPass(end+1) = prgDD.K.s == 2;

    prgD = prg.ReduceDual('d');
    PrintStats(prgD);
    testPass(end+1) = TestSolution(prg,prgD,eps);
    testPass(end+1)  =  prgD.K.s == 2;


    pass = runHorn(opts);
    testPass = [testPass,pass];


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('\n\n%d of %d tests passed \n',sum(testPass),length(testPass))

end


function pass = runHorn(opts)
    TestDisplay('Checking primal and dual reductions (horn form)');
    pass = [];
    eps = 10^-4;
 
    directory =  fileparts(which('frlibTests.m'));
    for i=2:5
        
        load([directory,'/copos/horn',num2str(i),'.mat']);
        p = frlibPrg(A,b,[],K);
        pred = p.ReducePrimal('dd',opts);
        PrintStats(pred);
        
        pass(end+1) = TestSolution(p,pred,eps);
  
        load([directory,'/copos/hornD',num2str(i),'.mat']);
        d = frlibPrg(A,[],c(:),K);
        dred = d.ReduceDual('dd',opts);
        PrintStats(dred);
        pass(end+1) = TestSolution(d,dred,eps);
        
      
        pass(end+1) = all(dred.K.s == pred.K.s) & (dred.K.s < d.K.s);
      

        
    end
    
end

function pass = TestSolution(prg,prgR,eps)

    if ~isempty(which('sedumi'));
        pars.fid = 0;
        fprintf(' Solving reduced problem...\n');
        [xr,yr] = prgR.Solve(pars);
        [x,y] = prgR.Recover(xr,yr);
        
        if isa(prgR,'reducedPrimalPrg')
            pass = prg.CheckPrimal(x,eps);
            pass = pass & prgR.faces{end}.InLinearSpan(x,eps);
        end

        if isa(prgR,'reducedDualPrg')
            pass = prg.CheckDual(y,eps);
            pass = pass & prgR.faces{end}.InLinearSpan(prg.c-y'*prg.A,eps);
        end
 
    else
        pass = 1; 
    end
    
    pass = pass & prgR.VerifyRedCert(eps);
    
end

function PrintStats(prg)
    
    prg.PrintStats();

end

function TestDisplay(msg)
   fprintf(['\n']);
    display('-------------------------------------------------------------------------')
   fprintf(['\t',msg,'\n']);
   display('-------------------------------------------------------------------------')
end


