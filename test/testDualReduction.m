function testDualReduction
   
    pass = [];
    load testDual.mat; 
    b = [1,0,0,0]';
    
    % syms y1 y2 y3 y4
    % T= [[ y1,  y1,     y3,  0,  0]
    % [ y1, -y1,     y2,  0,  0]
    % [ y3,  y2, y3 + 1,  0,  0]
    % [  0,   0,      0, y4,  0]
    % [  0,   0,      0,  0, y3] ]  ;

    %everything vanishes except 2x2 block
    numEq = nchoosek(5+1,2)-3;
    
    exactlyZero = 1;
    CheckTest(numEq,exactlyZero,A,b,c,K);

    Af = [0,0,0,1]';
    cf = [2];
    A = [Af,A];
    c = [cf,c(:)'];
    K.f = 1;
    exact = 1;
    CheckTest(numEq+1,exact,A,b,c,K);

end


function CheckTest(numEq,exactlyZero,A,b,c,K)

    eps = 10^-8;
    pass = [];
  

    b = [1,0,0,0]';
    A = sparse(A); b = sparse(b); c = sparse(c);
    f = frlibPrg(A,b,c,K);

    opts.removeDualEq = 1;
    rnoeq = f.ReduceDual('dd',opts);

    opts.removeDualEq = 0;
    req = f.ReduceDual('dd',opts); 
    if (req.K.f ~= numEq)
        xerror('test failed');
    else
        pass = [pass,1];
    end

    [xreq,yreq] = sedumi(req.A,req.b,req.c,req.K);
    xreq = sparse(xreq); yreq = sparse(yreq);
    
    [xrnoeq,yrnoeq] = sedumi(rnoeq.A,rnoeq.b,rnoeq.c,rnoeq.K);
    %for baby problem, sedumi doesn't output sparse answers?
    xrnoeq = sparse(xrnoeq); yrnoeq = sparse(yrnoeq);
    
    [xrecov,yrecov] = req.Recover(xreq,yreq);
    f.CheckSolution(xrecov,yrecov,eps);

    [xrecovnoeq,yrecovnoeq] = rnoeq.Recover(xrnoeq,yrnoeq);
    f.CheckSolution(xrecovnoeq,yrecovnoeq,eps);

    if  norm(yrecov(1:end-1)) > eps 
        xerror('test fail')
    end

    if (exactlyZero)
        if ~all(yrecovnoeq(1:end-1) == 0 )
            xerror('test fail')
        end
    else
        if  norm(yrecovnoeq(1:end-1)) > eps 
            xerror('test fail')
        end
    end

end


function xerror(string)
    error(string);
end


