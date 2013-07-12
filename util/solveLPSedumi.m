function [x,infeas,numErr] = solveLPSedumi(c,Aineq,bineq,Aeq,beq,lbnd,ubnd)

    Meq = size(Aeq,1);
    N = size(Aeq,2);
    
    AineqBnd = speye(N);
    AineqUBnd = AineqBnd( ~isinf(ubnd), :);
    bineqUBnd = ubnd(~isinf(ubnd));

    AineqBnd = speye(N);
    AineqLBnd = -AineqBnd( ~isinf(lbnd), :);
    bineqLBnd = -lbnd(~isinf(lbnd));

    AineqEq = [Aineq;AineqUBnd;AineqLBnd];
    bineqEq  = [bineq;bineqUBnd(:);bineqLBnd(:)];

    numSlack = size(AineqEq,1);

    K.l = numSlack;
    K.f = size(Aeq,2);

    A = [Aeq,sparse(Meq,numSlack); ...
         AineqEq, speye(numSlack) ];

    b = [beq;bineqEq];
    c = [c;sparse(numSlack,1)];
    
    pars.fid = 0;
    [xopt,~,info] = sedumi(A,b,[],K,pars);
    x = sparse( xopt(1:N,1));
    
    infeas = info.pinf == 1;
    numErr = info.numerr > 1;

