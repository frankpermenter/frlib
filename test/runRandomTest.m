function runRandomTest

    pass = [];
    K.s =  3;
    rankS = 1;%floor(2*rand(1))+2;
    numEq = floor(2*rand(1))+2;
    eps = 10^-8;

    [A,b,c,S] = RandomPrimalTestPSDCase(K,numEq,rankS,1);
    fr = frlibPrg(A,b,c,K);
    frRed = fr.ReducePrimal('d');
    pass = [pass;CheckPrimalTest(fr,frRed,S,rankS,eps)];

    if (all(pass) == 1) 
        display('All tests passed!');
    else
        error('Test failed!');
    end

    numEq = 4;
    [A,b,c,S] = RandomDualTestPSDCase(K,numEq,rankS,1);
    fr = frlibPrg(A,b,c,K);
    frRed = fr.ReduceDual('d');
    pass = [pass;CheckDualTest(fr,frRed,S,rankS,eps)];

    if (all(pass) == 1) 
        display('All tests passed!');
    else
        error('Test failed!');
    end



end

function [A,b,c,S] = RandomPrimalTestPSDCase(K,numEq,rankS,numIter)
    
    N = K.s;
    if (N < rankS) 
        error('N less than rank of reducing certificate S')
    end

    A = rand(numEq-1,N*N);
    S = diag([zeros(N-rankS,1);rand(rankS,1)]);
    S = S(:)';

    y = randn(numEq-1,1);
    A(numEq,:) = S - (A'*y)';

    y = [y;1];
    yPerp = null(y');
    b = yPerp(:,1);

    K.s = [N];
    c = rand(N*N,1);

end


function [A,b,c,S] = RandomDualTestPSDCase(K,numEq,rankS,numIter)
    
    N = K.s;
    if (N < rankS) 
        error('N less than rank of reducing certificate S')
    end

    A = rand(numEq-1,N*N);
    S = diag([zeros(N-rankS,1);rand(rankS,1)]);
    
    V = NullQR(S);
    
    Z = coneBase(K);
    Stemp = S(:)';    
    Stemp(Z.indxDiag{1}) = Stemp(Z.indxDiag{1})/2;
    Stemp = Z.UpperTri(Stemp);
    Snull = null(Stemp);

    for (i=1:numEq+1)
        
        Sperp = Snull*randn(size(Snull,2),1);
        temp = zeros(K.s,K.s);
        temp(Z.UpperTriIndx()) = Sperp;
        temp(Z.LowerTriIndx()) = Sperp;

        if (i<numEq+1)
            A(i,:) = temp(:)';
        else
            c = temp(:)';
        end
    end

    b = rand(numEq,1);

end


function pass = CheckDualTest(frOrig,frRed,S,rankS,eps)
    [~,y,info] = frRed.Solve();
    [~,yorig,infoOrig] = frOrig.Solve();

    pass = frRed.K.s+rankS == frOrig.K.s;
    pass = [pass;norm( frOrig.b(:)'*y-frOrig.b(:)'*yorig) < eps];
    
    if (info.pinf == 1 || info.dinf == 1)
        pass = [pass;infoOrig.pinf == info.pinf];
        pass = [pass;infoOrig.dinf == info.dinf];
    else 
        pass = [pass;frOrig.CheckDual(y,eps)];
        pass = [pass;CheckRedCert(mat(frOrig.c(:)-frOrig.A'*yorig),S,eps)];
    end
    
    if (all(pass) ~= 1)
        xerror('fail');
    end
    
    pass = all(pass);
end



function pass = CheckPrimalTest(frOrig,frRed,S,rankS,eps)
    [x,y,info] = frRed.Solve();
    [xOrig,~,infoOrig] = frOrig.Solve();

    pass = frRed.K.s+rankS == frOrig.K.s;

    if (info.pinf == 1 || info.dinf == 1)
        pass = [pass;infoOrig.pinf == info.pinf];
        pass = [pass;infoOrig.dinf == info.dinf];
    else
        display('checking')
        xRecover = frRed.RecoverPrimal(x);
        pass = frOrig.CheckPrimal(xRecover,eps);
        pass = [pass,CheckRedCert(xOrig,S,eps)];
    end
    if (all(pass) ~= 1)
        xerror('fail');
    end
    pass = all(pass);
end

function pass = CheckRedCert(x,S,eps)
    try
    if (norm(trace(x*S)) > eps)
        xerror(['Inner product with reducing certificate not zero: trace(XS) = ',num2str(norm(trace(x*S)))]);
        pass = 0;
    else
        pass = 1;
    end
    catch
        warning('dfd')
        pass = 0;
    end
end

function xerror(str)
    error(str)
end

