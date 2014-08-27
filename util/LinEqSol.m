function [x0,nullS]=LinEqSol(A,b)

    eps = 10^-9;
    nullS = NullQR(A);
    if issparse(A)

        %Suppress warnings from \ operator
        s1 = warning('off', 'MATLAB:rankDeficientMatrix');
        s2 = warning('off', 'MATLAB:singularMatrix');

        R = qr(A);
        x = R\(R'\(A'*b));
        %iterative refinement
        r = b - A*x;
        e = R\(R'\(A'*r));
        x0 = x + e;

        use_pinv = norm(A*x0-b) > eps | any(isnan(x0));

        % Restore warnings 
        s1.state = 'on';
        s2.state = 'on';
        warning(s1);
        warning(s2);

    else

        %linsolve is faster, but doesn't accept sparse inputs
        [x0,rank] = linsolve(A,b); %get rank to suppress warnings
        use_pinv = norm(A*x0-b) > eps | any(isnan(x0));

    end

    %switch to svd method 
    if (use_pinv)
        if issparse(A)
            x0 = pinv(full(A))*b; 
        else
            x0 = pinv(A)*b; 
        end
    end


