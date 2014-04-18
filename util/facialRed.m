classdef facialRed

    methods(Static)

        function [x,numErr,infeas] = solveLP(c,Aineq,bineq,Aeq,beq,lbnd,ubnd)

            gurobiExists = 0; ~isempty(which('gurobi'));
           % glpkExists = ~isempty(which('glpk'));
            sedumiExists = 0; ~isempty(which('sedumi'));

            numErr = 0;
            infeas = 0;
            if (gurobiExists)

                model.A = sparse([Aineq;Aeq]);
                model.obj = full(c);
                model.rhs = full([bineq;beq]);
                model.lb = full(lbnd);
                model.ub = full(ubnd);
                model.sense = char( ['<'*ones(length(bineq),1);'='*ones(length(beq),1)]);

                params.method = 1;
                params.crossover = 0;
                params.outputflag = 0;
                result = gurobi(model,params);
                flag = result.status;

                if (strcmp(flag,'INFEASIBLE'))
                    infeas = 1;
                    x = [];
                    return;
                else
                    x = result.x;
                end

                return
            end

            %if (glpkExists)
               %x = linprog(c,Aineq,bineq,Aeq,beq,lbnd,ubnd);
            %   return
            %end

            if (sedumiExists)
                [x,infeas,numErr] = solveLPSedumi(c,Aineq,bineq,Aeq,beq,lbnd,ubnd);
                return;
            else
                [x,~,flag] = linprog(c,Aineq,bineq,Aeq,beq,lbnd,ubnd);
                infeas = infeas == -2 | infeas == -5;
                numErr = ~infeas && flag ~= 1;
            end

        end

        function [c,Aineq,bineq,Aeq,beq,ubnd,lbnd] = BuildDualLP(A,c,Feq,geq,W)

            numGens = size(W,1);
            numLam = size(Feq,1);
            numEq = size(A,1)+1;

            if numGens == 0
                error('No generators specified. size(W,1) = 0')
            end

            Aeq1 = [A*W',-Feq'];
            Aeq2 = [c(:)'*W',-geq'];
            Aeq = [Aeq1;Aeq2];
            Aeq = [Aeq,sparse(numEq,numGens)];
            beq = sparse(size(Aeq,1),1);

            if size(Aeq,1) > size(Aeq,2)
                [Aeq,beq] = cleanLinear(Aeq,beq);
            else
                [r,c] = find(Aeq);
                rKeep = unique(r);
                Aeq = Aeq(rKeep,:);
                beq = beq(rKeep,1);
            end

            Aineq = [-speye(numGens),sparse(numGens,numLam),speye(numGens)];
            bineq = zeros(size(Aineq,1),1);
            c = [zeros(numGens+numLam,1);-ones(numGens,1)];
            ubnd = [Inf*ones(numGens+numLam,1);ones(numGens,1)];
            lbnd = [-Inf*ones(numGens+numLam,1);zeros(numGens,1)];

        end

        function [c,Aineq,bineq,Aeq,beq,ubnd,lbnd] = BuildPrimLP(Z,W)

            %equation: \sum_i A_i mu_i = \sum_i s_i W_i
            % s_i \ge t_i \ge 0
            %cols: (mu,y,t)


            numGens = size(W,1);
            if numGens == 0
                error('No generators specified. size(W,1) = 0')
            end

            numMu = size(Z.A,1);


            A = Z.lowerTri(Z.A);
            W = Z.lowerTri(W);

            Aeq1 = [-W',A'];
            Aeq2 = [sparse(1,numGens),Z.b'];

            Aeq = [Aeq1;Aeq2];
            Aeq = [Aeq,sparse(size(Aeq,1),numGens)];
            beq = sparse(size(Aeq,1),1);

            %Remove dependent equations.
            if size(Aeq,1) > size(Aeq,2)
                [Aeq,beq] = cleanLinear(Aeq,beq);
            else
                [r,~] = find(Aeq);
                rKeep = unique(r);
                Aeq = Aeq(rKeep,:);
                beq = beq(rKeep,1);
            end

            Aineq = [-speye(numGens),sparse(numGens,numMu),speye(numGens)];
            bineq = sparse(size(Aineq,1),1);
            c = [sparse(numGens+numMu,1);-ones(numGens,1)];
            ubnd = [Inf*ones(numGens+numMu,1);ones(numGens,1)];
            lbnd = [-Inf*ones(numGens+numMu,1);sparse(numGens,1)];

        end

        function [success,A,c,K,T,S] = PolyhedralPrimIter(Z,W)

            numGens = size(W,1);
            K = Z.K;
            A = Z.A;
            b = Z.b;
            c = Z.c;
            T =[];
            S = [];
            success = 0;

            if (numGens == 0)
                return;
            end

            [cost,Aineq,bineq,Aeq,beq,ubnd,lbnd] = facialRed.BuildPrimLP(Z,W);

            [x,numerr,infeas] = facialRed.solveLP(cost,Aineq,bineq,Aeq,beq,lbnd,ubnd);
            success = numerr == 0 & infeas == 0 & -cost'*x > .1;

            if success
                xtemp = x(1:numGens,1);
                xtemp = sparse( xtemp > max(xtemp)*.01);
                S = x(numGens+1:numGens+size(Z.A,1));
                spanS = W'*xtemp;
                [A,c,K,T] = facialRed.ReducePrimal(Z,spanS);
            else
                K = Z.K;
                A = Z.A;
                b = Z.b;
                c = Z.c;
                T =[];
            end

        end

        function [success,A,c,K,T,S] = DiagDomPrimIter(A,b,c,K)

            Z = coneHelp(A,b,c,K);
            W = Z.extRaysDD();
            
            [success,A,c,K,T,S] = facialRed.PolyhedralPrimIter(Z,W);

        end

        function [success,A,c,K,T,S] = DiagPrimIter(A,b,c,K)

            Z = coneHelp(A,b,c,K);
            W = Z.extRaysD();
            [success,A,c,K,T,S] = facialRed.PolyhedralPrimIter(Z,W);

        end

        function [success,A,c,K,Deq,feq,S] = PolyhedralDualIter(Z,Deq,feq,W)

            numGens = size(W,1);
            K = Z.K;
            A = Z.A;
            b = Z.b;
            c = Z.c;
            success = 0;
            S = [];

            if (numGens == 0)
                return;
            end

            [cost,Aineq,bineq,Aeq,beq,ubnd,lbnd] = facialRed.BuildDualLP(Z.A,Z.c,Deq,feq,W);
            [x,numerr,infeas] = facialRed.solveLP(cost,Aineq,bineq,Aeq,beq,lbnd,ubnd);
            success = numerr == 0 & infeas == 0 & -cost'*x >= .2;

            if success
                xtemp = sparse(x(1:numGens,1));
                S = W'*xtemp;
                xtemp = xtemp > max(xtemp)*.0001;
                spanS = W'*xtemp;
                [A,c,K,DeqN,feqN] = facialRed.ReduceDual(Z,spanS);
                Deq = [Deq;DeqN];
                feq = [feq;feqN];
            else
                K = Z.K;
                A = Z.A;
                b = Z.b;
                c = Z.c;
            end

        end

        function [success,A,c,K,Deq,feq,S] = DiagDomDualIter(A,c,K,Deq,feq)

            Z = coneHelp(A,[],c,K);
            S = [];
            if (Z.anyConicVars() == 0)
                success = 0;
            else
                W = Z.extRaysDD();
                [success,A,c,K,Deq,feq,S] = facialRed.PolyhedralDualIter(Z,Deq,feq,W);
            end

        end

        function [success,A,c,K,Deq,feq,S] = DiagDualIter(A,c,K,Deq,feq)

            Z = coneHelp(A,[],c,K);
            S = [];
            if (Z.anyConicVars() == 0)
                success = 0;
            else
                W = Z.extRaysD();
                [success,A,c,K,Deq,feq,S] = facialRed.PolyhedralDualIter(Z,Deq,feq,W);
            end

        end

        function [A,c,K,Deq,feq] = ReduceDual(Z,SblkDiag)

            Deq = []; feq = [];
            K = Z.K;

            [s,e] = Z.GetIndx('f',1);
            A = Z.A(:,s:e);
            c = Z.c(s:e)';

            [s,e] = Z.GetIndx('l',1);
            S = SblkDiag(s:e);
            
            if (any(S))
                
                colsRemove = find(S)+s-1;
                colsKeep = setdiff(s:e,colsRemove);
                A = [A,Z.A(:,colsKeep)];
                c = [c,Z.c(colsKeep)'];
                Deq = Z.A(:,colsRemove)';
                feq = Z.c(colsRemove);
                K.l = K.l-length(colsRemove);
                
            else
                
                A = [A,Z.A(:,s:e)];
                c = [c,Z.c(s:e)'];
                
            end

            newLinearA=[];
            newLinearC=[];

            for i=1:length(Z.K.q)
                
                [s,e] = Z.GetIndx('q',i);
                S = SblkDiag(s:e);
                
                if (any(S))
                    [DeqTemp,feqTemp,ineqA,ineqC] = facialRed.ReduceLorentzDual(S,Z.A(:,s:e),Z.c(s:e));
                    Deq = [Deq;DeqTemp];
                    feq = [feq;feqTemp];
                    newLinearA = [newLinearA,ineqA];
                    newLinearC = [newLinearC,ineqC];
                    K.q(i) = 0;
                else
                    A = [A,Z.A(:,s:e)];
                    c = [c,Z.c(s:e)'];
                end
                
            end

            for i=1:length(Z.K.r)
                
                [s,e] = Z.GetIndx('r',i);
                S = SblkDiag(s:e);
                
                if (any(S) )
                    [DeqTemp,feqTemp,ineqA,ineqC] = facialRed.ReduceLorentzDual(S,Atemp,c);
                    Deq = [Deq;DeqTemp];
                    feq = [feq;feqTemp];
                    newLinearA = [newLinearA,ineqA];
                    newLinearC = [newLinearC,ineqC];
                    K.r(i) = 0;
                else
                    A = [A,Z.A(:,s:e)];
                    c = [c,Z.c(s:e)'];
                end
                
            end

            As = []; cs = [];
            
            indicesFlag = zeros(size(Z.A,2),1);
            
            for i = 1:length(K.s)

                [s,e] = Z.GetIndx('s',i);
                S = mat(SblkDiag(s:e));
                Aconj = Z.A(:,s:e);
                Cconj = Z.c(s:e);

                if (nnz(diag(S)) == nnz(S))
                   diagonalS = true; 
                else
                   diagonalS = false; 
                end
  
                if (diagonalS)
                    
                    rowColFlag = zeros(K.s(i),K.s(i));
                    rowColKeep = find(diag(S)==0);
                    rowColFlag(rowColKeep,rowColKeep) = 1;
                     
                    Aconj = Aconj(:,rowColFlag(:) > 0);
                    Cconj = Cconj(rowColFlag(:) > 0 )';
                    
                    rowColFlag = zeros(K.s(i),K.s(i));
                    rowColFlag(:,diag(S)~=0) = 1;
                    indicesFlag(find(rowColFlag(:)) + s-1) = 1;
                 
                    K.s(i) = length(rowColKeep);  
                    
                else
                    
                    V = S(:,any(S~=0,2));
                    U = nullqr(S);
                   
                    [Aconj,Cconj] = Z.ConjAandC(U,i);
                    
                    for j=1:size(V,2)
                        Deq = [Deq;Z.AtimesV(V(:,j),i)];
                        feq = [feq;Z.CtimesV(V(:,j),i)];
                    end
                    
                    K.s(i) = size(U,2);
                    
                end
                
                As = [As,Aconj];
                cs = [cs,Cconj];  
            end
  
            if any(indicesFlag)
                Deq = [Deq;Z.A(:,indicesFlag > 0)'];
                feq = [feq;Z.c(indicesFlag > 0)];
            end


            A = [A, As];
            c = [c(:)', cs(:)'];

        end


        function [A,c,K,T] = ReducePrimal(Z,SblkDiag)

            A = Z.flqrCols(Z.A);
            c = Z.flqrCols(Z.c(:)');

            K = Z.K;
            T = [];

            As = [];
            cs = [];

            for i = 1:length(K.s)

                [s,e] = Z.GetIndx('s',i);
                S = mat(SblkDiag(s:e));

                if ( any( eig(S) <= -10^-12))
                    error('S must be PSD');
                end

                if any(any(abs(S-S')) > 10^-12)
                    error('S must be symmetric');
                end

                U = nullqr(S);
                if (K.s(i) > 0)
                    K.s(i) = size(U,2);
                    As = [As,Z.ConjA(U,i)];
                    cs = [cs(:)',Z.ConjC(U,i)];
                end
                T{i} = U;
                
            end

            A = [A, As];
            c = [c, cs];

        end


        function [success,A,c,K,T] = SDDPrimIter(A,b,c,K)

            T = [];
            Z = coneHelp(A,b,c,K);
            prg = spotsosprg;
            m = size(Z.A,1);
            tr = 0;

            [prg,u] = prg.newFree(m);
            prg = prg.withEqs(u'*Z.b);

            for i=1:length(K.s)
                if (K.s(i) > 1 )
                    [prg,Stemp] = prg.newSDD(K.s(i));
                    S2 = mat(Z.AdjA(u,i));
                    S{i} = S2;
                    prg = prg.withEqs(Stemp-S2);
                    tr = tr + trace(Stemp);
                end
            end

            prg = prg.withEqs(tr-1);

            try
                sol = prg.optimize();
                success = sol.info.pinf == 0 & sol.info.numerr == 0;
            catch
                success = 0;
            end

            S =[];
            if success
                for i=1:length(S)
                    Sn = sol.eval(S{i});
                    S = [S,Sn(:)']
                end
                [A,c,K,T] = facialRed.ReducePrimal(Z,S);
            end

        end


        function [Deq,feq,A,c] = ReduceLorentzDual(S,A,c)

            S = S(:);
            c = c(:);
            S = S/S(1);

            if norm(S(2:end)) < 1
                %c-A'y = 0
                Deq = A';
                feq = c;
                A = []; c = [];
            else

                sPerp = [1;-S(2:end)]
                N = length(c);

                norm_sPerpSqr = sum(sPerp.*sPerp);

                % (c-A'y) must be parallel to sPerp
                % (c-A'y) - P(c-A'y) = 0
               
                Deq = (eye(N) - sPerp/norm_sPerpSqr*sPerp')*A';
                feq = (eye(N) - sPerp/norm_sPerpSqr*sPerp')*c;

                %c(1) \ge A(:,1)'*y
                A = A(:,1);
                c = c(1);

            end
        end

        function [success,A,c,K,Deq,feq] = SDDDualIter(A,c,K,Deq,feq)

            Z = coneHelp(A,[],c,K);
            prg = spotsosprg;

            tr = 0;
            SdotM = 0;
            Cdot = 0;

            for i=1:length(K.s)
                [prg,Stemp] = prg.newSDD(K.s(i));
                S{i} = Stemp;
                SdotM = SdotM + Z.AdotX(Stemp,i);
                SdotC = SdotC + Z.CdotX(Stemp,i);
                tr = tr + trace(Stemp);
            end
            [prg,lam]  = prg.newFree(length(feq));
            prg = prg.withEqs(SdotM+Deq'*lam);
            prg = prg.withEqs(SdotC+Deq'*lam);
            prg = prg.withEqs(tr-1);
            try
                sol = prg.optimize();
                success = sol.info.pinf == 0;
            catch
                success = 0;
            end

            if success
                [A,c,K,Deq,feq] = facialRed.ReduceDual(Z,S);
            end

        end

    end
end
