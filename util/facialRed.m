classdef facialRed
 
    methods(Static)

        function [x,numErr,infeas] = solveLP(c,Aineq,bineq,Aeq,beq,lbnd,ubnd)
        
            gurobiExists = ~isempty(which('gurobi'));
           % glpkExists = ~isempty(which('glpk'));
            sedumiExists = ~isempty(which('sedumi'));
            
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

            numGens = size(W,2);
            numLam = size(Feq,1);   
            numEq = size(A,1)+1;

            Aeq1 = [A*W,-Feq'];
            Aeq2 = [c(:)'*W,-geq'];
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

            numGens = size(W,1);
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
                [r,c] = find(Aeq);
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
        
        function [success,A,c,K,T] = PolyhedralPrimIter(Z,W)
    
            numGens = size(W,1); 
            
            [cost,Aineq,bineq,Aeq,beq,ubnd,lbnd] = facialRed.BuildPrimLP(Z,W);
            
            W = W';
            [x,numerr,infeas] = facialRed.solveLP(cost,Aineq,bineq,Aeq,beq,lbnd,ubnd);
            success = numerr == 0 & infeas == 0 & -cost'*x > .4;

            if success
                xtemp = x(1:numGens,1);
                xtemp = sparse( xtemp > max(xtemp)*.0001);
                spanS = W*xtemp; 
                [A,c,K,T] = facialRed.ReducePrimal(Z,spanS);
            else
                K = Z.K;
                A = Z.A;
                b = Z.b;
                c = Z.c;
                T =[];
            end
      
        end
        
        function [success,A,c,K,T] = DiagDomPrimIter(A,b,c,K)
           
            Z = coneHelp(A,b,c,K);
            W = Z.extRaysDD();           
            
            [success,A,c,K,T] = facialRed.PolyhedralPrimIter(Z,W);

        end
        
        function [success,A,c,K,T] = DiagPrimIter(A,b,c,K)
            
            Z = coneHelp(A,b,c,K);
            W = Z.extRaysD();           
          
            [success,A,c,K,T] = facialRed.PolyhedralPrimIter(Z,W);

        end

        function [success,A,c,K,Deq,feq] = PolyhedralDualIter(Z,Deq,feq,W)

            W = W';
            numGens = size(W,2);
            
            [cost,Aineq,bineq,Aeq,beq,ubnd,lbnd] = facialRed.BuildDualLP(Z.A,Z.c,Deq,feq,W);
            [x,numerr,infeas] = facialRed.solveLP(cost,Aineq,bineq,Aeq,beq,lbnd,ubnd);
            success = numerr == 0 & infeas == 0 & -cost'*x >= .2;
    
            if success
                xtemp = x(1:numGens,1);       
                S = W*xtemp;
                xtemp = xtemp > max(xtemp)*.0001;
                spanS = W*xtemp;
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
        
        function [success,A,c,K,Deq,feq] = DiagDomDualIter(A,c,K,Deq,feq)

            Z = coneHelp(A,[],c,K);  
            W = Z.extRaysDD(); 
            [success,A,c,K,Deq,feq] = facialRed.PolyhedralDualIter(Z,Deq,feq,W);
            
        end
        
        function [success,A,c,K,Deq,feq] = DiagDualIter(A,c,K,Deq,feq)
            
            Z = coneHelp(A,[],c,K); 
            W = Z.extRaysD(); 
            [success,A,c,K,Deq,feq] = facialRed.PolyhedralDualIter(Z,Deq,feq,W);
            
        end
              
        function [A,c,K,Deq,feq] = ReduceDual(Z,SblkDiag)
                
                A = Z.A(:,1:Z.GetIndx('s',1)-1); 
                c = Z.c(1:Z.GetIndx('s',1)-1); 
                K = Z.K;
                Deq = []; feq = [];
                
                for i = 1:length(K.s) 
                    
                    offset = Z.GetIndx('s',1)-1;
                    [s,e] = Z.GetIndx('s',i);
                    S = mat(SblkDiag(s-offset:e-offset));
               
                    if ( any( eigs(S) <= -10^-12))
                        error('S must be PSD');
                    end

                    if  any(any(abs(S-S') > 10^-12))
                        error('S must be symmetric');
                    end

                    V = S;
                    U = nullqr(V);
             
                    if (K.s(i) > 0)
                        K.s(i) = size(U,2);
                        A = [A,Z.ConjA(U,i)];
                        c = [c,Z.ConjC(U,i)];
                        for j=1:size(V,2)
                            Deq = [Deq;Z.AtimesV(V(:,j),i)];
                            feq = [feq;Z.CtimesV(V(:,j),i)];
                        end
                    end
                end
        end


        function [A,c,K,T] = ReducePrimal(Z,SblkDiag)

                
                A = Z.A(:,1:Z.GetIndx('s',1)-1); 
                c = Z.c(1:Z.GetIndx('s',1)-1); 
                K = Z.K;
                T = [];
 
                for i = 1:length(K.s) 
                 
                    [s,e] = Z.GetIndx('s',i);
                    S = mat( SblkDiag(s:e));
                  
                   % if ( any( eigs(S) <= -10^-12))
                   %     error('S must be PSD');
                   % end

                   % if any(any(abs(S-S')) > 10^-12)
                   %     error('S must be symmetric');
                   % end
                    
                    U = nullqr(S);
                    if (K.s(i) > 0)
                        K.s(i) = size(U,2);
                        A = [A,Z.ConjA(U,i)];
                        c = [c(:)',Z.ConjC(U,i)];
                    end
                    T{i} = U;
                end
                
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




