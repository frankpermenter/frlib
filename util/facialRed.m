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

             %   params.method = 1;
              %  params.crossover = 0;
              %  params.outputflag = 0;
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

        function [c,Aineq,bineq,Aeq,beq,ubnd,lbnd,solMap] = BuildDualLP(A,c,U,V,W,cone)

            numGens = size(W,1);
            numEq = size(A,1)+1;
           
            
            if numGens == 0
                error('No generators specified. size(W,1) = 0')
            end

            if ~all(cellfun(@isempty,U))
                
               
                Tuu = cone.BuildMultMap(U,U);
                Tuv = cone.BuildMultMap(U,V);
                Tvv = cone.BuildMultMap(V,V);
                
                [s,e] = cone.flqrIndx();
                Tvv(s:e,:) = 0;
                Tuv(s:e,:) = 0;
                
                Aeq_c1 = [W*(Tuu*A'),W*(Tuu*c(:))]';
                Aeq_c2 = [Tvv*A',Tvv*c(:)]';
                Aeq_c3 = [Tuv*A',Tuv*c(:)]';
                
            else
                
                Aeq_c1 = [W*(A'),W*(c(:))]';
                Aeq_c2 = []';
                Aeq_c3 = [];
                
            end
                    
            Aeq = [Aeq_c1,Aeq_c2,Aeq_c3];
            numFree = size(Aeq,2)-numGens;
            Aeq = [Aeq,sparse(numEq,numGens)];
            beq = sparse(size(Aeq,1),1);

            solMap.coneGens.s = 1;
            solMap.coneGens.e = numGens;
            solMap.shat_vvt.s = solMap.coneGens.e + 1;
            solMap.shat_vvt.e = solMap.shat_vvt.s + size(Aeq_c2,2)-1;
            solMap.beta_uvt.s = solMap.shat_vvt.e + 1; 
            solMap.beta_uvt.e = solMap.beta_uvt.s + size(Aeq_c3,2)-1;
            
            [Aeq,beq] = cleanLinear(Aeq,beq);
        
            Aineq = [-speye(numGens),sparse(numGens,numFree),speye(numGens)];
            bineq = zeros(size(Aineq,1),1);
            c = [zeros(numGens+numFree,1);-ones(numGens,1)];
            ubnd = [Inf*ones(numGens+numFree,1);ones(numGens,1)];
            lbnd = [-Inf*ones(numGens+numFree,1);zeros(numGens,1)];

        end

        function [c,Aineq,bineq,Aeq,beq,ubnd,lbnd,solMap] = BuildPrimLP(A,b,U,W,cone,Kface)

            %equation: \sum_i A_i mu_i = \sum_i s_i W_i
            % s_i \ge t_i \ge 0
            %cols: (mu,y,t)
            numGens = size(W,1);
            if numGens == 0
                error('No generators specified. size(W,1) = 0')
            end

            numY = size(A,1);
            solMap.y.s = numGens + 1;
            solMap.y.e = numY+solMap.y.s-1;

            
            if ~all(cellfun(@isempty,U))
                Tuu = cone.BuildMultMap(U,U);
                face = ConeBase(Kface);
                A = [Tuu*A']';
            else
                face = cone;
            end
            
            A = face.LowerTri(A);
            W = face.LowerTri(W);
            
            Aeq1 = [-W',A'];
            Aeq2 = [sparse(1,numGens),b'];

            Aeq = [Aeq1;Aeq2];
            Aeq = [Aeq,sparse(size(Aeq,1),numGens)];
            beq = sparse(size(Aeq,1),1);

            %Remove dependent equations.
            [Aeq,beq] = cleanLinear(Aeq,beq);

            Aineq = [-speye(numGens),sparse(numGens,numY),speye(numGens)];
            bineq = sparse(size(Aineq,1),1);
            c = [sparse(numGens+numY,1);-ones(numGens,1)];
            ubnd = [Inf*ones(numGens+numY,1);ones(numGens,1)];
            lbnd = [-Inf*ones(numGens+numY,1);sparse(numGens,1)];

        end

    function [success,U,Kface,redCert] = PolyhedralPrimIter(self,U,gensOrType,cone,Kface)

            W = [];  redCert = [];
            if ischar(gensOrType)
                W = facialRed.GetGenerators(Kface,gensOrType);
            else
                W = gensOrType;
            end
          
            numGens = size(W,1);
            success = 0;
   
            if (numGens == 0)
                return;
            end

            [cost,Aineq,bineq,Aeq,beq,ubnd,lbnd,solMap] = facialRed.BuildPrimLP(self.A,self.b,U,W,cone,Kface);
            [x,numerr,infeas] = facialRed.solveLP(cost,Aineq,bineq,Aeq,beq,lbnd,ubnd);
            
            success = numerr == 0 & infeas == 0 & -cost'*x > .1;

            x = sparse(x);
            
            if success

                x = sparse(x);
                xtemp = sparse(x(1:numGens,1));
                
                y = x(solMap.y.s:solMap.y.e);
                redCert = y'*self.A;
                xtemp = xtemp > max(xtemp)*.0001;
                sBarExtRay = W'*xtemp;         
                
                [U,~,Kface] = facialRed.Reduce(sBarExtRay,U,[],Kface);
                
            end

        end

        function W = GetGenerators(Kface,type)
            
            Z = ConeApprox(Kface);
                
            switch type

                case 'd'
                     W = Z.extRaysD();

                case 'dd'
                     W = Z.extRaysDD();

            end

        end

        function [success,U,V,Kface,redCert,solMap] = PolyhedralDualIter(self,U,V,gensOrType,cone,Kface)

            W = []; redCert = [];
            if ischar(gensOrType)
                Z = ConeApprox(Kface);
                
                switch gensOrType
                    
                    case 'd'
                         W = Z.extRaysD();
                         
                    case 'dd'
                         W = Z.extRaysDD();
                    
                end
                
            else
                W = gensOrType;
            end
          
            numGens = size(W,1);
            success = 0;
   
            if (numGens == 0)
                return;
            end
            

            
            [cost,Aineq,bineq,Aeq,beq,ubnd,lbnd,solMap] = facialRed.BuildDualLP(self.A,self.c,U,V,W,cone);
          
                        
           % if gensOrType == 'd' &&  all(cellfun(@isempty,U))   
                
                %indxDiag = cell2mat(cone.indxDiag); 
                %temp = [self.A;self.c(:)'];
                %x = all(temp(:,indxDiag)==0,1)';
               % success = nnz(x) > 0;
                
           % else
           %     return
                [x,numerr,infeas] = facialRed.solveLP(cost,Aineq,bineq,Aeq,beq,lbnd,ubnd);
                success = numerr == 0 & infeas == 0 & -cost'*x >= .2;
                
          %  end
            
            if success

                x = sparse(x);
                xtemp = sparse(x(1:numGens,1));

                sBar = W'*xtemp; 
                xtemp = xtemp > max(xtemp)*.0001;
                sBarExtRay = W'*xtemp;
                
                sHat = x(solMap.shat_vvt.s:solMap.shat_vvt.e);
                beta = x(solMap.beta_uvt.s:solMap.beta_uvt.e);
                
                if ~all(cellfun(@isempty,U))     
                    redCert = cone.ConjBlock2by2(sBar,sHat,beta/2,U,V);
                else
                    redCert = sBar;
                end
                
                [U,V,Kface] = facialRed.ReduceDual(sBarExtRay,U,V,Kface);
                
            end

        end

        function [U,V,K] = ReduceDual(SblkDiag,U,V,Kface)
            
            K = Kface;
            Z = ConeBase(Kface);
            for i = 1:length(Kface.s)
                [s,e] = Z.GetIndx('s',i);
                S = mat(SblkDiag(s:e));
                [B,rangeS] = nullqr(S);
                if (~isempty(U{i}))
                    V{i} = [V{i},U{i} * rangeS];
                    U{i} = U{i} * B; 
                else
                    V{i} = rangeS;
                    U{i} = B;  
                end
                K.s(i) = size(U{i},2); 
            end

        end


        function [U,V,K] = Reduce(SblkDiag,U,V,Kface)
            
            K = Kface;
            Z = ConeBase(Kface);
            for i = 1:length(Kface.s)
                [s,e] = Z.GetIndx('s',i);
                S = mat(SblkDiag(s:e));
                [B,rangeS] = nullqr(S);
                if (any(size(U{i}) ~= 0))
                    if ~isempty(V)
                        V{i} = [V{i},U{i} * rangeS];
                    end
                    U{i} = U{i} * B; 

                else
                    if ~isempty(V)
                        V{i} = rangeS;
                    end
                    U{i} = B;  
                end
                K.s(i) = size(U{i},2); 
            end

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
