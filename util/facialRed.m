classdef facialRed

    methods(Static)

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
            
            [Aeq,beq] = CleanLinear(Aeq,beq);
        
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
                face = coneBase(Kface);
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
            [Aeq,beq] = CleanLinear(Aeq,beq);

            Aineq = [-speye(numGens),sparse(numGens,numY),speye(numGens)];
            bineq = sparse(size(Aineq,1),1);
            c = [sparse(numGens+numY,1);-ones(numGens,1)];
            ubnd = [Inf*ones(numGens+numY,1);ones(numGens,1)];
            lbnd = [-Inf*ones(numGens+numY,1);sparse(numGens,1)];

        end

        function [success,U,V,Kface,redCert,yRed,timeRed] = PolyhedralPrimIter(self,U,V,gensOrType,cone,Kface)

            redCert = []; yRed = []; timeRed = [];
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
            
            tic;
            [x,numerr,infeas] = LPSolver.SolveLP(cost,Aineq,bineq,Aeq,beq,lbnd,ubnd);
            timeRed = toc;
           
            success = numerr == 0 & infeas == 0 & -cost'*x > .1;

            x = sparse(x);
            
            if success

                x = sparse(x);
                xtemp = sparse(x(1:numGens,1));
                
                yRed = x(solMap.y.s:solMap.y.e);
                redCert = yRed'*self.A;
                xtemp = xtemp > max(xtemp)*.0001;
                sBarExtRay = W'*xtemp;         
                
                [U,V,Kface] = facialRed.Reduce(sBarExtRay,U,V,Kface);
                
            end

        end

        function [success,U,V,Kface,redCert,timeRed] = PolyhedralDualIter(self,U,V,gensOrType,cone,Kface)

            W = [];  redCert = [];
            if ischar(gensOrType)
                W = facialRed.GetGenerators(Kface,gensOrType);
            else
                W = gensOrType;
            end
          
            numGens = size(W,1);
            success = 0;
            timeRed = [];
            
            if (numGens == 0)
                return;
            end
            
            [cost,Aineq,bineq,Aeq,beq,ubnd,lbnd,solMap] = facialRed.BuildDualLP(self.A,self.c,U,V,W,cone);
            
            tic;
            [x,numerr,infeas] = LPSolver.SolveLP(cost,Aineq,bineq,Aeq,beq,lbnd,ubnd);
            timeRed = toc;
            
            success = numerr == 0 & infeas == 0 & -cost'*x >= .2;
                   
            if success

                x = sparse(x);
                xtemp = sparse(x(1:numGens,1));

                sBar = W'*xtemp; 
                xtemp = xtemp > max(xtemp)*.0001;
                sBarExtRay = W'*xtemp;
                
                sBar = sBar(cone.Kstart.s(1):end);
                sHat = x(solMap.shat_vvt.s:solMap.shat_vvt.e);
                beta = x(solMap.beta_uvt.s:solMap.beta_uvt.e);
                
                if ~all(cellfun(@isempty,U))     
                    redCert = cone.ConjBlock2by2(sBar,sHat,beta/2,U,V);
                else
                    s = cone.GetIndx('s',1);
                    temp = sparse(1,s-1);
                    redCert = [temp,sBar'];
                end
                
                [U,V,Kface] = facialRed.Reduce(sBarExtRay,U,V,Kface);
                
            end

        end
                 
        function [success,U,V,Kface,redCert,yRed,timeRed] = SDDPrimIter(self,U,V,cone,Kface)
                 
            coneFace = coneBase(Kface);
            prg = spotsosprog;
            m = size(self.A,1);
            redCert = []; yRed = []; timeRed = [];
            
            if ~all(cellfun(@isempty,U))
                
                Tuu = cone.BuildMultMap(U,U);
                A = self.A*Tuu';
 
            else
                
                A = self.A;
            end
            
            
            [prg,u] = prg.newFree(m);
            prg = prg.withEqs(u'*self.b);
            S = sparse(1,coneFace.Kstart.s(1)-1);
            for i=1:length(coneFace.K.s)  
                if (coneFace.K.s(i) > 1 )
                    [prg,Stemp] = prg.newSDD(coneFace.K.s(i));   
                    S = [S,Stemp(:)'];
                else
                    [prg,Stemp] = prg.newPos(coneFace.K.s(i)); 
                     S = [S,Stemp(:)'];
                end
            end
            
            prg = prg.withEqs(S-u'*A);
            tr = sum(S(coneFace.indxNNeg));
            prg = prg.withEqs(tr-1);
            
            try
                tic
                sol = prg.minimize(0,@spot_mosek);
                timeRed = toc;
                redCert = double(sol.eval(S));
                yRed =  double(sol.eval(u));
                success = 1;
            catch
                success = 0;
                redCert = [];
                
            end

            if (success == 0)
               return 
            end
            %Heuristic rounding
            redCert(find(round(redCert*10^6)==0)) = 0;
           
            [U,V,Kface] = facialRed.Reduce(redCert,U,V,Kface);
        
        
        end 
        
        
        
        
        
        function [U,V,K] = Reduce(SblkDiag,U,V,Kface)
            
            K = Kface;
            Z = coneBase(Kface);
            for i = 1:length(Kface.s)
                
                [s,e] = Z.GetIndx('s',i);

                S = reshape(SblkDiag(s:e),Kface.s(i),Kface.s(i));

                [B,rangeS] = NullQR(S);
                
                if (any(size(U{i}) ~= 0))
                  
                    V{i} = [V{i},U{i} * rangeS];
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

        function W = GetGenerators(Kface,type)
            
            Z = coneApprox(Kface);
            if Z.NumVar == 0
               W = [];
               return
            end
            
            
            switch type

                case 'd'
                     W = Z.extRaysD();

                case 'dd'
                     W = Z.extRaysDD();

            end

        end

        
        function PrintStats(K,Korig)
            display([sprintf('\t'),'Before: Size of PSD constraint (K.s):',sprintf('\t'),sprintf('%d ',Korig.s)])
            display([sprintf('\t'),'After:  Size of PSD constraint (K.s):',sprintf('\t'),sprintf('%d ',K.s)])
        end
        
        
    end
end
