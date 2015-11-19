classdef facialRed

    methods(Static)

        function [c,Aineq,bineq,Aeq,beq,ubnd,lbnd,solMap] = BuildDualLP(A,c,face,W)

            numGens = size(W,1);
            numEq = size(A,1)+1;

            if numGens == 0
                error('No generators specified. size(W,1) = 0')
            end

            if face.isProper
                 
                Aeq_c1 = [(face.coneToFace*A')'*W';(face.coneToFace*c')'*W'];
                Aeq_c2 = [A*face.spanConjFace';c*face.spanConjFace'];
                Aeq_c3 = [A*face.resSubspace';c*face.resSubspace'];
                
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

        function [c,Aineq,bineq,Aeq,beq,ubnd,lbnd,solMap] = BuildPrimLP(A,b,face,W)

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

            if face.isProper
                A = [face.coneToFace*A']';
            end
            
            A = face.parser.LowerTri(A);
            W = face.parser.LowerTri(W);
            
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

        function W = GetGeneratorsOfApprox(face,gensOrType)
            
            W = [];  
            if ischar(gensOrType)
                W = facialRed.GetGenerators(face.K,gensOrType);
            else
                W = gensOrType;
            end
         
        end
        
        function [success,newFace,timeRed]  = PolyhedralPrimIter(self,face,gensOrType)

            success = 0; newFace = []; timeRed = 0; 
            
            W =  facialRed.GetGeneratorsOfApprox(face,gensOrType);
            numGens = size(W,1);
            if numGens == 0, return, end
            
            [cost,Aineq,bineq,Aeq,beq,ubnd,lbnd,solMap] = facialRed.BuildPrimLP(self.A,self.b,face,W);
            
           
            [x,numerr,infeas,time] = LPSolver.SolveLP(cost,Aineq,bineq,Aeq,beq,lbnd,ubnd);
            
            timeRed = timeRed + time;
            success = numerr == 0 & infeas == 0 & -cost'*x >= .2;
   
            if success

                x = sparse(x);
                xtemp = sparse(x(1:numGens,1));
                
                yRed = x(solMap.y.s:solMap.y.e);
                
                xtemp = xtemp > max(xtemp)*.0001;
                sBarExtRay = W'*xtemp;         
                
                [newFace] = face.Intersect(sBarExtRay);
                newFace.redCert.y = yRed;
                newFace.redCert.S = yRed'*self.A;
                newFace.redCert.extRays = W(xtemp,:)';
                newFace.time = time;
                
            end

        end
             
        function [success,newFace,timeRed] = PolyhedralDualIter(self,face,gensOrType)

            success = 0; newFace = []; timeRed = 0; 
            W =  facialRed.GetGeneratorsOfApprox(face,gensOrType);
            numGens = size(W,1);
            if numGens == 0, return, end
            
            [cost,Aineq,bineq,Aeq,beq,ubnd,lbnd,solMap] = facialRed.BuildDualLP(self.A,self.c,face,W);
            
            [x,numerr,infeas,time] = LPSolver.SolveLP(cost,Aineq,bineq,Aeq,beq,lbnd,ubnd);
            timeRed = timeRed + time;
            
            success = numerr == 0 & infeas == 0 & -cost'*x >= .2;
                   
            if success

                x = sparse(x);
                xtemp = sparse(x(1:numGens,1));

                sBar = W'*xtemp; 
                xtemp = xtemp > max(xtemp)*.0001;
                sBarExtRay = W'*xtemp;
                
                sHat = x(solMap.shat_vvt.s:solMap.shat_vvt.e);
                beta = x(solMap.beta_uvt.s:solMap.beta_uvt.e);
                
                if face.isProper   
                    redCert = face.coneToFace'*sBar(:)+face.resSubspace'*beta(:)+face.spanConjFace'*sHat(:);
                else
                    redCert = sBar;
                end
                
                [newFace] = face.Intersect(sBarExtRay);
                newFace.redCert.S = redCert;
                newFace.time = time;
                
            end

        end
                 
         function [success,newFace,timeRed] = SDDPrimIter(self,face,cone)
            success = 0; newFace = [];      
            coneFace = coneBase(face.K);
            prg = spotsosprog;
            m = size(self.A,1);
            redCert = []; yRed = []; timeRed = [];
            
            A = self.A;
            if face.isProper
                A = [face.coneToFace*A']';
            end
            
            %A = face.parser.LowerTri(A);
            
            
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
           
            if success


                [newFace] = face.Intersect(redCert);
                newFace.redCert.y = yRed;
                newFace.redCert.S = redCert;
                newFace.redCert.extRays = redCert;
                newFace.time = [];
                
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
   
    end
end
