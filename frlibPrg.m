classdef frlibPrg

    properties
        A
        b
        c
        K
    end

    properties(GetAccess=protected)
        Z
    end

    methods

        function self = frlibPrg(A,b,c,K)

            if ~isstruct(K)
               error('Invalid input. 4th argument must be a struct')
            end


            if isempty(c)
                c = sparse(size(A,2),1);
            end

            if isempty(b)
                b = sparse(size(A,1),1);
            end

            c = c(:)';
            [A,c] = makeSymmetric(A,c,K);

            self.Z = coneHelp(A,b,c,K); 

            self.K = self.Z.K;
            self.A = self.Z.A;
            self.b = self.Z.b(:);
            self.c = self.Z.c(:);

        end

        function [x,y,info] = Solve(self,opts)
            
            if ~exist('opts','var')
                opts = [];
            end
            
            if isfield(opts,'useQR')
                useQR = opts.useQR;
            else
                useQR = 0;
            end
            
            [A,b,T] = cleanLinear(self.A,self.b,useQR); 
            [x,y,info] = sedumi(A,b,self.c,self.K);
            y = T*y;

        end


        function [y] = SolveDual(self)

            [A,b,T] = cleanLinear(self.A,self.b); 

            if (self.K.f > 0) 

                [s,e] = self.Z.GetIndx('f',1);
                Af = self.A(:,s:e);
                cf = self.c(s:e);
%                 [y0,varsPresolved,tRemoveVars,tAddVars] = PreSolveLinearEq(Af',cf(:));
%                 cReduced = self.c(:)-A(varsPresolved,:)'*y0;
%                 AReduced = tRemoveVars*A; 
%                 bReduced = tRemoveVars*b;
%[~,yReduced] = sedumi(AReduced,bReduced,cReduced,self.K);
%                y = tAddVars*yReduced;
%                y(varsPresolved) = y0;

                [Ar,br,cr,Kr,Tr,y0] = RemoveDualEquations(A,b,self.c,self.K);
                
                
                [~,yReduced] = sedumi(Ar,br,cr,Kr);
                
                y = Tr*yReduced+y0


            else
                [x,y,info] = sedumi(A,b,self.c,self.K);
            end
            y = T*y;

        end

       function [x] = SolvePrimal(self)

            [A,b] = cleanLinear(self.A,self.b); 
            
            [x0,varsPresolved,tRemoveVars,tAddVars,AReduced,bReduced] = PreSolveLinearEq(A,b);
            cReduced = tRemoveVars*self.c(:);
      
            [x,~] = sedumi(AReduced,bReduced,cReduced,self.K);

            x = tAddVars*x;
            x(varsPresolved) = x0;

       end


       
        function [y] = SolveDualMosek(self)
                
            [A,b,T] = cleanLinear(self.A,self.b); 
            [Ar,br,cr,Kr,Tr,y0] = RemoveDualEquations(A,b,self.c,self.K); 
            Z = coneHelp(Ar,br,cr,Kr);
            A = Z.desymmetrize(Ar);
            c = Z.desymmetrize(cr(:)');
            info = [];

            K = Kr;
            if ~isfield(K,'f') 
                K.f = 0;
            end

            if ~isfield(K,'l') 
                K.l = 0;
            end

            if ~isfield(K,'q') || all(K.q == 0)
                K.q = [];
            end

            if ~isfield(K,'s') || all(K.s == 0)
               K.s = [];
            else
               K.s = K.s(K.s ~= 0); 
            end

            if ~isfield(K,'r') || all(K.r == 0)
                K.r = [];
            end

            [~,y] = spot_mosek(A,br,c(:),K);
            y = T*(Tr*y+y0);
        end 
       
       
       
        function [x,y,info] = SolveMosek(self)

   
            A = self.Z.desymmetrize(self.A);
            c = self.Z.desymmetrize(self.c(:)');
            info = [];

            K = self.K;
            if ~isfield(K,'f') 
                K.f = 0;
            end

            if ~isfield(K,'l') 
                K.l = 0;
            end

            if ~isfield(K,'q') || all(K.q == 0)
                K.q = [];
            end

            if ~isfield(K,'s') || all(K.s == 0)
               K.s = [];
            else
               K.s = K.s(K.s ~= 0); 
            end

            if ~isfield(K,'r') || all(K.r == 0)
                K.r = [];
            end

            [x,y] = spot_mosek(A,self.b,c(:),K,struct('verbose',1));

        end
        
        function pass = CheckPrimal(self,x,eps)
           pass = SolUtil.CheckPrimal(x,self.A,self.b,self.c,self.K,eps); 
           
        end
        
        function pass = CheckDual(self,y,eps)
           pass = SolUtil.CheckDual(y,self.A,self.b,self.c,self.K,eps);  
        end

        function [prg] = ReducePrimal(self,method)
            method = lower(method);

            if (strcmp(method,'sdd'))
                procReduce = @facialRed.SDDPrimIter;
            end

            if (strcmp(method,'d'))
                procReduce = @facialRed.DiagPrimIter;
            end

            if (strcmp(method,'dd'))
                procReduce = @facialRed.DiagDomPrimIter;
            end

            A = self.A; b = self.b; c = self.c; K = self.K;
            T = [];
            firstPass = 1; 
            Ua = [];
            Sa = [];
            
            while (1)

                [success,A,c,K,Tform,S] = feval(procReduce,A,b,c,K);
               % [A,b] = cleanLinear(A,b);

                if success == 0
                   break;
                end
                
                Ua{end+1} = Tform;
                Sa{end+1} = S;

                if iscell(Tform) 
                    for i=1:length(Tform)
                        if (firstPass)
                            U{i} = Tform{i};
                        else
                            U{i} = U{i}*Tform{i};
                        end
                    end
                else
                    if (firstPass)
                         U = Tform;
                     else
                         U = U*Tform;
                     end
                end

                firstPass = 0;

                T = U; 
              
            end

            prg = reducedPrimalPrg(A,b,c,K,Ua,Sa,self.Z);

        end

        function [prg] = ReduceDual(self,method)

            method = lower(method); 
            if (strcmp(method,'sdd'))
                procReduce = @facialRed.SDDDualIter;
            end

            if (strcmp(method,'d'))
                procReduce = @facialRed.DiagDualIter;
            end

            if strcmp(method,'dd')
                procReduce = @facialRed.DiagDomDualIter;
            end

            A = self.A; b = self.b; c = self.c; K = self.K;
            Deq = ones(0,size(A,1));feq=[];

            %convert these into equations on dual variables 
            if (self.K.f > 0)
                Deq = A(:,1:self.K.f)';
                feq = c(1:self.K.f);
                feq = feq(:);
                c = c(self.K.f+1:end);
                K.f = 0;
                A = A(:,self.K.f+1:end);
            end

            Scellarray= [];
            while (1)
                [success,A,c,K,Deq,feq,S] = feval(procReduce,A,c,K,Deq,feq);
                [Deq,feq] = cleanLinear(Deq,feq);
                Scellarray{end+1} = S;
                if success == 0
                    break;
                end
            end

            %convert linear constraints into free variables
            A = [Deq',A];
            c = [feq;c(:)]; 
            K.f = K.f + length(feq);
            prg = reducedDualPrg(A,b,c,K,[],Scellarray,'dual');

        end
    end

end
