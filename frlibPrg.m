classdef frlibPrg

    properties
        A
        b
        c
        K
    end

    properties(GetAccess=protected)
        cone
        defaultSolveOpts
    end

    methods

        function self = frlibPrg(A,b,c,K)
          
            if ~isstruct(K)
               error('Invalid input. 4th argument must be a struct')
            end
 
            self.cone = ConeBase(K); 
             
            if  isempty(b)
                b = zeros(size(A,1),1);
            end

            if  isempty(c)
                c = zeros(size(A,2),1);
            end
            
            if self.cone.NumVar ~= size(A,2)
                error(['Cone sizes do not match columns of A. Num vars: ',num2str(self.cone.NumVar),' Num Col: ', num2str(size(A,2))]);
            end

            if length(b) ~= size(A,1)
                error(['Number of rows of A do not match length of b. Num rows: ',num2str(size(A,1)),', Length b: ', num2str(length(b))]);
            end
             
            if self.cone.NumVar ~= length(c)
                error(['Number of variables do not match length of c. Num vars: ',num2str(self.cone.NumVar),', Length c: ', num2str(length(c))]);
            end
          
            A = self.cone.Symmetrize(A); 
            c = self.cone.Symmetrize(c(:)'); 

            self.A = A;
            self.c = c;
            self.b = b;
            self.K = self.cone.K;
            self.defaultSolveOpts = [];

        end

        function [x,y,info] = Solve(self,opts)
            
            if ~exist('opts','var')
                opts = self.defaultSolveOpts;
            end
            
            if isfield(opts,'useQR')
                useQR = opts.useQR;
            else
                useQR = 0;
            end
            
            if isfield(opts,'removeDualEq')
                removeDualEq = opts.removeDualEq;
            else
                removeDualEq = 0;
            end
            
            [A,b,Ty1] = cleanLinear(self.A,self.b,useQR); 
            y0 = 0;

            if removeDualEq & self.K.f > 0
                [A,b,c,K,Ty2,y0] = RemoveDualEquations(A,b,self.c,self.K);
                y0 = Ty1*y0;
                Ty = Ty1*Ty2; 
            else
                c = self.c;
                K = self.K;
                Ty = Ty1;
            end
            
            if size(A,1) ~= 0
                [x,y,info] = sedumi(A,b,c,K);
                y = Ty*y + y0;
            else
                x = sparse(size(A,2),1); y = 0;
                y = y0;
            end

            if removeDualEq
                indxNotFree = (self.K.f+1):size(self.A,2);
                indxFree = 1:self.K.f;
                costNotFree = self.c(indxNotFree);
                costError = (c(:)-costNotFree(:))'*x;
                eqError = self.b-self.A(:,indxNotFree)*x;

                xf = lsol([self.A(:,indxFree);self.c(indxFree)],[eqError;costError]);
                x = [xf;x];
            end

        end

        function [y] = SolveDual(self)

            [A,b,T] = cleanLinear(self.A,self.b); 

            if (self.K.f > 0) 

                [Ar,br,cr,Kr,Tr,y0] = RemoveDualEquations(A,b,self.c,self.K);
                [~,yReduced] = sedumi(Ar,br,cr,Kr);
                
                y = Tr*yReduced+y0;

            else
                [x,y,info] = sedumi(A,b,self.c,self.K);
            end
            
            y = T*y;

        end

        function [x,y] = SolveMosek(self)
                
            [A,b,T] = cleanLinear(self.A,self.b); 
            [Ar,br,cr,Kr,Tr,y0] = RemoveDualEquations(A,b,self.c,self.K); 
            Z = ConeBase(Kr);
            A = Z.Desymmetrize(Ar);
            c = Z.Desymmetrize(cr(:)');
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

            [x,y] = spot_mosek(A,br,c(:),K,struct('verbose',1));

            y = T*(Tr*y+y0);
            
        end 
      
        function pass = CheckSolution(self,x,y,eps)
           pass = 1;
           pass = pass & self.CheckPrimal(x,eps);
           pass = pass & self.CheckDual(y,eps);
           pass = pass & norm(self.c(:)'*x-self.b'*y) < eps;
        end

        function pass = CheckPrimal(self,x,eps)
           pass = SolUtil.CheckPrimal(x,self.A,self.b,self.c,self.K,eps); 
        end
        
        function pass = CheckDual(self,y,eps)
           pass = SolUtil.CheckDual(y,self.A,self.b,self.c,self.K,eps);  
        end

        function [prg] = ReducePrimal(self,method)
            
            if (strcmp(method,'d'))
                procReduce = @(self,U,cone,Kface) facialRed.PolyhedralPrimIter(self,U,'d',cone,Kface);
            end

            if strcmp(method,'dd')
                 procReduce = @(self,U,cone,Kface) facialRed.PolyhedralPrimIter(self,U,'dd',cone,Kface);
            end

            A = self.A; b = self.b; c = self.c; Kface = self.K;
            U = cell(1,length(self.K.s)); V = U;
            Uarry = {}; Sarry = {}; Karry = {}; 

            while (1)
                
                [success,U,Kface,S] = procReduce(self,U,self.cone,Kface);
                if (success)
                    Uarry{end+1} = U;
                    Sarry{end+1} = S;
                    Karry{end+1} = Kface;
                end

                if success == 0
                    break;
                end
                
            end

            prg = reducedPrimalPrg(A,b,c,self.K,Karry,Uarry,Sarry);

        end

        function [prg] = ReduceDual(self,method)

            if (strcmp(method,'d'))
                procReduce = @(self,U,V,cone,Kface) facialRed.PolyhedralDualIter(self,U,V,'d',cone,Kface);
            end

            if strcmp(method,'dd')
                 procReduce = @(self,U,V,cone,Kface) facialRed.PolyhedralDualIter(self,U,V,'dd',cone,Kface);
            end

            A = self.A; b = self.b; c = self.c; Kface = self.K;
            U = cell(1,length(self.K.s)); V = U;
            Uarry = {}; Sarry = {}; Karry = {}; Varry = {};

            while (1)
                
                [success,U,V,Kface,S] = procReduce(self,U,V,self.cone,Kface);
                if (success)
                    Uarry{end+1} = U;
                    Sarry{end+1} = S;
                    Varry{end+1} = V;
                    Karry{end+1} = Kface;
                end

                if success == 0
                    break;
                end
                
            end

            prg = reducedDualPrg(A,b,c,self.K,Karry,Uarry,Varry,Sarry);

        end
        
    end

end
