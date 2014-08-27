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
        defaultRedOpts
    end

    methods

        function self = frlibPrg(A,b,c,K)
          
            if ~isstruct(K)
               error('Invalid input. 4th argument must be a struct')
            end
 
            self.cone = coneBase(K); 
             
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
            self.defaultRedOpts = [];
            
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
                     
            if isfield(opts,'verbose')
                verbose = opts.verbose;
            else
                verbose = 1;
            end
             
            if isfield(opts,'removeDualEq') 
                removeDualEq = opts.removeDualEq;
            else
                removeDualEq = 0;
            end
            
            [A,b,Ty1] = CleanLinear(self.A,self.b,useQR); 
            y0 = 0;

            if removeDualEq && self.K.f > 0 && any(self.K.s > 0)
                
                [A,b,c,K,Ty2,y0] = RemoveDualEquations(A,b,self.c,self.K);
                y0 = Ty1*y0;
                Ty = Ty1*Ty2; 
                
            else

                removeDualEq = 0;
                c = self.c;
                K = self.K;
                Ty = Ty1;
                
            end
            
            if size(A,1) ~= 0
                
                if (verbose == 1)
                    pars.fid = 1;
                else
                    pars.fid = 0;
                end

                [x,y,info] = sedumi(A,b,c,K,pars);
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

                xf = LinEqSol([self.A(:,indxFree);self.c(indxFree)],[eqError;costError]);
                x = [xf;x];
                
            end

        end

        function [x,y] = SolveMosek(self)
                
            [A,b,T] = CleanLinear(self.A,self.b); 
            [Ar,br,cr,Kr,Tr,y0] = RemoveDualEquations(A,b,self.c,self.K); 
            Z = coneBase(Kr);
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
            
           pass = solUtil.CheckPrimal(x,self.A,self.b,self.c,self.K,eps); 
           
        end
        
        function pass = CheckDual(self,y,eps)
            
           pass = solUtil.CheckDual(y,self.A,self.b,self.c,self.K,eps);  
           
        end

        function [prg] = ReducePrimal(self,method,opts)
            
            procReduce = [];
            if (strcmp(method,'d'))
                procReduce = @(self,U,V,cone,Kface) facialRed.PolyhedralPrimIter(self,U,V,'d',cone,Kface);
            end

            if strcmp(method,'dd')
                procReduce = @(self,U,V,cone,Kface) facialRed.PolyhedralPrimIter(self,U,V,'dd',cone,Kface);
            end
                                 
            if isempty(procReduce)
                error('Valid approximation not specified.');
            end

            if ~exist('opts','var')
                opts = self.defaultRedOpts;
            end
            
            maxIter =  ParseRedOpts(opts);
            
            Kface = self.K;
            U = cell(1,length(self.K.s)); V = U;
            Uarry = {};  Varry = {};  Sarry = {}; Karry = {}; yRedarry = {}; redTimeArry={};

            iter = 1;
            while (1)
                 
                [success,U,V,Kface,S,yRed,timeRed] = procReduce(self,U,V,self.cone,Kface);
                if (success)
                    Uarry{end+1} = U;
                    Varry{end+1} = V;
                    Sarry{end+1} = S;
                    yRedarry{end+1} = yRed;
                    Karry{end+1} = Kface;
                end

                redTimeArry{end+1} = timeRed;
                if success == 0 || iter >= maxIter
                    break;
                end
                
                iter = iter + 1;
                
            end

            prg = reducedPrimalPrg(self.A,self.b,self.c,self.K,Karry,Uarry,Varry,Sarry,yRedarry,redTimeArry,opts);

        end
        
        function [prg] = ReduceDual(self,method,opts)
    
            procReduce = [];
            if (strcmp(method,'d'))
                procReduce = @(self,U,V,cone,Kface) facialRed.PolyhedralDualIter(self,U,V,'d',cone,Kface);
            end

            if strcmp(method,'dd')
                procReduce = @(self,U,V,cone,Kface) facialRed.PolyhedralDualIter(self,U,V,'dd',cone,Kface);
            end

            if isempty(procReduce)
                error('Valid approximation not specified.');
            end

            if ~exist('opts','var')
                opts = self.defaultRedOpts;
            end
            
            maxIter = ParseRedOpts(opts);
                
            Kface = self.K;
            U = cell(1,length(self.K.s)); V = U;
            Uarry = {}; Sarry = {}; Karry = {}; Varry = {}; redTimeArry={};

            iter = 1;
            while (1)
                
                [success,U,V,Kface,S,timeRed] = procReduce(self,U,V,self.cone,Kface);
                if (success)
                    Uarry{end+1} = U;
                    Sarry{end+1} = S;
                    Varry{end+1} = V;
                    Karry{end+1} = Kface;
                end

                redTimeArry{end+1} = timeRed;
                if success == 0 || iter >= maxIter
                    break;
                end
                
                iter = iter + 1;
                
            end

            prg = reducedDualPrg(self.A,self.b,self.c,self.K,Karry,Uarry,Varry,Sarry,redTimeArry,opts);
            
        end
 
    end

end

function [maxIter] = ParseRedOpts(opts)

    if isfield(opts,'maxIter')
        maxIter = opts.maxIter;
    else
        maxIter = 10^8;
    end

end

