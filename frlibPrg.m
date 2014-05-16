classdef frlibPrg

    properties
        A
        b
        c
        K
    end

    properties(GetAccess=protected)
        cone
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
          
           
            [A,c] = makeSymmetric(A,c(:)',K);  
            
            self.A = A;
            self.c = c;
            self.b = b;
            self.K = self.cone.K;
            
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

                [s,e] = self.cone.GetIndx('f',1);
                Af = self.A(:,s:e);
                cf = self.c(s:e);


                [Ar,br,cr,Kr,Tr,y0] = RemoveDualEquations(A,b,self.c,self.K);
                [~,yReduced] = sedumi(Ar,br,cr,Kr);
                
                y = Tr*yReduced+y0;


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
       
       
       
        function [x,y,info] = SolveMosek2(self)

            A = self.cone.Desymmetrize(self.A);
            c = self.cone.Desymmetrize(self.c(:)');
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

            method = lower(method); 
            if (strcmp(method,'sdd'))
                procReduce = @facialRed.SDDDualIter;
            end

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
