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
            self.c = c(:)';
            self.b = b(:);
            self.K = self.cone.K;
            self.defaultSolveOpts = [];
            self.defaultRedOpts = [];
            
        end

        function [x,y,info] = Solve(self,pars)

            if ~exist('pars','var')
                pars = [];
            end

            x = []; y = []; info = [];
            
            if (size(self.A,1) > 0)
            
                if ~isempty(which('sedumi'))
                   options.verbose = 1;
                   options.solver_options = [];
                    A = self.cone.Desymmetrize(self.A);
                    c = self.cone.Desymmetrize(self.c(:)');
                    try
                        [x,y,~,info] = spot_mosek(A,self.b,c,self.K,options);
                    catch
                        [x,y,info] = sedumi(self.A,self.b,self.c,self.K,pars);
                    end
                else
                   error('SeDuMi not found.'); 
                end
            
            else

                x = sparse(size(self.A,2),1); y = 0;
                y = 0;

            end
            
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
                procReduce = @(self,face) facialRed.PolyhedralPrimIter(self,face,'d');
            end

            if strcmp(method,'dd')
                procReduce = @(self,face) facialRed.PolyhedralPrimIter(self,face,'dd');
            end
            
            if ~exist('opts','var')
                opts = self.defaultRedOpts;
            end                          
                    
            frlibPrg.CheckInputs(procReduce,method);
            faces = self.Reduce(procReduce,opts);
            prg = reducedPrimalPrg(self,faces,opts);

        end
        
        function [prg] = ReduceDual(self,method,opts)
    
            procReduce = [];
            if (strcmp(method,'d'))
                procReduce = @(self,face) facialRed.PolyhedralDualIter(self,face,'d');
            end

            if strcmp(method,'dd')
                procReduce = @(self,face) facialRed.PolyhedralDualIter(self,face,'dd');
            end
            
            if ~exist('opts','var')
                opts = self.defaultRedOpts;
            end
              
            frlibPrg.CheckInputs(procReduce,method);
            faces = self.Reduce(procReduce,opts);
            prg = reducedDualPrg(self,faces,opts);
            
        end
          

        function x = FindDualFace(self)
         
            e = sparse(self.cone.indxNNeg,1,1,self.cone.NumVar,1)';
            Aaux = [self.A;self.c;e];
            baux = [zeros(size(Aaux,1)-1,1);1];
            temp = frlibPrg(Aaux,baux,[],self.K);
            [x,y] = temp.Solve();

        end
        
        function y = FindPrimalFace(self)
         
            K = self.K;
            K.f = self.K.f + 1;
            Caux = [zeros(self.cone.NumVar+1,1)];
            Aaux = [self.b,self.A];
            temp = frlibPrg(Aaux,0*self.b,Caux,K);
            [x,y] = temp.Solve();

        end


    end
     
    methods(Access=protected)
        
       function faces = Reduce(self,procReduce,opts) 
                
            maxIter =  frlibPrg.ParseRedOpts(opts);
           
            currentFace = faceBase(self.cone,self.cone.K);
            faces = {currentFace};
            
            iter = 1;

            while (1)

                [success,currentFace] = procReduce(self,currentFace);
                if (success)
                    faces{end+1} = currentFace;
                end

                if success == 0 || iter >= maxIter
                    break;
                end

                iter = iter + 1;
               

            end

        end
       
        
    end
    
    methods(Static,Access=protected)
        
        function [maxIter] = ParseRedOpts(opts)

            if isfield(opts,'maxIter')
                maxIter = opts.maxIter;
            else
                maxIter = 10^8;
            end

        end


        function CheckInputs(procReduce,method)

            if isempty(procReduce)
                if isa(method,'char')
                    method = ['''',method,''''];
                end
                error(['Specified approximation ',num2str(method),' not recognized. Use ''d'' or ''dd''.']);
            end

        end


        

    end
    
end







