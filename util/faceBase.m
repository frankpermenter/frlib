classdef faceBase 
    
    properties
       coneToFace
       resSubspace 
       spanConjFace 
       U
       V
       K
       isProper
       redCert
       parser
    end

    properties(Access=private)
       parserCone
    end

    methods
        
        function self = faceBase(cone,K,U,V)  

            self.isProper = any(K.s < cone.K.s);
            self.parserCone = cone;
            self.parser = coneBase(K);
            self.K = self.parser.K;

            if (self.isProper) && nargin > 2
                
                self.U = U;
                self.V = V;
                [self.coneToFace,self.spanConjFace,self.resSubspace] = self.ComputeMaps();
             
            end

        end
      
        function dim = spanConjFaceDim(self)
            dim = size(self.spanConjFace,1);                
        end 

        function dim = resSubspaceDim(self)
            dim = size(self.resSubspace,1);                
        end 

        function newFace = Intersect(self,SblkDiag)
            
            K = self.K;
   
            for i = 1:length(self.K.s)
                
                [s,e] = self.parser.GetIndx('s',i);

                S = reshape(SblkDiag(s:e),self.K.s(i),self.K.s(i));

                [B,rangeS] = NullQR(S);
                
                if self.isProper
                  
                    V{i} = [self.V{i},self.U{i} * rangeS];
                    U{i} = self.U{i} * B; 
                    
                else
                      
                    V{i} = rangeS;  
                    U{i} = B;  
                    
                end
                
                K.s(i) = size(U{i},2); 
                
            end
            
            newFace = faceBase(self.parserCone,K,U,V);
          
        end
        
        function conjFace = Conjugate(self)
            K.s = cellfun(@(V)  size(V,2),self.V);
            conjFace = coneBase(K);
        end
 
        function [coneToFace,spanConjFace,resSubspace] = ComputeMaps(self)
         
            x = cell(length(self.U),1);
            for i=1:length(self.U)
                x{i}.U = self.U{i};
                x{i}.V = self.V{i};
            end

            resSubspaceNotSym = cellfun(@(x)  kron(x.U,x.V),x,'UniformOutput',0);
            resSubspaceNotSym = blkdiag(resSubspaceNotSym{:})';
            resSubspaceNotSym = self.AddZeroFlrqCols(resSubspaceNotSym);
            resSubspace = self.parserCone.Symmetrize(resSubspaceNotSym);

            spanConjFace = cellfun(@(x)  kron(x.V,x.V),x,'UniformOutput',0);
            spanConjFace = blkdiag(spanConjFace{:})';
            spanConjFace = self.AddZeroFlrqCols(spanConjFace);
            
            conjFace = self.Conjugate();
            spanConjFace = 1/2*spanConjFace(conjFace.LowerTriIndx(),:)...
                         + 1/2*spanConjFace(conjFace.UpperTriIndx(),:);

            
            coneToFace = cellfun(@(x)  kron(x.U,x.U),x,'UniformOutput',0);
            coneToFace = blkdiag(coneToFace{:})';
            coneToFace = blkdiag(speye(self.parserCone.NumFlqrVars),coneToFace);

        end       

        function pass = InLinearSpan(self,x,eps)
            
            if ~exist('eps','var') || isempty(eps)
                eps = 10^-8;
            end
            
            pass = norm(self.resSubspace*x(:))  < eps;
            pass = pass & norm(self.spanConjFace*x(:)) < eps;
            
        end

        function pass = InDualCone(self,x,eps)
            
            if ~exist('eps','var') || isempty(eps)
                eps = 10^-8;
            end
           
            if self.isProper
                xface = full(self.coneToFace*x(:));
            else
                xface = full(x); 
            end
            
            for i=1:length(self.K.s)
                [s,e] = self.parser.GetIndx('s',i);
                xtest = solUtil.mat(xface(s:e));

                if (min(eig(xtest)) > -eps) 
                    pass(i) = 1;
                else
                    pass(i) = isempty(xtest);
                end
            
            end
            
            pass = all(pass);
            
        end

        function [DualEqA,DualEqC] = FacePerpEqs(self,A,c)
 
            DualEqA = A*self.spanConjFace';
            DualEqA = [DualEqA,A*self.resSubspace'];
            DualEqC = c*self.spanConjFace';
            DualEqC = [DualEqC,c*self.resSubspace'];
            
        end

    end

    methods(Access=protected)
        
        function mat = AddZeroFlrqCols(self,mat)
            
            [s,e] = self.parserCone.flqrIndx();
            numZeroCol = e-s+1;
            numZeroRow = size(mat,1);
            mat = [sparse(numZeroRow,numZeroCol),mat];
            
        end

    end
 
end