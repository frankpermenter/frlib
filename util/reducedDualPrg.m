classdef reducedDualPrg < reducedPrg
   
    properties(Access=protected)
        
        prgWithEq;

    end

    methods

        function self = reducedDualPrg(unreducedPrg,faces,opts)
            
            if ~isfield(opts,'useQR')
                opts.useQR = 0;
            end

            if ~isfield(opts,'removeDualEq') 
                opts.removeDualEq = 0;
            end

            if faces{end}.isProper
                 
                y0 = 0; Ty = speye(size(unreducedPrg.A,1));
                
                Auu = [faces{end}.coneToFace*unreducedPrg.A']';
                cuu = [faces{end}.coneToFace*unreducedPrg.c']';

                [DualEqA,DualEqC] = faces{end}.FacePerpEqs(unreducedPrg.A,unreducedPrg.c);

                prgWithEq.A = [DualEqA,Auu];
                prgWithEq.b = unreducedPrg.b;
                prgWithEq.c = [DualEqC,cuu];
                prgWithEq.K = faces{end}.K;
                
                prgWithEq.K.f = prgWithEq.K.f + size(DualEqA,2);

                if ~opts.removeDualEq

                    redPrg = prgWithEq;
   
                else

                    [redPrg,Ttemp,y0temp] = PreSolveDualEqs(DualEqA,DualEqC,Auu,unreducedPrg.b,cuu);
                    redPrg.K = faces{end}.K;
                    y0 = y0temp;
                    Ty = Ttemp;

                end

            else
               
                y0 = [];
                Ty = [];

                redPrg = unreducedPrg;
                prgWithEq = unreducedPrg;
    
            end
           
            self@reducedPrg(redPrg.A,redPrg.b,redPrg.c,redPrg.K);
            self.faces = faces;
            self.unreducedPrg = unreducedPrg;    
            self.prgWithEq = prgWithEq;
            self.defaultSolveOpts = [];
            self.opts = opts;
            self.Ty = Ty;
            self.y0 = y0;
           
        end
        
        function [xr,yr,primal_recov_success,x0] = Recover(self,x,y,eps)

            if ~exist('eps','var')
                eps = 10^-4;
            end
            
            yr = self.RecoverDual(y);
            
            [xr,x0,primal_recov_success] = self.RecoverPrimal(x,eps);
            if (primal_recov_success ~= 1)
                xr = [];
            end
  
        end
             
        function yr = RecoverDual(self,y)
                    
            if (self.noReductions)
                yr = y;
                return
            end
            
            if (size(self.Ty,2) > 0)
                yr = self.Ty*y + self.y0;
            else
                yr = self.y0;
            end
            
        end
                
        function [xr,x0,success] = RecoverPrimal(self,x,eps)
            
            if ~exist('eps','var')
                eps = 10^-4;
            end
           
            if (self.noReductions)
                xr = x; x0 = x;  success = 1;
                return
            end
                 
            if ~exist('eps','var')
                eps = 10^-4;
            end
            
            freeRecovFailed = 0; 
            if self.opts.removeDualEq

                n = self.faces{end}.spanConjFaceDim + self.faces{end}.resSubspaceDim;
               
                costRef = self.c*x;
                costError = costRef-self.prgWithEq.c(n+1:end)*x;
                eqError = self.prgWithEq.b-self.prgWithEq.A(:,n+1:end)*x;
                
                xf = LinEqSol([self.prgWithEq.A(:,1:n);self.prgWithEq.c(:,1:n)],[eqError;costError]);
                x = [xf;x];

                freeRecovFailed = norm(self.prgWithEq.A*x-self.prgWithEq.b) > eps;
                freeRecovFailed = freeRecovFailed & norm(self.prgWithEq.c*x-costRef) > eps;        
                
            end

            dimSpanConj = size(self.faces{end}.spanConjFace,1);
            dimRes = size(self.faces{end}.resSubspace,1);

			R = x(1:dimSpanConj);
            Q = self.faces{end}.coneToFace'*x(dimSpanConj+dimRes+1:end);
            Z = x(dimSpanConj+1:dimRes+dimSpanConj);

            xperp = self.faces{end}.resSubspace'*Z(:)+self.faces{end}.spanConjFace'*R(:);
            x = Q(:)' + xperp(:)';
          
            xr = x; 
            x0 = xr;
     
            if (freeRecovFailed)
               success = 0; 
               return 
            end
            
            ComputeDelta = @(t,redCert) t*redCert.S;
            [xr,success] = solUtil.LineSearch(x0,self.faces,ComputeDelta,eps);  
            
        end
   
    end

end

function [redPrgData,T,y0] = PreSolveDualEqs(Deq,feq,A,b,c)
    
    [y0,T] = LinEqSol(Deq',feq');
   
    redPrgData.A = T'*A;
    redPrgData.c = c(:) - A'*y0;
    redPrgData.b = T'*b;

end




