classdef reducedDualPrg < frlibPrg

    properties (SetAccess=protected)
        
        Uarry
        Varry
        Sarry
        Karry
        Aorig
        borig
        corig
        Korig
        noReductions
        lpSolveTime
        opts
        
    end

    properties(Access=protected)
        
        Ty
        y0
        AwithEq;
        bwithEq;
        cwithEq;
        KwithEq;
        dualEqMap;
       
    end
    
    methods

        function self = reducedDualPrg(A,b,c,K,Karry,U,V,S,timeRed,opts)
            
            c = c(:)';
            Aorig = A; borig = b; corig = c; Korig = K;
            
            if ~isfield(opts,'useQR')
                opts.useQR = 0;
            end
                               
            if ~isfield(opts,'removeDualEq') 
                opts.removeDualEq = 0;
            end

            noReductions = isempty(U);

            if ~(noReductions)
                 
                [A,b,Ty] = CleanLinear(A,b,opts.useQR);
                y0 = 0;
                
                cone = coneBase(K);
                Kr = coneBase.cleanK(Karry{end});

                [DualEqA,DualEqC,Tuu,~,~,dualEqMap] = GetDualEqs(cone,U{end},V{end},A,c);
                Auu = [A*Tuu'];
                cuu = [c(:)'*Tuu'];

                AwithEq = [DualEqA,Auu];
                cwithEq = [DualEqC,cuu];
                KwithEq = coneBase.cleanK(Karry{end});
                KwithEq.f = KwithEq.f + size(DualEqA,2);

                if ~opts.removeDualEq

                    Ar = AwithEq;
                    br = b;
                    cr = cwithEq;
                    Kr = KwithEq;
   
                else

                    [Ar,br,cr,Ttemp,y0temp] = PreSolveDualEqs(DualEqA,DualEqC,Auu,b,cuu);
                    y0 = y0temp;
                    Ty = Ttemp;

                end

            else
               
                y0 = [];
                Ty = [];
                Kr = K;
                Ar = A;
                br = b;
                cr = c;
                AwithEq = A;
                bwithEq = b;
                cwithEq = c;
                KwithEq = K;
                dualEqMap = {};
               
            end
           
            self@frlibPrg(Ar,br,cr,Kr);
            self.Uarry = U;
            self.Sarry = S;
            self.Karry = Karry;
            self.Varry = V;
            self.Aorig = Aorig;
            self.borig = borig;
            self.corig = corig;
            self.Korig = Korig;

            self.noReductions = noReductions;
            self.lpSolveTime = timeRed;
            
            self.AwithEq = AwithEq;
            self.bwithEq = b;
            self.cwithEq = cwithEq;
            self.KwithEq = KwithEq;
            self.dualEqMap = dualEqMap;
            self.defaultSolveOpts = [];
            self.opts = opts;
            self.Ty = Ty;
            self.y0 = y0;
           
        end
        
        function PrintStats(self)
            facialRed.PrintStats(self.K,self.Korig);
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

                e1 = self.dualEqMap.beta_uvt.e;
                e2 = self.dualEqMap.shat_vvt.e;
                n = max(e1,e2);
               
                costError = (self.cwithEq(n+1:end)-self.c(:)')*x;
                eqError = self.bwithEq-self.AwithEq(:,n+1:end)*x;
                
                xf = LinEqSol([self.AwithEq(:,1:n);self.cwithEq(:,1:n)],[eqError;costError]);
                x = [xf;x];

                if (norm(self.AwithEq*x-self.bwithEq) > eps)
                    warning('Recovery of free dual variables failed');
                    freeRecovFailed = 1; 
                end               
            end
            
            U = self.Uarry; V = self.Varry;
            solMap = self.dualEqMap;
            coneOrig = coneBase(self.Korig);
                        
            cone = coneBase(self.KwithEq);
            s = cone.GetIndx('s',1);
            sBar = x(s:end);
            sHat = x(solMap.shat_vvt.s:solMap.shat_vvt.e);

            Kconj.s = cellfun(@(V)  size(V,2),V{end});
            coneConj = coneBase(Kconj);
            sHat = coneConj.InitSymmetric(sHat(:)');
            
            beta = x(solMap.beta_uvt.s:solMap.beta_uvt.e);
              
            xflqr = x(solMap.beta_uvt.e+1:s-1);            
            if ~all(cellfun(@isempty,U))     
                [x,x11m,x12m,~] = coneOrig.ConjBlock2by2(sBar,sHat,beta/2,U{end},V{end});
            else
                x = sBar;
            end

            x(1:length(xflqr)) = xflqr;
            xr = x; 
            x0 = xr;
     
            if (freeRecovFailed)
               success = 0; 
               return 
            end
            
            
            
            if length(U) == 1
                for i=1:length(x12m)
                    if norm(x12m{i}'*NullQR(x11m{i}),'fro') > eps
                        success = 0;
                        return
                    end         
                end
            end

            ComputeDelta = @(t,dir) dir(:)*t;
            [xr,success] = solUtil.LineSearch(x0,U,self.Sarry,self.Korig,ComputeDelta,eps); 
            
        end
   
    end

end

function [DualEqA,DualEqC,Tuu,Tuv,Tvv,solMap] = GetDualEqs(cone,U,V,A,c)

    Tuu = cone.BuildMultMap(U,U);
    Tuv = cone.BuildMultMap(U,V);
    Tvv = cone.BuildMultMap(V,V);

    [~,e] = cone.flqrIndx();
    Tvv = Tvv(e+1:end,:);
    Tuv = Tuv(e+1:end,:);
   
    Kconj.s = cellfun(@(V)  size(V,2),V);
    coneConj = coneBase(Kconj);
    DualEq1 = coneConj.UpperTri(A*Tvv');
    DualEq2 = A*Tuv';

    solMap.shat_vvt.s = 1;
    solMap.shat_vvt.e = solMap.shat_vvt.s + size(DualEq1,2)-1;
    solMap.beta_uvt.s = solMap.shat_vvt.e + 1; 
    solMap.beta_uvt.e = solMap.beta_uvt.s + size(DualEq2,2)-1;

    DualEqA = [DualEq1,DualEq2];
    DualEqC = [coneConj.UpperTri(c(:)'*Tvv'),c(:)'*Tuv'];

end

function [A_presolve,b_presolve,c_presolve,T,y0] = PreSolveDualEqs(Deq,feq,A,b,c)
    
    [y0,T] = LinEqSol(Deq',feq');
   
    A_presolve = T'*A;
    c_presolve = c(:) - A'*y0;
    b_presolve = T'*b;

end




