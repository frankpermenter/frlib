classdef reducedDualPrg < frlibPrg

    properties 
        Uarry
        Varry
        Sarry
        Karry
        solMap
        Aorig
        corig
        Korig
        noReductions
        lpSolveTime
    end

    methods

        function self = reducedDualPrg(A,b,c,K,Karry,U,V,S,timeRed)

            if (length(U) > 0)
                
                cone = coneBase(K);
              
                Tuu = cone.BuildMultMap(U{end},U{end});
                Tuv = cone.BuildMultMap(U{end},V{end});
                Tvv = cone.BuildMultMap(V{end},V{end});
           
                [s,e] = cone.flqrIndx();
                Tvv(s:e,:) = 0;
                Tuv(s:e,:) = 0;
                
                DualEq1 = (A*Tvv');
                DualEq2 = A*Tuv';
           
                DualEqA = [DualEq1,DualEq2];
                DualEqC = [(c(:)'*Tvv'),c(:)'*Tuv'];

                solMap.shat_vvt.s = 1;
                solMap.shat_vvt.e = solMap.shat_vvt.s + size(DualEq1,2)-1;
                solMap.beta_uvt.s = solMap.shat_vvt.e + 1; 
                solMap.beta_uvt.e = solMap.beta_uvt.s + size(DualEq2,2)-1;

                Ar = [DualEqA,A*Tuu'];
                cr = [DualEqC,c*Tuu'];
                Kr = coneBase.cleanK(Karry{end});
                Kr.f = Kr.f + size(DualEqA,2);
                noReductions = 0;

            else
               
                Kr = K;
                Ar = A;
                cr = c;
                solMap ={};
                noReductions = 1;

            end

            
            self@frlibPrg(Ar,b,cr,Kr);
            self.Uarry = U;
            self.Sarry = S;
            self.Karry = Karry;
            self.Varry = V;
            self.Aorig = A;
            self.corig = c;
            self.Korig = K;
            self.solMap = solMap;
            self.noReductions = noReductions;
            self.lpSolveTime = timeRed;
            opts.removeDualEq = 1;
            self.defaultSolveOpts = opts;

        end
        
        function PrintStats(self)
            facialRed.PrintStats(self.K,self.Korig);
        end
        
        function [xr,yr,primal_recov_success,x0] = Recover(self,x,y,eps)

            if (self.noReductions)
                xr = x; yr = y; primal_recov_success = 1; x0 = x;
                return
            end

            if ~exist('eps','var')
                eps = 10^-4;
            end

            yr = y;
            [xr,x0,primal_recov_success] = self.RecoverPrimal(x,eps);
            if (primal_recov_success ~= 1)
                xr = [];
            end
  
        end
      
        function [xr,x0,success] = RecoverPrimal(self,x,eps)

            success = 1;
            if (self.noReductions)
                xr = x; x0 = x; 
                return
            end
                 
            if ~exist('eps','var')
                eps = 10^-4;
            end
            
            U = self.Uarry; V = self.Varry;
            solMap = self.solMap;
            coneOrig = coneBase(self.Korig);
            
            s = self.cone.GetIndx('s',1);
            sBar = x(s:end);
            sHat = x(solMap.shat_vvt.s:solMap.shat_vvt.e);
            beta = x(solMap.beta_uvt.s:solMap.beta_uvt.e);
              
            xflqr = x(solMap.beta_uvt.e+1:s-1);            
            if ~all(cellfun(@isempty,U))     
                [x,x11m,x12m,x22m] = coneOrig.ConjBlock2by2(sBar,sHat,beta/2,U{end},V{end});
            else
                x = sBar;
            end

            x(1:length(xflqr)) = xflqr;
            xr = x; 
            x0 = xr;
     
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
