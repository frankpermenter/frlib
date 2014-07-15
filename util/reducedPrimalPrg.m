classdef reducedPrimalPrg < frlibPrg

    properties

        Uarry
        Varry
        Sarry
        yRedarry
        Karry
        Aorig
        corig
        Korig
        noReductions

    end

    methods

        
        function self = reducedPrimalPrg(A,b,c,K,Karry,U,S,yRedArry)
         
            if (length(U) > 0)
                
                cone = coneBase(K);

                Tuu = cone.BuildMultMap(U{end},U{end});
                Ar = [A*Tuu'];
                cr = [c*Tuu'];
                Kr = coneBase.cleanK(Karry{end});
                
                noReductions = 0;

            else
               
                Kr = K;
                Ar = A;
                cr = c;
               
                noReductions = 1;

            end

            
            self@frlibPrg(Ar,b,cr,Kr);
            self.Uarry = U;
            self.Sarry = S;
            self.yRedarry = yRedArry;
            self.Karry = Karry;
            self.Aorig = A;
            self.corig = c;
            self.Korig = K;
            self.noReductions = noReductions;
            
            
        end

        
        function [xr,yr,dual_recov_success] = Recover(self,x,y,eps)

            if (self.noReductions)
                xr = x; yr = y;
                return
            end

            if ~exist('eps','var')
                eps = 10^-4;
            end

            
            xr = self.RecoverPrimal(x);
            [yr,y0,dual_recov_success] = self.RecoverDual(y,eps);
            if (dual_recov_success ~= 1)
                yr = [];
            end
  
        end
        
        
        
        function [x] = RecoverPrimal(self,x)

            if (self.noReductions)
                return
            end

            face = coneBase(self.Karry{end});
            transposeU = 0;
            Tuu = face.BuildMultMap(self.Uarry{end},self.Uarry{end},transposeU);
            x = Tuu*x;

        end
        
        
        function [yr,y0,success] = RecoverDual(self,y0,eps)
 

            success = 1;
            if (self.noReductions)
                yr = y0; y0 = yr; 
                return
            end
                 
            if ~exist('eps','var')
                eps = 10^-4;
            end
            
            U = self.Uarry; 

            
            ComputeDelta = @(t,dir) t*dir'*self.Aorig;
            s0 = self.corig-y0'*self.Aorig;
            [sr,success,deltas] = solUtil.LineSearch(s0,U,self.yRedarry,self.Korig,ComputeDelta); 
            
            yr = y0;
            for i=1:length(U)
                yr = yr - deltas(i)*self.yRedarry{i};
            end
        
            schk =  self.corig-yr'*self.Aorig;
            
        end
        
        
        function PrintStats(self)
            facialRed.PrintStats(self.K,self.Korig);
        end
          

    end

end
