classdef reducedPrimalPrg < frlibPrg

    properties(SetAccess=protected)

        Uarry
        Varry
        Sarry
        yRedarry
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
       
    end

    methods
 
        function self = reducedPrimalPrg(A,b,c,K,Karry,U,V,S,yRedArry,timeRed,opts)
     
            c = c(:)';
            Aorig = A; borig = b; corig = c; Korig = K;

            if ~exist('opts','var')
                opts.useQR = 0;
            end

            if ~isfield(opts,'useQR')
                opts.useQR = 0;
            end
            
            if ~isfield(opts,'quiet')
                opts.quiet = 0;
            end
            
            noReductions = isempty(U);
            
            if ~(noReductions)

                cone = coneBase(K);
                Tuu = cone.BuildMultMap(U{end},U{end});
                Ar = [A*Tuu'];
                cr = [c*Tuu'];
                Kr = coneBase.cleanK(Karry{end});

                [Ar,br,Ty] = CleanLinear(Ar,b,opts.useQR);
                
                if  opts.useQR == 0 && size(Ar,2) < size(Ar,1)
                    if  opts.quiet == 0
                        warning(['frlib: reduced problem has more equations '...
                        'than primal variables.  Call Reduce(<approx>,opts) '..., 
                        'with opts.useQR = 1 to remove linearly dependent equations.']); 
                    end
                end
                
            else

                Kr = K;
                Ar = A;
                br = b;
                cr = c;
                Ty  = speye(size(A,1));
                 
            end

            self@frlibPrg(Ar,br,cr,Kr);
            self.Uarry = U;
            self.Varry = V;
            self.Sarry = S;
            self.yRedarry = yRedArry;
            self.Karry = Karry;
            
            self.Aorig = Aorig;
            self.borig = borig;
            self.corig = corig;
            self.Korig = Korig;
            self.lpSolveTime = timeRed;
            self.noReductions = noReductions;
            self.Ty = Ty;
            self.opts = opts;      
     
        end

        
        function [xr,yr,dual_recov_success] = Recover(self,x,y,eps)

            if (self.noReductions)
                xr = x; yr = y; dual_recov_success = 1;
                return
            end

            if ~exist('eps','var')
                eps = 10^-4;
            end
           
            xr = self.RecoverPrimal(x);
            [yr,~,dual_recov_success] = self.RecoverDual(y,eps);
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
            x = Tuu*x(:);

        end
        
        
        function [yr,y0,success] = RecoverDual(self,yinput,eps)

            y0 = self.Ty * yinput;
             
            if (self.noReductions)
                yr = y0; success = 1;
                return
            end
                 
            if ~exist('eps','var')
                eps = 10^-4;
            end
            
            U = self.Uarry; 
    
            ComputeDelta = @(t,dir) t*dir'*self.Aorig;
            s0 = self.corig-y0'*self.Aorig;
         
            [sr,success,deltas] = solUtil.LineSearch(s0,U,self.yRedarry,self.Korig,ComputeDelta,eps);    
            if (success)
                yr = y0;
                for i=1:length(U)
                    yr = yr - deltas(i)*self.yRedarry{i};
                end
            else
                yr = [];
            end

        end
              
        function PrintStats(self)
            facialRed.PrintStats(self.K,self.Korig);
        end
          

    end

end
