classdef reducedPrimalPrg < reducedPrg

    methods
 
        function self = reducedPrimalPrg(unreducedPrg,faces,opts)
     
            if ~exist('opts','var')
                opts.useQR = 0;
            end

            if ~isfield(opts,'useQR')
                opts.useQR = 0;
            end
            
            if ~isfield(opts,'quiet')
                opts.quiet = 0;
            end
            
            if faces{end}.isProper
                
                coneToFace = faces{end}.coneToFace;
                Ar = unreducedPrg.A*coneToFace';
                cr = unreducedPrg.c*coneToFace';
                Kr = faces{end}.K;

                [Ar,br,Ty] = CleanLinear(Ar,unreducedPrg.b,opts.useQR);
                                
            else

                Ar = unreducedPrg.A;
                br = unreducedPrg.b;
                cr = unreducedPrg.c;
                Kr = unreducedPrg.K;
                Ty  = speye(size(Ar,1));
                 
            end

            self@reducedPrg(Ar,br,cr,Kr);
            self.unreducedPrg = unreducedPrg;
            
           
            self.Ty = Ty;
            self.opts = opts;     
            self.faces = faces;
     
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

            face = self.faces{end};
            x = face.coneToFace'*x(:);

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

            ComputeDelta = @(t,redCert) t*redCert.S;
            s0 = self.unreducedPrg.c-y0'*self.unreducedPrg.A;
         
            [~,success,deltas] = solUtil.LineSearch(s0,self.faces,ComputeDelta,eps);    
            if (success)
                yr = y0;
                for i=2:length(self.faces)
                    yr = yr - deltas(i)*self.faces{i}.redCert.y;
                end
            else
                yr = [];
            end

        end
              
    end

end
