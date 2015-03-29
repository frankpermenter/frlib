classdef blkdiagPrg < reducedPrg
 
    properties
        indx
        clique
    end

    methods
 
        function self = blkdiagPrg(unreducedPrg)
     
            [clique,Ar,cr,Kr,indx] = BuildMask(unreducedPrg.A,unreducedPrg.b,...
                                        unreducedPrg.c,unreducedPrg.K);
                                                             
            self@reducedPrg(Ar,unreducedPrg.b,cr,Kr);
            self.unreducedPrg = unreducedPrg;
            self.clique = clique;
            
            self.indx = indx;
           
     
        end

        
        function [xr,yr,dual_recov_success] = Recover(self,x,y,eps)

     

            if ~exist('eps','var')
                eps = 10^-4;
            end
           
            xr = self.RecoverPrimal(x);
            [yr,~,dual_recov_success] = self.RecoverDual(y,eps);
            if (dual_recov_success ~= 1)
                yr = [];
            end
  
        end
      
        function [xr] = RecoverPrimal(self,x)

            xf = x(1:self.K.f);
            xs(self.indx) = x(self.K.f+1:end);
           
            xsmat = mat(xs);
            Ks = [self.unreducedPrg.K.l,self.unreducedPrg.K.s];
            Ks = Ks(Ks~=0);
            xs = [];
            for i=1:length(Ks)
                s = sum(Ks(1:i-1))+1;
                e = s + Ks(i)-1;
                xtemp = xsmat(s:e,s:e); xs = [xs;xtemp(:)];
            end

            xr = [xf(:);xs(:)];
        end
           
        function [yr,y0,success] = RecoverDual(self,yinput,eps)

            y0 = self.Ty * yinput;
             
        end
        
        function SolveDualLP(self)
            
            [bnd,bndlp]= checkoptimality(self.unreducedPrg.A,self.unreducedPrg.b,self.unreducedPrg.c,self.unreducedPrg.K);
            bnd - bndlp
            
        end
        
        function SolvePrimalLP(self)
            
            [M]= buildip2(self.unreducedPrg.A,self.unreducedPrg.b,self.unreducedPrg.c,self.unreducedPrg.K);
            sizes = cellfun(@(x) size(x,1),self.clique); 
            cliM = sum(sizes.*sizes);
            optM = sum(sum(M));
            optM-cliM
            
        end        
        
              
    end

end
