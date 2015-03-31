classdef blkdiagPrg < reducedPrg
 
    properties
        indx
        clique
        dim
    end

    methods
 
        function self = blkdiagPrg(unreducedPrg)
     
            [clique,Ar,cr,Kr,indx,M] = BuildMask(unreducedPrg.A,unreducedPrg.b,...
                                        unreducedPrg.c,unreducedPrg.K);
           
            
            [Ar,br] = CleanLinear(Ar,unreducedPrg.b);
            self@reducedPrg(Ar,br,cr,Kr);
            self.unreducedPrg = unreducedPrg;
            self.clique = clique;
            self.dim(1) =  1/2 * ( nnz(diag(M))+nnz(M));
            K = unreducedPrg.K;
            self.dim(2) =  sum( (K.s.^2 + K.s)/2   ) + K.l;
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
            xsmat = solUtil.mat(xs);
            
            %extract K.l variables
            Kl = self.unreducedPrg.K.l;
            xl = diag(xsmat(1:Kl,1:Kl));
            xsmat = xsmat(Kl+1:end,Kl+1:end);
            
            Ks = self.unreducedPrg.K.s;
            Ks = Ks(Ks~=0);
            xs = []; 
            for i=1:length(Ks)
                s = sum(Ks(1:i-1))+1;
                e = s + Ks(i)-1;
                xtemp = xsmat(s:e,s:e); xs = [xs;xtemp(:)];
            end

            xr = [xf(:);xl(:);xs(:)];
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
        

      function PrintStats(self)
            
            K = self.K;
            Korig = self.unreducedPrg.K;
   
            if  length(K.s) == length(Korig.s) & all(K.s == Korig.s)
                display('-------------------------------------------------------------------------')
                display('frlib: No reductions found.')
                display('-------------------------------------------------------------------------')
            else
                display('-------------------------------------------------------------------------')
                display('frlib: reductions found!')
                display('-------------------------------------------------------------------------')
                display(['  Dim PSD constraint(s) (original):  ',sprintf('%d ',sort(Korig.s,'descend') )])
                display(['  Dim PSD constraint(s) (reduced):   ',sprintf('%d ',sort(K.s,'descend'))])              

            end

        end
              
    end

end
