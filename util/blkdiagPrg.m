classdef blkdiagPrg < reducedPrg
 
    properties
        indx
        clique
        dim
        M
    end

    methods
 
        function self = blkdiagPrg(unreducedPrg)
     
            [clique,Ar,cr,Kr,indx,M] = BuildMask(unreducedPrg.A,unreducedPrg.b,...
                                        unreducedPrg.c,unreducedPrg.K);
            
            [Ar,br,T] = CleanLinear(Ar,unreducedPrg.b);
            self@reducedPrg(Ar,br,cr,Kr);
            self.unreducedPrg = unreducedPrg;
            self.clique = clique;
           
            self.dim(1) = length(indx); 
            self.dim(2) = size(unreducedPrg.A,2);
            self.indx = indx;
            self.M = M;
            self.Ty = T;
            
        end

        
        function [xr,yr] = Recover(self,x,y)

            xr = self.RecoverPrimal(x);
            yr = self.RecoverDual(y);
            
        end
      
        function [xr] = RecoverPrimal(self,x)

            cone = coneBase(self.unreducedPrg.K);
            xr = sparse(ones(1,length(self.indx)),self.indx,x,1,cone.NumVar)';
           
        end
           
        function [yr] = RecoverDual(self,yinput,eps)

            yr = self.Ty * yinput;
             
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
