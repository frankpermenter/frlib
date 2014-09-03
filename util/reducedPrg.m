classdef reducedPrg < frlibPrg

    properties(SetAccess=protected)

        faces
        opts
        unreducedPrg
        
    end
          
    properties(Access=protected)
        
        Ty 
        y0

    end

    methods(Abstract)

         Recover(self,x,y,eps)
         RecoverPrimal(self,x)
         RecoverDual(self,yinput,eps)

    end

    methods

        function self = reducedPrg(A,b,c,K)
            self = self@frlibPrg(A,b,c,K);
        end

        function pass = VerifyRedCert(self)

            pass = 1;
            for i=1:length(self.faces)-1
                if ~self.faces{i}.InDualCone(self.faces{i+1}.redCert.S)
                    pass = 0;
                end
            end

        end
            
        
        function PrintStats(self)
            
            K = self.K;
            Korig = self.unreducedPrg.K;
   
            if  all(K.s == Korig.s)
                display(sprintf('***frlib: No reductions found.***'));
            else
                display(sprintf('***frlib: reductions found!*** Size of PSD constraint(s):'));
                display([sprintf('\t'),'Before:',sprintf('\t'),sprintf('%d ',Korig.s)])
                display([sprintf('\t'),'After:',sprintf('\t'),sprintf('%d ',K.s)])
            end

        end

        function y = noReductions(self)
            y = ~self.faces{end}.isProper;
        end

    end

end
