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
                display('-------------------------------------------------------------------------')
                display('frlib: No reductions found.')
                display('-------------------------------------------------------------------------')
            else
                display('-------------------------------------------------------------------------')
                display('frlib: reductions found!')
                display('-------------------------------------------------------------------------')
                display(['  Dim PSD constraint(s) (original):  ',sprintf('%d ',Korig.s)])
                display(['  Dim PSD constraint(s) (reduced):   ',sprintf('%d ',K.s)])              
                if  self.opts.useQR == 0 && size(self.A,2) < size(self.A,1)
                   display(' ')
                   display('  Warning: Reduced problem has more equations than primal variables. This may ')
                   display('  cause solver errors.  Perform reductions with opts.useQR = 1 to remove');
                   display('  linearly dependent equations.');
                end   
            end

        end

        function y = noReductions(self)
            y = ~self.faces{end}.isProper;
        end

    end

end
