classdef reducedPrg < frlibPrg

    properties
        T
        S
    end

    methods
    
        function self = reducedPrg(A,b,c,K,T,S)
           
            self@frlibPrg(A,b,c,K); 

            if ~exist('T','var')
                T = [];
            end
            
            if ~exist('S','var')
                S = [];
            end
            

            self.T = T;
            self.S = S;
        end

        function xr = RecoverPrimal(self,x);
            T = self.T;
            Z = self.Z;
            if (iscell(T))            
               [startPos,endPos]=Z.GetIndx('s',1);
               xr = sparse(x(1:startPos-1));
               for i=1:length(T)
                    [startPos,endPos]=Z.GetIndx('s',i);
                    xrN = T{i}*mat(x(startPos:endPos))*T{i}';
                    xr = [xr;xrN(:)];
               end 
            else
                xr = T*x;
            end
        end

    end

end
