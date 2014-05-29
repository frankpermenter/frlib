classdef reducedDualPrg < frlibPrg

    properties (GetAccess=protected)
        Uarry
        Varry
        Sarry
        Karry
        solMap
        Aorig
        corig
        Korig
        noReductions
    end

    methods

        function self = reducedDualPrg(A,b,c,K,Karry,U,V,S)

            if (length(U) > 0)
                
                cone = ConeBase(K);
              
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
                Kr = ConeBase.cleanK(Karry{end});
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

            opts.removeDualEq = 1;
            self.defaultSolveOpts = opts;

        end
        
        function PrintStats(self)
            facialRed.PrintStats(self.K,self.Korig);
        end
        

    end

end
