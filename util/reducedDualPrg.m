classdef reducedDualPrg < frlibPrg

    properties
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
            
        end

        function [xr,x0] = RecoverRamana(self,x,eps)
            
            if (self.noReductions)
                xr = x; x0 = x; 
                return
            end
            
            
            if ~exist('eps','var')
                eps = 10^-4;
            end
            
            U = self.Uarry; V = self.Varry;
            solMap = self.solMap;
            coneOrig = ConeBase(self.Korig);
            
            s = self.cone.GetIndx('s',1);
            sBar = x(s:end);
            sHat = x(solMap.shat_vvt.s:solMap.shat_vvt.e);
            beta = x(solMap.beta_uvt.s:solMap.beta_uvt.e);
            
            
            if ~all(cellfun(@isempty,U))     
                [x,x11,x22,x21] = coneOrig.ConjBlock2by2(sBar,sHat,beta/2,U{end},V{end});
            else
                x = sBar;
            end          
            
            
            
            V = self.Varry{end}{1};
            U = self.Uarry{end}{1};       
            
            x = mat(x);
            
            xpsd = U*U'*x * U*U';
            xperp = x-xpsd;
            
            
            Gpsd{1} = mat(self.Sarry{1});
            Gperp{1} = zeros( size(Gpsd{1},1));
            for i=2:length(self.Sarry)
                G = 0;
                
               
                G =  mat(self.Sarry{i});   
             
                
                
                V = self.Varry{i-1}{1};
                U = self.Uarry{i-1}{1};
                
                
                
                Gpsd{i} =  Gpsd{i-1} + U * U' * G  * U * U';
                Gperp{i} = Gperp{i-1} + G-U * U' * G  * U * U';
                
                if (trace(Gperp{i}*U*U'))
                  error('projection failed')
                end

                
            end
            
            
            
            full([Gpsd{end},xperp;xperp',eye(4)])
            
  
            
            
        end
        
        
        
        function [xr,x0] = RecoverPrimal(self,x,eps)

            if (self.noReductions)
                xr = x; x0 = x; 
                return
            end
            
            
            if ~exist('eps','var')
                eps = 10^-4;
            end
            
            U = self.Uarry; V = self.Varry;
            solMap = self.solMap;
            coneOrig = ConeBase(self.Korig);
            
            s = self.cone.GetIndx('s',1);
            sBar = x(s:end);
            sHat = x(solMap.shat_vvt.s:solMap.shat_vvt.e);
            beta = x(solMap.beta_uvt.s:solMap.beta_uvt.e);
            
            
            xflqr = coneOrig.GetIndx('s',1);
            xflqr = x(solMap.beta_uvt.e+1:s-1);
            
            if ~all(cellfun(@isempty,U))     
                [x,x11,x22,x21] = coneOrig.ConjBlock2by2(sBar,sHat,beta/2,U{end},V{end});
            else
                x = sBar;
            end

            x(1:length(xflqr)) = xflqr;
            xr = x; 
            x0 = xr;
            fail = [];
            
            for i=length(U):-1:1

                fail(i) = 1;
                if (i >= 2)

                    Uface = U{i-1};

                else

                    Uface = [];

                end

               delta = 0.1;

               for k = 1:100
                    xr = x+self.Sarry{i}(:)'*delta;
                    feas = SolUtil.CheckInFace(xr,Uface,self.Korig,eps);
                    if (feas == 0)
                        delta = delta+1;
                    else
                        x = xr;
                        fail(i) =  0;
                        break;
                    end
                end



            end
                      
                if (any(fail))
                   error('failed') 
                end
        
        end

    end

end
