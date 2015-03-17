classdef solUtil

    methods(Static)
        
        function [nneg,lor,rlor,psd] = GetDualSlacks(y,A,c,K)
            K = coneBase.cleanK(K);
            nneg = []; lor = []; rlor = []; psd = [];
            c = c(:);
            offset = K.f;
            indx = 1:K.f;

            K.q = K.q(K.q>0);
            K.r = K.r(K.r>0);
            K.s = K.s(K.s>0);


            indx = offset+[1:K.l];

            if (K.l > 0)
                nneg = c(indx)' - y'*A(:,indx);
            else
                nneg = [];
            end

            offset = offset + K.l;

            for i=1:length(K.q)
                indx = offset + [1:K.q(i)];
                lor{i} = c(indx)' - y'*A(:,indx);
                offset = offset + K.q(i);
            end

            for i=1:length(K.r) 
                indx = offset + [1:K.r(i)];
                rlor{i} = c(indx)' - y'*A(:,indx);
                offset = offset + K.r(i);
            end

            for i=1:length(K.s)
                indx = offset + [1:K.s(i).^2];
                psd{i} = solUtil.mat(c(indx)' - y'*A(:,indx));
                offset = offset + K.s(i).^2;
            end
            
        end
        
        function [nneg,lor,rlor,psd]  = GetPrimalVars(x,K)
            K = coneBase.cleanK(K);
            offset = K.f;

            indx = offset+[1:K.l];
            nneg = x(indx);
            offset = offset + K.l;

            for i=1:length(K.q)
                indx = offset + [1:K.q(i)];
                lor{i} = x(indx);
                offset = offset + K.q(i);
            end

            for i=1:length(K.r) 
                indx = offset + [1:K.r(i)];
                rlor{i} = x(indx);
                offset = offset + K.r(i);
            end
            
            for i=1:length(K.s)
                indx = offset + [1:K.s(i).^2];
                psd{i} = solUtil.mat(x(indx));
                offset = offset + K.s(i).^2;
            end
            
        end
    
        function success = CheckPrimal(x,A,b,c,K,eps)

            if ~exist('eps','var') || isempty(eps)
                eps = 10^-8;
            end

            [l,q,r,s] = solUtil.GetPrimalVars(x,K);
            pass = [];

            pass(end+1) = all(l >= 0);

            for i=1:length(q) 
                pass(end+1) = solUtil.CheckLor( q{i},eps);
            end

            for i=1:length(s)
                pass(end+1) = solUtil.CheckPSD( s{i},eps);
            end

            pass(end+1) = norm(A*x-b) < eps;
            success = all(pass==1);
       
        end


        function success = CheckDual(y,A,b,c,K,eps)

            if ~exist('eps','var') || isempty(eps)
                eps = 10^-8;
            end
            
            [l,q,r,s] = solUtil.GetDualSlacks(y,A,c(:),K);
            pass = [];
            for i=1:length(q) 
                pass(end+1) = solUtil.CheckLor( q{i},eps);
            end

            for i=1:length(s)
                pass(end+1) = solUtil.CheckPSD( s{i},eps);
            end
            success = all(pass == 1);

        end
        
         
        function [dimOut] = GetSubSpaceDim(self,Primal)

            A = self.A; 
            b = self.b;
            c = self.c;
            K = self.K;
            cone = coneBase(K);
            useQR = 1;
            A = CleanLinear(A,b,useQR);
            numIndEq = size(A,1);
            
            dimCone = cone.K.f + cone.K.l + cone.K.r + cone.K.q;
            for i=1:length(cone.K.s)
                if (cone.K.s(i) > 0)
                    dimCone = dimCone  + nchoosek(cone.K.s(i) + 1,2);
                end
            end
            
            dimP = dimCone - numIndEq;
            A = CleanLinear(self.A,self.b*0,useQR); 
            numIndGen = size(A,1);

            dualEqs = CleanLinear(self.A(:,1:K.f)',self.c(1:K.f)',useQR); 
            numIndDualEq = size(dualEqs,1);
            dimD = numIndGen - numIndDualEq;

            if (dimP + dimD ~= dimCone-numIndDualEq)
                warning('dim calc error')       
            end

            if Primal
                dimOut = dimP;
            else
                dimOut = dimD;
            end

        end
                 
        function pass = CheckLor(x,eps)
            
            pass = 1;
            if length(x) > 0
                if x(1) - norm(x(2:end)) > -eps
                    pass = 1;
                else
                    pass = 0;
                end
            end
            
        end

        function pass = CheckPSD(x,eps)
            
            x = full(x);
            if (min(eig(x)) > -eps) 
                pass = 1;
            else
                pass = isempty(x);
            end
            
        end
        
                                          
        function [xr,success,deltas] = LineSearch(x,faces,ComputeDelta,eps)
            
            fail = []; xr = x; deltas = [];
            for i=length(faces):-1:2

                fail(i) = 1;
                delta = 0.0;

                %rewrite this as a gevp?
                for k = 1:100
                    deltas(i) = delta;
                    deltaX = ComputeDelta(delta,faces{i}.redCert);
                    
                    xr = x(:)+deltaX(:);
                  
                    feas = faces{i-1}.InDualCone(xr,eps);
                    if (feas == 0)
                        delta = delta+1;
                    else
                        x = xr;  
                        fail(i) =  0;
                        break;
                    end
                end
               
            end
            
            success = ~any(fail == 1);
            
        end


        function y = mat(x)
            %emulates SeDuMi's function mat()
            n = floor(sqrt(length(x)));
            y = reshape(x,n,n);
        end
      
    end
            
end
