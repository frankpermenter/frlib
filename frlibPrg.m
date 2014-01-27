classdef frlibPrg

    properties
        A
        b
        c
        K
    end

    properties(GetAccess=protected)
        Z
    end

    methods

        function self = frlibPrg(A,b,c,K)

            if ~isstruct(K)
               error('Invalid input. 4th argument must be a struct')
            end


            if isempty(c)
                c = sparse(size(A,2),1);
            end

            if isempty(b)
                b = sparse(size(A,1),1);
            end

            c = c(:)';
            [A,c] = makeSymmetric(A,c,K);

            self.Z = coneHelp(A,b,c,K); 

            self.K = self.Z.K;
            self.A = self.Z.A;
            self.b = self.Z.b(:);
            self.c = self.Z.c(:);

        end

        function [x,y,info] = Solve(self)

            [A,b,T] = cleanLinear(self.A,self.b); 
            [x,y,info] = sedumi(A,b,self.c,self.K);
            y = T*y;

        end


        function [dimOut] = GetSubSpaceDim(self,Primal)

            A = cleanLinear(self.A,self.b);
            dimC = self.K.l + self.K.q + self.K.r;
            dimC = dimC + sum(self.K.s.^2/2+self.K.s/2);

            dimP = self.K.f + dimC - size(A,1);

            A = cleanLinear(self.A,self.b*0); 
            dimD = size(A,1);

            dualEqs = cleanLinear(self.A(:,1:self.K.f),self.c(1:self.K.f)); 
            dimD = dimD - size(dualEqs,2);

            if (dimP + dimD ~= dimC)
                %error('dim calc error')
            end

            if Primal
                dimOut = dimP;
            else
                dimOut = dimD;
            end

        end

        function success = CheckPrimal(self,x,eps)

            if ~exist('eps','var') || isempty(eps)
                eps = 10^-8;
            end

            [l,q,r,s]=self.GetPrimalVars(x);
            pass = [];

            pass(end+1) = all(l >= 0);

            for i=1:length(q) 
                pass(end+1) = self.CheckLor( q{i},eps);
            end

            for i=1:length(s)
                pass(end+1) = self.CheckPSD( s{i},eps);
            end

            pass(end+1) = norm(self.A*x-self.b) < eps;
            success = all(pass==1);
        end


        function success = CheckDual(self,y,eps)

            if ~exist('eps','var') || isempty(eps)
                eps = 10^-8;
            end
            
            [l,q,r,s]=self.GetDualSlack(y);
            pass = [];
            for i=1:length(q) 
                pass(end+1) = self.CheckLor( q{i},eps);
            end

            for i=1:length(s)
                pass(end+1) = self.CheckPSD( s{i},eps);
            end
            success = all(pass == 1);

        end

        function [nneg,lor,rlor,psd] = GetDualSlack(self,y)

            K = self.K;
            c = self.c;
            A = self.A;
            offset = K.f;
            indx = 1:K.f;
            freeDual = c(indx)' - y'*A(:,indx);

            indx = offset+[1:K.l];
            nneg = c(indx)' - y'*A(:,indx);
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
                psd{i} = mat(c(indx)' - y'*A(:,indx));
                offset = offset + K.s(i).^2;
            end

        end

        function [nneg,lor,rlor,psd] = GetPrimalVars(self,x)

            K = self.K;
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
                psd{i} = mat(x(indx));
                offset = offset + K.s(i).^2;
            end

        end


        function [prg] = ReducePrimal(self,method)
            method = lower(method);

            if (strcmp(method,'sdd'))
                procReduce = @facialRed.SDDPrimIter;
            end

            if (strcmp(method,'d'))
                procReduce = @facialRed.DiagPrimIter;
            end

            if (strcmp(method,'dd'))
                procReduce = @facialRed.DiagDomPrimIter;
            end

            A = self.A; b = self.b; c = self.c; K = self.K;
            T = [];
            firstPass = 1; 

            while (1)

                [success,A,c,K,Tform] = feval(procReduce,A,b,c,K);
%                [A,b] = cleanLinear(A,b);

                if success == 0
                   break;
                end

                if iscell(Tform) 
                    for i=1:length(Tform)
                        if (firstPass)
                            U{i} = Tform{i};
                        else
                            U{i} = U{i}*Tform{i};
                        end
                    end
                else
                    if (firstPass)
                         U = Tform;
                     else
                         U = U*Tform;
                     end
                end

                firstPass = 0;

                T = U; 
            end

            prg = reducedPrg(A,b,c,K,T);

        end

        function [prg] = ReduceDual(self,method)

            method = lower(method); 
            if (strcmp(method,'sdd'))
                procReduce = @facialRed.SDDDualIter;
            end

            if (strcmp(method,'d'))
                procReduce = @facialRed.DiagDualIter;
            end

            if strcmp(method,'dd')
                procReduce = @facialRed.DiagDomDualIter;
            end

            A = self.A; b = self.b; c = self.c; K = self.K;
            Deq = ones(0,size(A,1));feq=[];

            %convert these into equations on dual variables 
            if (self.K.f > 0)
                Deq = A(:,1:self.K.f)';
                feq = c(1:self.K.f);
                feq = feq(:);
                c = c(self.K.f+1:end);
                K.f = 0;
                A = A(:,self.K.f+1:end);
            end

            Scellarray= [];
            while (1)
                [success,A,c,K,Deq,feq,S] = feval(procReduce,A,c,K,Deq,feq);
                [Deq,feq] = cleanLinear(Deq,feq);
                Scellarray{end+1} = S;
                if success == 0
                    break;
                end
            end

            %convert linear constraints into free variables
            A = [Deq',A];
            c = [feq;c(:)]; 
            K.f = K.f + length(feq);
            prg = reducedPrg(A,b,c,K,[],Scellarray);

        end
    end

    methods(Static)

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
            if (min(eig(x)) > -eps) 
                pass = 1;
            else
                pass = isempty(x);
            end
        end

    end

end
