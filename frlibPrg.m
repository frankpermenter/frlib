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
 
            self.K = cleanK(K);
            self.A = A;

            if isempty(c)
                c = sparse(size(A,2),1);
            end

            if isempty(b)
                b = sparse(size(A,1),1);
            end

            self.b = b(:);
            self.c = c(:);
      
            self.Z = coneHelp(A,b,c,K);
    
        end

        function [x,y] = Solve(self)
            
            [A,b,T] = clean_linear(self.A,self.b); 
            [x,y] = sedumi(A,b,self.c,self.K);                
            y = T*y;

        end

        function success = CheckPrimal(self,x)
            eps = 10^-8; 
            [l,q,r,s]=self.GetPrimalVars(x);
            pass = [];

            for i=1:length(q) 
                pass(end+1) = self.CheckLor( q{i},eps);
            end

            for i=1:length(s)
                pass(end+1) = self.CheckPSD( s{i},eps);
            end

            pass(end+1) = norm(self.A*x-self.b) < eps;
            norm(self.A*x-self.b)
            success = pass;
        end


        function success = CheckDual(self,y)
           
            eps = 10^-12;
            [l,q,r,s]=GetDualSlack(y,self.A,self.c,self.K);

            for i=1:length(q) 
                pass(end+1) = self.CheckLor( q{i},eps);
            end

            for i=1:length(s)
                pass(end+1) = self.CheckPSD( s{i},eps);
            end

        end


        function [nneg,lor,rlor,psd] = GetPrimalVars(self,x)

            K = self.K;
            offset = K.f;
            indx = 1:K.f;

            indx = offset+[1:offset+K.l];
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
                procDiag = @FacialRed.SDDPrimIter;
            else

                if (strcmp(method,'d'))
                    procDiag = @FacialRed.DiagPrimIter;
                else
                    procDiag = @FacialRed.DiagDomPrimIter;
                end

            end

            A = self.A; b = self.b; c = self.c; K = self.K;
            T = [];
            firstPass = 1; 
            
            while (1)

                [success,K,A,c,Tform] = feval(procDiag,A,b,c,K);
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
           
            %T = speye(K.f);
            %for i=1:length(U)
            %    i
            %    T = blkdiag(T,UtoT(U{i}));
            %end
            prg = ReducedPrg(A,b,c,K,T);

        end

        function [prg] = ReduceDual(self,method)

            method = lower(method); 
            if (strcmp(method,'sdd'))
                procDiag = FacialRed.SDDDualIter;
            else

                if (strcmp(method,'d'))
                    procDiag = @FacialRed.DiagDualIter;
                else
                    procDiag = @FacialRed.DiagDomDualIter;
                end
       
            end

            A = self.A; b = self.b; c = self.c; K = self.K;
            Deq = ones(0,size(A,1));feq=[];

            while (1)
                [success,A,c,K,Deq,feq] = feval(procDiag,A,c,K,Deq,feq); 
                if success == 0
                    break;
                end
            end
            [Deq,feq] = clean_linear(Deq,feq);
            A = [Deq',A];
            c = c(:);
            c = [feq;c]; 
            K.f = K.f + length(feq);

            prg = ReducedPrg(A,b,c,K);
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
