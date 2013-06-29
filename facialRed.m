classdef facialRed
 
    methods(Static)

        function [success,K,A,c,T] = SDDPrimIter(A,b,c,K)
            
            T = [];
            Z = coneHelp(A,b,c,K);
            prg = spotsosprg;
            m = size(Z.A,1);
            tr = 0;
            
            [prg,u] = prg.newFree(m);
            prg = prg.withEqs(u'*Z.b);
            
            for i=1:length(K.s)  
                if (K.s(i) > 1 )
                    [prg,Stemp] = prg.newSDD(K.s(i));   
                    S2 = mat(Z.AdjA(u,i));
                    S{i} = S2;
                    prg = prg.withEqs(Stemp-S2);
                    tr = tr + trace(Stemp);
                end
            end
            prg = prg.withEqs(tr-1);
            try
                sol = prg.optimize();
                success = sol.info.pinf == 0 & sol.info.numerr == 0;
            catch
                success = 0;
            end
            if success

                A = A(:,1:Z.GetIndx('s',1)-1); 
                c = c(1:Z.GetIndx('s',1)-1)'; 
                c = c(:)';
                for i = 1:length(K.s) 
                
                    V = double(sol.eval(S{i}));
                    [U,sd] = nullqr(V);
                    if (K.s(i) > 0)
                        K.s(i) = size(U,2);
                        A = [A,Z.ConjA(U,i)];
                        c = [c,Z.ConjC(U,i)];
                    end
                    T{i} = U;

                end
           end 
        end

        function [success,A,c,K,Deq,feq] = SDDDualIter(A,c,K,Deq,feq)

            Z = coneHelp(A,[],c,K);
            prg = spotsosprg;

            tr = 0;
            for i=1:length(K.s)  
                [prg,Stemp] = prg.newSDD(K.s(i));
                S{i} = Stemp;
                prg = prg.withEqs(Z.AdotX(Stemp,i));
                tr = tr + trace(Stemp);
            end

            prg = prg.withEqs(tr-1);
            try
                sol = prg.optimize();
                success = sol.info.pinf == 0;
            catch
                success = 0;
            end

            if success
                A = A(:,1:Z.GetIndx('s',1)-1); 
                c = c(1:Z.GetIndx('s',1)-1); 
                for i = 1:length(K.s) 
                    V = double(sol.eval(S{i}));
                    U = nullqr(V);
                    if (K.s(i) > 0)
                        K.s(i) = size(U,2);
                        A = [A,Z.ConjA(U,i)];
                        c = [c,Z.ConjC(U,i)];
                        for j=1:size(V,2)
                            Deq = [Deq;Z.AtimesV(V(:,j),i)];
                            feq = [feq;Z.CtimesV(V(:,j),i)];
                        end
                    end
                end
            end 

        end

        function [success,K,A,c,T] = DiagPrimIter(A,b,c,K)

            Z = coneHelp(A,[],[],K);
            A_ineq = Z.innerPrdDiag(); 
            T = []; 

           [activeIneq,infeas]= lpact(A_ineq,0,A,b);
            if (infeas)
                success = 0;
               % K.s = -1;
                return 
            end

            activeIneq = find(activeIneq);
            [~,zeroVar] = find(A_ineq(activeIneq,:));    
            
            [zeroVar,K] = Z.FindMustVanish(zeroVar);   
            T=speye(Z.NumVar);
            T(zeroVar,:) = 0; 
            indxKeep = setdiff(1:Z.NumVar,zeroVar);

            T = T(:,indxKeep);
            c = c(indxKeep);
            A = A*T;

            if ~isempty(zeroVar)
                success = 1;
            else
                success = 0;
            end

        end

        function [success,K,A,c,T] = DiagDomPrimIter(A,b,c,K)


            %the extreme rays of DD include the extreme rays of D. Try 
            %reducing with D first since it is cheaper 
            %[success,K,A,c,T] = FacialRed.DiagPrimIter(A,b,c,K);
            %if (success == 0 )

                T = [];
                Z = coneHelp(A,[],c,K);
                A_ineq = Z.innerPrdDD();
                [activeIneq,infeas]= lpact(A_ineq,0,A,b);
                if (infeas)
                    success = 0;
               %     K.s = -1;
                    return 
                end

                activeIneq = find(activeIneq);

                A = A(:,1:Z.GetIndx('s',1)-1); 
                c = c(1:Z.GetIndx('s',1)-1)'; 
                c = c(:)';
                for i = 1:length(K.s) 
                    [U,V] = Z.computeU(A_ineq(activeIneq,:),i);

                    if (K.s(i) > 0)
                        K.s(i) = size(U,2);
                        A = [A,Z.ConjA(U,i)];
                        c = [c,Z.ConjC(U,i)];
                    end
                    T{i} = U;

                end
                success = length(activeIneq) > 0;
            %end

        end
        
        
        function [success,A,c,K,Deq,feq] = DiagDomDualIter(A,c,K,Deq,feq)
            %the extreme rays of DD include the extreme rays of D. Try 
            %reducing with D first since it is cheaper 
            [success,A,c,K,Deq,feq] = facialRed.DiagDualIter(A,c,K,Deq,feq);

            if (success == 0 && all(K.s > 1))

                Z = coneHelp(A,[],c,K);

                [A_ineq,b_ineq,V] = Z.ineqDiagDomDual();                                                   

                activeIneq = find(lpact(A_ineq,b_ineq,Deq,feq));
                success = length(activeIneq) > 0;

                if success  

                    Deq = [Deq;A_ineq(activeIneq,:)];
                    feq = [feq;b_ineq(activeIneq,:)];
                   
                    A = A(:,1:Z.GetIndx('s',1)-1); 
                    c = c(1:Z.GetIndx('s',1)-1); 
                    
                    for i = 1:length(K.s) 
                        [U,V] = Z.computeU(V(activeIneq,:),i);
                        if (K.s(i) > 0)
                            K.s(i) = size(U,2);
                            A = [A,Z.ConjA(U,i)];
                            c = [c,Z.ConjC(U,i)];
                            for j=1:size(V,2)
                                Deq = [Deq;Z.AtimesV(V(:,j),i)];
                                feq = [feq;Z.CtimesV(V(:,j),i)];
                            end
                        end
                    end

                end
            end 

        end


        function [success,A,c,K,Deq,feq] = DiagDualIter(A,c,K,Deq,feq)
         
            Z = coneHelp(A,[],c,K);
            [A_ineq,b_ineq] = Z.innerPrdDiagDual(); 
           
            activeIneq = find( lpact(A_ineq,b_ineq,Deq,feq) );
            Deq = [Deq;A_ineq(activeIneq,:)];
            feq = [feq;b_ineq(activeIneq,:)];

            %find columns of A to delete 
            zeroCol = Z.indxNNeg(activeIneq);
            [zeroCol,K] = Z.FindMustVanish(zeroCol);   
          
            %Q V = 0 
            Deq2 = A(:,zeroCol)'; 
            feq2 = zeros(length(zeroCol),1);
            Deq = [Deq;Deq2];
            feq = [feq;feq2];

            if ~isempty(zeroCol)
                success = 1;
                indxKeep = setdiff( 1:size(A,2),zeroCol);
                A = A(:,indxKeep);
                c = c(indxKeep,1);              
             else
                success = 0;
            end
        end


        function [varRmv,Knew,Anew,c,T] = doIter(K,A,b,c,DD)

            varRmv=[];Knew=[];Anew=[];T=[];col=[];

            NumVars = K.f + K.l + sum(K.q) + sum(K.r) + sum(K.s.^2);
            
            %get variables that are non-negative if conic
            %constraints in K hold
            nneg = sedumiNonNeg(K);
            A_ineq = sparse(length(nneg),NumVars);
            for i=1:length(nneg)
                A_ineq(i,nneg(i)) = -1;
            end

            if (DD == 1)
                V = ones(2);
                A_ineq1 = -sedumiNbyN(K,V);
                V(2,1) = -1*V(2,1);
                V(1,2) = -1*V(1,2);

                A_ineq2 = -sedumiNbyN(K,V);
                A_ineq1 = [A_ineq1;A_ineq2];
                A_ineq = [A_ineq;A_ineq1];
            else
                A_ineq1 = [];
            end

            b_ineq = zeros(size(A_ineq,1),1);
            y = lpact(A_ineq,b_ineq,A,b);
                         
            N = size(A_ineq,1)-size(A_ineq1,1);
                
            if (1)%DD == 0)
               % Knew = K;
                %save test.mat Knew A_ineq A_ineq1 A_ineq2 y A b
                %zeroVar2 = removeDiag(K,A_ineq1,y(N+1:end));
            else
                zeroVar2 = [];
            end

            zeroVar2 = [];
            %active constraints
            activeCnst = find(y(1:N));

            BuildU(A_ineq(activeCnst,:));


            [~,zeroVar1] = find(A_ineq(activeCnst,:));

            indxArray = [zeroVar1];
            %other vars that must vanish
            [varRmv,Knew] = FindMustVanish(K,zeroVar1);
            %build map: oldVar = T*newvar
            T=speye(NumVars);
            T(varRmv,:) = 0;
            
            indxKeep = setdiff(1:NumVars,varRmv);
            T = T(:,indxKeep);
            %move lorentz/rotated lorentz constraints
            %with one variable to K.l
            convertToLinear = find(Knew.q == 1);
            offsetQ = Knew.f + Knew.l;
            
            T2 = speye(length(indxKeep));
            for i = 1:length(convertToLinear)
                indxPaste = offsetQ + sum(Knew.q(1:convertToLinear(i)-1))+1;
                indxCut = Knew.l + 1;
                T2 = CutAndPaste(T2,indxCut,indxPaste);
                Knew.l = Knew.l + 1;
                Knew.q(convertToLinear(i)) = 0;
            end

            convertToLinear = find(Knew.r == 1);
            offsetR = Knew.f + Knew.l+sum(Knew.q);
            for i = 1:length(convertToLinear)
                indxPaste = offsetR + sum(Knew.r(1:convertToLinear(i)-1))+1;
                indxCut = Knew.l + 1;
                T2 = CutAndPaste(T2,indxCut,indxPaste);
                Knew.l = Knew.l + 1;
                Knew.r(convertToLinear(i)) = 0;
            end
            T = T*T2;
            Anew = A * T;
            c = c(indxKeep);

            
            
        end


        function T = CutAndPaste(T,indxCut,indxPaste)


            if (indxCut == indxPaste)
                return
            end

            if (indxPaste > indxCut)
                T = [T(1:indxCut-1,:);T(indxCut+1:indxPaste,:);T(indxCut,:);T(indxPaste+1:end,:)];
            else
                T = [T(1:indxPaste,:);T(indxCut,:);T(indxPaste+1:indxCut-1,:);T(indxCut+1:end,:)];
            end

        end

    end
end




