classdef coneHelp

    properties 
        K
        A
        b
        c
        Kstart
        Kend
        indxDiag
        indxNNeg
        NumVar
    end
    
    methods(Static)

        function K = cleanK(K)

            if ~isfield(K,'f') || isempty(K.f)
                K.f = 0;
            end

            if ~isfield(K,'l') || isempty(K.l)
                K.l = 0;
            end

            if ~isfield(K,'q') || isempty(K.q)
                K.q = 0;
            end

            if ~isfield(K,'s') || isempty(K.s)
                K.s = 0;
            end
            
            if ~isfield(K,'r') || isempty(K.r)
                K.r = 0;
            end

            K.s = K.s(:)';    
            K.r = K.r(:)';
            K.q = K.q(:)';

       end

    end

    methods(Access=private)

       function self = CalcIndices(self,K)

            K = self.K;
            self.NumVar = K.f+K.l+sum(K.q)+sum(K.r)+sum(K.s.^2);

            Kstart.s = [];
            Kstart.f  = [];
            Kstart.l = [];    
            Kstart.r = [];
            Kstart.q  = [];

            Kend.s = [];
            Kend.f  = [];
            Kend.l = [];    
            Kend.r = [];
            Kend.q  = [];

            if (K.f)
                Kstart.f = 1;
                Kend.f = Kstart.f + K.f - 1;
            end

            if (K.l)
                Kstart.l =  K.f + 1;
                Kend.l = Kstart.l + K.l - 1;
            end

            if (any(K.q))
                offset = K.f + K.l; 
                temp = [0,cumsum(K.q)]+1 + offset;
                Kstart.q = temp(1:end-1);
                Kend.q = Kstart.q + K.q - 1;
            end

            if (any(K.r))
                offset = K.f + K.l + sum(K.q);
                temp = [0,cumsum(K.r)]+1 + offset;
                Kstart.r = temp(1:end-1);
                Kend.r = Kstart.r + K.r - 1;
            end

            indxDiag=[];
            if (any(K.s))
                offset = K.f + K.l + sum(K.q) + sum(K.r); 
                temp = [0,cumsum(K.s.^2)]+1;
                Kstart.s = temp(1:end-1) + offset;
                Kend.s = Kstart.s + K.s.^2 - 1;
               

                for i=1:length(Kstart.s)
                    indxDiag{i} = [Kstart.s(i):K.s(i)+1:Kend.s(i)];
                end
            end
        
            self.Kstart = Kstart;
            self.Kend = Kend; 
            self.indxDiag = indxDiag;
            self.indxNNeg = [Kstart.l:Kend.l,Kstart.q,Kstart.r,Kstart.r+1,cell2mat(indxDiag)];

            end

    end
    
    methods


    function [A,b] = innerPrdSubMatDualRetired(self,V)
      
        for i=1:length(self.K.s) 
        
            n = self.K.s(i);
            k = size(V,1);
            offset = self.Kstart.s(i)-1;

            if n < k 
                A = []; return
            end

            subMat = nchoosek(1:n,k); 
            delta=diff(subMat,1,2);

            %each row contains the indices of a submatrix
            indices=[self.indxDiag{i}(subMat(:,1))',self.indxDiag{i}(subMat(:,1))'+delta,...
                 self.indxDiag{i}(subMat(:,2))'-delta,self.indxDiag{i}(subMat(:,2))'];
            M = size(self.A,1);            
            Vrep = sparse(repmat(V(:),1,M)');
            for j=1:size(indices,1)
                A(:,j) = sum(self.A(:,indices(j,:)).*Vrep,2); 
                b(j,1) = sum(self.c(indices(j,:)).*V(:));
            end     

            A = A';
        end
        
    end

    function [A,b] = DualIneq(self,U) 
        c = self.c(:);
        A = U*self.A';
        b = U*c; 
    end
 
    function [A,b,v] = ineqDiagDomDual(self)
       
        [Aineq1,b1,v1] = self.innerPrdSubMatDual(1);
        [Aineq2,b2,v2] = self.innerPrdSubMatDual([1,-1;-1,1]);
        [Aineq3,b3,v3] = self.innerPrdSubMatDual([1,1;1,1]);
        
        A = [Aineq1;Aineq2;Aineq3];
        b = [b1;b2;b3]; 
        v = [v1;v2;v3]; 

    end

    function [A,b,E] = innerPrdSubMatDual(self,V)
        b = []; E = []; A = []; 
        for i=1:length(self.K.s)
          Et = self.innerPrdSubMat(V,i);  
          [At,bt] = self.DualIneq(Et); 
          E = [E;Et];
          A = [A;At];
          b = [b;bt];
        end

    end

    function [A,b] = innerPrdDiagDual(self)
        A = self.A(:,cell2mat(self.indxDiag))';
        b = self.c(cell2mat(self.indxDiag));
    end

    function A = getA(self,cone,num)
        [startPos,endPos]=self.GetIndx(cone,num);
        A = self.A(:,startPos:endPos);
    end

    function self = coneHelp(A,b,c,K) 
        self.K = self.cleanK(K);
        self.A = A;
        
        if  isempty(b)
            self.b = zeros(size(A,1),1);
        else
            self.b = b(:);
        end

        if  isempty(c)
            self.c = zeros(size(A,2),1);
        else
            self.c = c(:);
        end
        
        self = self.CalcIndices(self.K); 
        if self.NumVar ~= size(self.A,2)
            error(['Cone sizes do not match columns of A. Num vars: ',num2str(self.NumVar),' Num Col: ', num2str(size(self.A,2))]);
        end
    end
 

    function A = innerPrd(self,U,cone,num)

        startPos = getfield(self.Kstart,cone);
        startPos = startPos(num);
        endPos = getfield(self.Kend,cone);
        endPos = endPos(num);
        
        A = sparse(1,startPos:endPos,U(:),1,self.NumVar);

    end

    function [startPos,endPos]= GetIndx(self,cone,num)
        startPos = getfield(self.Kstart,cone);
        startPos = startPos(num);
        endPos = getfield(self.Kend,cone);
        endPos = endPos(num);
    end 

    function UAU = ConjA(self,U,num)  

        [startPos,endPos]=self.GetIndx('s',num);
        A = self.A(:,startPos:endPos);

    %   [N,M] = size(U);
    %   UAU = spalloc(size(A,1),M*M,nnz(A));
        T = UtoT(U');
        UAU = (T * A')';
%        for i=1:size(A,1)
%            temp = U'*mat(A(i,:))*U;
%            UAU(i,:) = temp(:)';
%        end
    end


    function UCU = ConjC(self,U,num)  

        [startPos,endPos]=self.GetIndx('s',num);
        c = self.c(startPos:endPos);
        temp = U'*mat(c)*U;
        UCU = temp(:)';

    end

    function [U,V] = computeU(self,A,num)

        [startPos,endPos]=self.GetIndx('s',num);
        A = sum(A(:,startPos:endPos),1);

        [U,V] = nullqr(mat(A)'); 
    end


    function A = innerPrdSubMat(self,V,num)
     
        n = self.K.s(num);
        k = size(V,1);
        offset = self.Kstart.s(num)-1;

        if n < k 
            A = []; return
        end

        posArray = nchoosek(1:n,k); 
        numSubMat = size(posArray,1);
        
        nnzA = numSubMat*nnz(V);

        cnt = 1;
        cols = zeros(nnzA,1);
        rows = zeros(nnzA,1);
        val = zeros(nnzA,1);
        for j=1:numSubMat
            
            %get row/col defining the submat
            pos = posArray(j,:);   
            
            %get first entry of submat
            indxStart = n*(pos(1)-1)+pos(1);

            %how far apart are defining rows 
            deltaRow = [0,cumsum(diff(pos))];

            %indices of the first column of sub mat
            firstCol = indxStart(1) + deltaRow;
      
            rowPos = j*ones(length(k),1);
            for p = 1:length(pos) 
                cols(cnt:cnt+k-1) = firstCol + n*deltaRow(p)+offset;
                rows(cnt:cnt+k-1) = rowPos;
                val(cnt:cnt+k-1) = V(:,p);
                cnt = cnt+k;
            end
        end

        A = sparse(rows,cols,val,numSubMat,self.NumVar);
    end

    function A = innerPrdDiag(self)

        NumIneq = length(self.indxNNeg);
        A = sparse(1:NumIneq,self.indxNNeg,-1,NumIneq,self.NumVar);

    end


    function A = innerPrdDD(self)

        A = self.innerPrdDiag();

        for i = 1:length(self.K.s)
      
            V = ones(2); 
            A = [A;-self.innerPrdSubMat(V,i)];
            V(2,1) = -1*V(2,1);
            V(1,2) = -1*V(1,2);
            A = [A;-self.innerPrdSubMat(V,i)];

        end

    end


    function A = innerPrdFW(self,V)
        
        K = self.K;
        %inner product all semidefinite variables 
        %with matrix V
        offset = K.f + K.l + sum(K.q) + sum(K.r);
        NumVars = K.f + K.l + sum(K.q) + sum(K.r) + sum(K.s.^2);
        A = 0*speye(1,NumVars);

        cnt = 1;
        k = size(V,1);

        for i=1:length(K.s)
            
            n = K.s(i);
            
            if (n < k)
                continue
            end
            posArray = nchoosek(1:n,k); 
            for j=1:size(posArray,1)
                pos = posArray(j,:);   
                indxDiag = n*(pos-1)+pos;
                delta = [0,cumsum(diff(pos))];
                colOne = indxDiag(1) + delta;
         
                for k = 1:length(pos) 
                    col = colOne + n*delta(k)+offset;
                    A(cnt,col) = V(:,k);
                end

                cnt = cnt + 1;
            end

            offset = offset + K.s(i).^2;

        end
        [row,col] = find(A);
        row = unique(row);
        A = A(row,:);

    end

    function [Aeq] = AtimesV(self,V,num)

        [startPos,endPos]=self.GetIndx('s',num);
        A = self.A(:,startPos:endPos);
        for i=1:size(A,1)
            Aeq(:,i) = mat(A(i,:))*V; 
        end
        
    end

    function y = AdotX(self,X,num) 

        [startPos,endPos]=self.GetIndx('s',num);
        A = self.A(:,startPos:endPos);
        m = size(A,1);
        X = X(:)';  
        Xrep = repmat(X,m,1);
        y = sum(Xrep.*A,2); 

    end

    function y = AdjA(self,u,num) 

        [startPos,endPos]=self.GetIndx('s',num);
        A = self.A(:,startPos:endPos);
        m = size(A,1);
        y = 0; 
        for i=1:m
            y = y + A(i,:)*u(i);
        end
    end


    function y = CdotX(self,X,num) 

        [startPos,endPos]=self.GetIndx('s',num);
        c = self.c(startPos:endPos);
        for i=1:size(X,1)
            y(i,1) = sum(c(:)'.*X(i,:));
        end

    end


    function [beq] = CtimesV(self,V,num)

        [startPos,endPos]=self.GetIndx('s',num);
        C = self.c(startPos:endPos);
        beq = mat(C)*V;

    end


    %Find variables in cone that vanish if others vanish
    function [indxZero,Knew] = FindMustVanish(self,indx)
      
        indx = unique(indx);
        indx = indx(:);
        indxZero = indx';
        Knew = self.K;

        if isempty(indx)
            return   
        end

        if self.Kstart.f
            indxF = indx(indx <= self.Kend.f);
            Knew.f = Knew.f - length(indxF);
        end

        if self.Kstart.l
            indxL = indx(indx >= self.Kstart.l & indx <= self.Kend.l); 
            Knew.l = Knew.l - length(indxL);
        end

        if self.Kstart.q
        
            indxQ = indx(indx >= self.Kstart.q(1) & indx <= self.Kend.q(end));
            for i=indxQ
                cnstN = max(find(K.start.q <= i));
                Knew.q(cnstN) = Knew.q(cnstN) - 1;
            end

            %zero out other variables if first
            [~,cnstN] = intersect(self.Kstart.q,indxQ);       
            Knew.q(cnstN) = 0;
            for i = cnstN
                indxZero = [indxZero,self.Kstart.q(i)+1:self.Kend.q(i)];
            end 
        end

        %Rotated lorentz constraints
        %remove variable from constraint    
        if self.Kstart.r

            indxR = indx(indx >= self.Kstart.r(1) & indx <= self.Kend.r(end)); 
            for i=indxR
                cnstN = max(find(self.Kstart.r <= i));
                Knew.r(cnstN) = Knew.r(cnstN) - 1;
            end

            %zero out rotated constraints 
            [~,cnstR1] = intersect(self.Kstart.r,indxR);      
            [~,cnstR2] = intersect(self.Kstart.r+1,indxR);   
            cnstR = union(cnstR1,cnstR2);
            cnstR12 = intersect(cnstR2,cnstR1);
            Knew.r(cnstR) = 1;
            Knew.r(cnstR12) = 0;
            
            for i = cnstR
                indxZero = [indxZero,self.Kstart.r(i)+2:self.Kend.r(i)];
            end 

        end

        if self.Kstart.s
            
            indxS = indx(indx >= self.Kstart.s(1) & indx <= self.Kend.s(end));
            
            indxD = intersect(cell2mat(self.indxDiag),indxS);
           
            for indx=indxD
                  
               cnstN = max(find(self.Kstart.s <= indx));
               offset = self.Kstart.s(cnstN)-1;
               N = self.K.s(cnstN);

               %find variables on same row and col
               row = floor( (indx-offset-1)/N)+1;
               indxColStart = N*(row-1) + 1;
               indxCol = indxColStart:indxColStart+N-1;
               indxRow = row:N:N*N;
               indxRC = union(indxCol,indxRow)+offset;
               %indxRC = setdiff(indxRC,indx);
               %decrement dimension of LMI by 1
               Knew.s(cnstN) =  Knew.s(cnstN) - 1; 
               indxZero = [indxZero,indxRC];
                
            end
            
        end


    end

end



end
