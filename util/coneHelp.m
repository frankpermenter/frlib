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

        if self.NumVar ~= length(c)
            error(['Number of variables do not match length of c. Num vars: ',num2str(self.NumVar),', Length c: ', num2str(length(c))]);
        end

        if length(self.b) ~= size(A,1)
            error(['Number of rows of A do not match length of b. Num rows: ',num2str(size(A,1)),', Length b: ', num2str(length(b))]);
        end



    end

    function y = anyConicVars(self)
       y = length(self.indxNNeg) > 0;
    end


    function indx =  lowerTriIndx(self)

        startIndx = self.GetIndx('s',1);

        if isempty(startIndx)
            indx = 1:size(self.A,2);
            return
        end

        indx = [1:startIndx-1];

        for i=1:length(self.K.s)
            indxDiag = self.indxDiag{i};
            for j=1:self.K.s(i)
               numIndx = self.K.s(i)-j;
               indx = [indx,indxDiag(j):indxDiag(j)+numIndx];
            end
        end

    end

    function indx  =  upperTriIndx(self)

        [startIndx] = self.GetIndx('s',1);
        if isempty(startIndx)
            indx = 1:size(self.A,2);
            return
        end

        indx = [1:startIndx-1];

        for i=1:length(self.K.s)
            [~,endIndx] = self.GetIndx('s',i);
            indxDiag = self.indxDiag{i};
            N = self.K.s(i);
            for j=1:self.K.s(i)
               indx = [indx,indxDiag(j):N:endIndx];
            end
        end

    end

    function A =  symmetrize(self,A)

        indxL = self.lowerTriIndx();
        indxU = self.upperTriIndx();

        As = (A(:,indxL) + A(:,indxU) )/2;

        A(:,indxL) = As;
        A(:,indxU) = As;

    end

    function A =  lowerTri(self,A)

        indx = self.lowerTriIndx();
        A = A(:,indx);

    end

    function A =  upperTri(self,A)

        indx = self.upperTriIndx();
        A = A(:,indx);

    end

    function A = flqrCols(self,A)
       cols = 1:max([self.Kend.f;self.Kend.l;self.Kend.r(:);self.Kend.q(:);0]);
       A = A(:,cols);
    end


    function mats = matsFromSubMat(self,U)
        mats =[];
        for i=1:length(self.K.s)
            [mats] = [mats;self.matsFromSubMat_i(U,i)];
        end
    end

    function extR = extRaysDD(self)

        v1 = sparse(1:length(self.indxNNeg),self.indxNNeg,1,length(self.indxNNeg),self.NumVar);
        v2 = sparse(0,self.NumVar);
        for i=1:length(self.K.q)
            N = self.K.q;
            if (N > 0)
                vq1 = [ones(N-1,1),speye(N-1)];
                vq2 = [ones(N-1,1),-speye(N-1)];
                [s,e] = self.GetIndx('q',i);
                v2(end+1:end+2*N-2,s:e) = [vq1;vq2];
            end
        end


        v1 = self.matsFromSubMat([1]);
        v2 = [];
        v3 = self.matsFromSubMat([1,-1;-1,1]);
        v4 = self.matsFromSubMat([1,1;1,1]);

        extR = [v1;v2;v3;v4];

    end

    function extR = extRaysD(self)
        %extR = sparse(1:length(self.indxNNeg),self.indxNNeg,1,length(self.indxNNeg),self.NumVar);
        indxDiag = cell2mat(self.indxDiag);
        extR = sparse(1:length(indxDiag),indxDiag,1,length(indxDiag),self.NumVar);
        %extR = self.matsFromSubMat([1]);
    end


    function [startPos,endPos]= GetIndx(self,cone,num)

        startPos = getfield(self.Kstart,cone);
        endPos = [];

        if (startPos)
            startPos = startPos(num);
            endPos = getfield(self.Kend,cone);
            endPos = endPos(num);
        end

    end

    function UAU = ConjA(self,U,num)

        [startPos,endPos]=self.GetIndx('s',num);
        A = self.A(:,startPos:endPos);

        T = UtoT(U');
        UAU = (T * A')';

    end

    function UCU = ConjC(self,U,num)

        [startPos,endPos]=self.GetIndx('s',num);
        c = self.c(startPos:endPos);
        temp = U'*mat(c)*U;
        UCU = temp(:)';

    end

    function A = matsFromSubMat_i(self,V,num)

        if isempty(self.K.s)
           A = [];
           return
        end

        n = self.K.s(num);
        k = size(V,1);

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
        offset = self.Kstart.s(num)-1;

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
            %loop across columns
            for p = 1:length(pos)
                cols(cnt:cnt+k-1) = firstCol + n*deltaRow(p)+offset;
                rows(cnt:cnt+k-1) = rowPos;
                val(cnt:cnt+k-1) = V(:,p);
                cnt = cnt+k;
            end
        end

        A = sparse(rows,cols,val,numSubMat,self.NumVar);

    end


    function [Aeq] = AtimesV(self,V,num)

        [startPos,endPos]=self.GetIndx('s',num);
        A = self.A(:,startPos:endPos);
        for i=1:size(A,1)
            Aeq(:,i) = mat(A(i,:))*V;
        end

    end

    function [beq] = CtimesV(self,V,num)

        [startPos,endPos]=self.GetIndx('s',num);
        C = self.c(startPos:endPos);
        beq = mat(C)*V;

    end

    function y = AdotX(self,X,num)

        [startPos,endPos]=self.GetIndx('s',num);
        A = self.A(:,startPos:endPos);
        m = size(A,1);
        X = X(:);
        y = A*X;

    end

    function y = CdotX(self,X,num)

        [startPos,endPos]=self.GetIndx('s',num);
        C = self.c(startPos:endPos);
        X = X(:);
        y = C(:)'*X;

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
