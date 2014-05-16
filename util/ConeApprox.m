classdef ConeApprox < ConeBase


    methods

        function self = ConeApprox(K)

            self@ConeBase(K);

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
 
    end
    
end
