classdef coneApprox < coneBase


    methods

        function self = coneApprox(K)

            self@coneBase(K);

        end

        function mats = matsFromSubMat(self,U)

            carry = []; rarry = [];varry =[];
            for i=1:length(self.K.s)
                
                if (size(U,1) ~= size(U,2))
                   error('U must be square'); 
                end
                
                if (size(U,1) > self.K.s(i))
                   continue 
                end
                
                if (size(U,1) == 2)
                    [r,c,v] = find(self.matsFromSubMat2x2_i(U,i));
                else
                    [r,c,v] = find(self.matsFromSubMat_i(U,i));
                end
                
                if ~isempty(rarry)
                    r =  r+max(rarry);
                end
                
                rarry = [rarry;r(:)];
                carry = [carry;c(:)];
                varry = [varry;v(:)];
                
            end
    
            if ~isempty(rarry)
                mats = sparse(rarry,carry,varry,max(rarry),self.NumVar);
            else
                mats = [];
            end
            
        end

        function extR = extRaysDD(self)

            v1 = self.extRaysD();
            v2 = self.matsFromSubMat([1,-1;-1,1]);
            v3 = self.matsFromSubMat([1,1;1,1]);

            extR = [v1;v2;v3;];

        end

        function extR = extRaysD(self)
            indxDiag = cell2mat(self.indxDiag);
            extR = sparse(1:length(indxDiag),indxDiag,1,length(indxDiag),self.NumVar);
        end
        
        function A =  matsFromSubMat2x2_i(self,V,num)

            if isempty(self.K.s)
               A = [];
               return
            end        

            n = self.K.s(num);
            k = size(V,2);

            if (k ~= 2)
                error('submat must be 2x2')
            end

            if n < k
                A = []; return
            end

            offset = self.Kstart.s(num)-1;
            groupings = nchoosek(1:n,k);
            numSubMats = size(groupings,1);
            sizeSubMat = size(V,1);

            blk11 = sub2ind([n,n],groupings(:,1),groupings(:,1));
            blk12 = sub2ind([n,n],groupings(:,1),groupings(:,2));
            blk22 = sub2ind([n,n],groupings(:,2),groupings(:,2));
            blk21 = sub2ind([n,n],groupings(:,2),groupings(:,1));

            Vrep = repmat(V(:)',numSubMats,1)';
            
            colIndx = [blk11,blk12,blk21,blk22]';
            colIndx = colIndx(:)';
            rowIndx = repmat(1:numSubMats,sizeSubMat*sizeSubMat,1);
            rowIndx = rowIndx(:)';

            A = sparse(rowIndx,colIndx+offset,Vrep(:));

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
