classdef coneBase

    properties
        K
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

        function self = coneBase(K)

            self.K = self.cleanK(K);
            self = self.CalcIndices(self.K);

        end
        
        function y = AnyConicVars(self)
           y = length(self.indxNNeg) > 0;
        end

        function y = NumFlqrVars(self)
           [s,e] = self.flqrIndx();
           y = max(0,e-s+1);
        end


        function indx =  LowerTriIndx(self)

            startIndx = self.GetIndx('s',1);

            if isempty(startIndx)
                indx = 1:self.NumVar;
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

        function indx  =  UpperTriIndx(self)

            [startIndx] = self.GetIndx('s',1);
            if isempty(startIndx)
                indx = 1:self.NumVar;
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


        function A =  Symmetrize(self,A)

            indxL = self.LowerTriIndx();
            indxU = self.UpperTriIndx();

            As = (A(:,indxL) + A(:,indxU) )/2;

            A(:,indxL) = As;
            A(:,indxU) = As;

        end

        function A =  LowerTri(self,A)

            indx = self.LowerTriIndx();
            A = A(:,indx);

        end

        function A =  UpperTri(self,A)

            indx = self.UpperTriIndx();
            A = A(:,indx);

        end

        
        function A = InitSymmetric(self,vals)
                       
            A = sparse(size(vals,1), self.NumVar);
            A(:,self.LowerTriIndx()) = vals;
            A(:,self.UpperTriIndx()) = vals;
            
        end
        
        
        function A = Desymmetrize(self,A)

            indxDiag = cell2mat(self.indxDiag);
            temp = A;
            temp(:,indxDiag) =  temp(:,indxDiag)/2;

            indxRescale = self.LowerTriIndx();
            indxRescale = indxRescale( indxRescale >= indxDiag(1));
            temp(:,indxRescale) = temp(:,indxRescale)*2;
            temp = self.LowerTri(temp);

            A = sparse( size(A,1), size(A,2));
            A(:,self.LowerTriIndx()) = temp;

        end

        function A = flqrCols(self,A)

           cols = 1:max([self.Kend.f;self.Kend.l;self.Kend.r(:);self.Kend.q(:);0]);
           A = A(:,cols);

        end
        
        
        function [s,e] = flqrIndx(self)
                   
           s = 1;
           e = max([self.Kend.f;self.Kend.l;self.Kend.r(:);self.Kend.q(:);0]);
  
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


        function [T,s] = ConjByU(self,U,num)

            [T,s] = PrePostMultByUandVt(U,U)

        end

        function [T,s] = PrePostMultByUandVt(self,U,V,num)

            if exist('num','var')  
                T = BuildConjMap(U,V);
            else
                T = speye(self.GetIndx('s',1)-1);

                for i=1:length(U)
                    T = blkdiag(T,BuildPrePostMultMap(U{i},V{i}));
                    s(i) = size(U{i},2);
                end  
            end

        end


        function [y] = ColIndx(self,num,colIndx)

            [startPos,~] = self.GetIndx('s',num);
            n = self.K.s(num);
            s = (colIndx-1)*n+startPos;
            e = s+n-1;
            y = [s,e];

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


        function T = BuildMultMap(self,U1,U2,transposeU)
            %if tranpose U == 0): Computes map (Tq)(:) = (U1 q U2') (:)
            %else: Computes map (Tq)(:) = (U1' q U2) (:) 
            
            if (~exist('transposeU','var'))
               transposeU = 1;
            end
            
            if (transposeU)
               dimForSizeCalc = 2;
               dimForErrorCheck = 1;
            else
               dimForSizeCalc = 1; 
               dimForErrorCheck = 2;
            end
            

            if ~all(cellfun(@(x)  size(x,dimForErrorCheck),U2) == self.K.s)
               error('size of inputs do not match cone sizes') 
            end
            
            rarry = [];carry=[];varry=[];
            [s,e] = self.flqrIndx();

            numNonPSDVar = e-s+1;  
            numPSDVarAfterConj_1 = cellfun(@(x)  size(x,dimForSizeCalc),U1);   
            numPSDVarAfterConj_2 = cellfun(@(x)  size(x,dimForSizeCalc),U2);   
            numPSDVarAfterConj = sum(numPSDVarAfterConj_2.*numPSDVarAfterConj_1);
            
            %don't do anything to non PSD vars
            rarry = [s:e]';
            carry = [s:e]';
            varry = ones(size(rarry,1),1);
            

            indx = length(s:e)+1;
            for i=1:length(U1)
                s = self.GetIndx('s',i);
                if any(size(U1{i}) > 0) && any(size(isempty(U2)) > 0)
                    
                    if (transposeU == 1)
                        temp = BuildConjMap(U1{i}',U2{i}'); 
                    else
                        temp = BuildConjMap(U1{i},U2{i});
                    end

                    [r,c,v] = find(temp);
                    r = r(:)+indx-1;
                    c = c(:)+s-1;
                    v = v(:);
                    rarry = [rarry;r];
                    carry = [carry;c];
                    varry = [varry;v];

                    indx = indx+size(temp,1);
                    
                end
    
            end
           
            T = sparse(rarry,carry,varry,numNonPSDVar+numPSDVarAfterConj,self.NumVar);
            
        end
        
        function crossTerm = CrossTerms(self,x,U,V)
            
            Tuv = self.BuildMultMap(V,U);
            VtXU = Tuv*x';
            s = 1;
            crossTerm = {};
            
            for i=1:length(self.K.s)

                numColV = size(V{i},2);
                numColU = size(U{i},2);

                if (numColU + numColV == self.K.s(i))
                    xtemp = reshape(VtXU(s:s+numColV*numColU-1),numColV,numColU);
                    crossTerm{end+1} = xtemp;
                else
                    error(['numcol(U)+numcol(V) does not equal cone size(',num2str(self.K.s(i)),')']);
                end

            end
        
        end
        
        function [x,x11m,x12m,x22m] = ConjBlock2by2(self,x11,x22,x12,U,V)
            
            s = self.GetIndx('s',1);
            x = sparse(1,s-1);
            s1 = 1; s2 = 1; s3 = 1;
            
            x11m = {}; x12m = {}; x22m = {};
                         
            for i=1:length(self.K.s)
                
                numColV = size(V{i},2);
                numColU = size(U{i},2);
                
                if (numColU + numColV == self.K.s(i))
                    x11i = solUtil.mat(x11(s1:s1+numColU*numColU-1));
                    x22i = solUtil.mat(x22(s2:s2+numColV*numColV-1));
                    x22i = (x22i + x22i')/2;
                    x12i = reshape(x12(s3:s3+numColV*numColU-1),numColU,numColV);
                else
                    error(['numcol(U)+numcol(V) does not equal cone size(',num2str(self.K.s(i)),')']);
                end
     
                s1 = s1+length(x11i(:));
                s2 = s2+length(x22i(:));
                s3 = s3+length(x12i(:));
                  
                x11m{i} = x11i;
                x22m{i} = x22i;
                x12m{i} = x12i;
                
                vx12u = U{i}*x12i*V{i}';
                temp = U{i}*x11i*U{i}'+V{i}*x22i*V{i}'+vx12u+vx12u';
                x = [x,temp(:)'];
                if max(max(abs((temp)'-(temp)))) > 10^-3
                    error('not sym')
                end
                
            end
            
        end


end



end
