function [indxZeroOut,Knew] = FindMustVanish(K,indxArray)
  

    %If variables in indxarray vanish, find all other variables that
    %vanish as an implication.

    indxArray = unique(indxArray);
    indxZeroOut = indxArray';
    Knew = K;
    if isempty(indxArray)
        return   
    end

    for i=1:length(indxArray)
        indx = indxArray(i);

        if ( indx <= K.f) 
            Knew.f = Knew.f-1;
            continue;
        end

        if (indx <= K.l+K.f)
            Knew.l = Knew.l - 1;
            continue;
        end

        if (indx <= (sum(K.q) + K.l + K.f) )
            
                offset = K.f + K.l;
                cnstStart = [0,cumsum(K.q)]+1 + offset;
                cnstStart = cnstStart(1:end-1);
                
                cnstNum = max(find(cnstStart <= indx));
                
                if indx ~= cnstStart(cnstNum)
                    Knew.q(cnstNum) = Knew.q(cnstNum) - 1
                else
                    %if a start variable, delete all the others
                    indxZero = cnstStart(cnstNum)+1:(cnstStart(cnstNum)+K.q(cnstNum)-1);
                    Knew.q(cnstNum) = 0;
                    indxZeroOut = [indxZeroOut,indxZero];
                end
                
                continue;
                
        end
        
        if (indx <= (sum(K.q) + K.l + K.f + sum(K.r)))
         
                offset = K.f + K.l + sum(K.q);
                cnstStart = [0,cumsum(K.r)]+1 + offset;
                cnstStart = cnstStart(1:end-1);
                
                cnstNum = max(find(cnstStart <= indx));
                cnstStart = cnstStart(cnstNum);
                

                if indx ~= cnstStart && indx ~= cnstStart+1
                    Knew.r(cnstNum) = Knew.r(cnstNum) - 1;
                else
                    
                    %if a start variable, delete all but the first two
                    indxZero = cnstStart+2:(cnstStart+K.r(cnstNum)-1);
                    indxZeroOut = [indxZeroOut,indxZero];
                    
                    
                    %If both starting vars vanish, remove constraint, else
                    %flag that the constraint should be set to a linear
                    %inequality
                    if any(cnstStart == indxZeroOut) && any( cnstStart+1 == indxZeroOut)
                         Knew.r(cnstNum) = 0;
                    else
                         Knew.r(cnstNum) = 1;
                    end
                    
                end
                
                
                continue;
                
        else
          
               offset = K.f + K.l + sum(K.q) + sum(K.r);
               
               %where do new SDP vars start
               sdpStart = [0,cumsum(K.s.^2)]+1;
               sdpStart = sdpStart(1:end-1) + offset;
               
               %which SDP constraint are we in
               cnstNum = max(find(sdpStart <= indx));
    
               N = K.s(cnstNum);
               offset = sdpStart(cnstNum)-1;
               indxCnstr = indx-offset;
               diagPos = [1:N+1:N*N]';

               findDiag = diagPos(diagPos == indxCnstr);
               
               if isempty(findDiag)
                   continue;
               end

               %find variables on same row and col
               row = floor( (findDiag-1)/N)+1;
               startCol = N*(row-1) + 1;
               indxZeroCol = startCol:startCol+N-1;
               indxZeroRow = row:N:row+N*N;
               indxZeroRow = indxZeroRow(indxZeroRow <= N*N);
               indxZero = union(indxZeroCol,indxZeroRow)+offset;
  
                 

               %decrement dimension of LMI by 1
               Knew.s(cnstNum) =  Knew.s(cnstNum) - 1; 
               indxZeroOut = [indxZeroOut,indxZero];
        end
        
        

    end
end
