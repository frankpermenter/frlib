function [cliques,Ar,cr,Kr,indx,M] = BuildMask(A,b,c,K)

if K.q + K.r  > 0
    error('Lorentz cone constraints not yet supported.');
end

M = full(c~=0);
M = M(:);
nnz_M = nnz(M);
cone = coneBase(K);
Kr = K;

cliques = {};

while(1)

    [M] = SubspaceClosureCoordDisjointSupport(M,A,b);
    for i=1:length(K.s)
        [s,e] = cone.GetIndx('s',i);
        [temp,cliques_i] = BinaryPsdCompletion(solUtil.mat(M(s:e)));
        M(s:e) = temp(:);
        Kr.s(i) = length(cliques_i);
        cliques{i} = cliques_i;
    end
    
    if (nnz_M == nnz(M))
       break; 
    end
        
    nnz_M = nnz(M);
    
end


[s,e] = cone.GetIndx('f',1);
indx = find(M(s:e));
Kr.f = length(indx);

[s,e] = cone.GetIndx('l',1);
indx = find(M(s:e));
Kr.l = length(indx);

[s,e] = cone.GetIndx('s',1);
indx =  find(M(1:s-1));

Kr.s = [];
for i=1:length(K.s)
    for j=1:length(cliques{i})
        indx = [indx;cone.SubMatToIndx(cliques{i}{j},i  )];
        Kr.s(end+1) = length(cliques{i}{j});
    end
end

Ar = A(:,indx);
br = b;
cr = c(indx);

end


function [S] = SubspaceClosureCoordDisjointSupport(S,A,b)

    isSparse = issparse(S);
    M = S';
    M = any(A(b~=0,:),1) | M; %the stuff we must pass through
    tau = any(A(:,M>0)~=0,2); %all rows at least partially passed through
    S = any(A(tau,:),1)';
    
    %S will be sparse if A is sparse.  Undo to match input.
    if ~isSparse
        S = full(S); 
    end
    
end


