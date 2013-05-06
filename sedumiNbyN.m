function A= sedumiNbyN(K,V)
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


