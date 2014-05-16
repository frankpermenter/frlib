%computes T for which Tq(:) = [U1 q U2'] (:)
function T =  BuildConjMap(U1,U2)

    tranposeMap = 1;
    
    %For efficiency, swap U1 and U2 depending on dimensions 
    if (size(U1,2) < size(U2,2))
        temp = U2;
        U2 = U1;
        U1 = temp;
        tranposeMap = 0;
    end
    
    [r1,~] = size(U1);
    [r2,c2] = size(U2);
    
    %compute T for which Tq(:) = [U1 q U2']' (:)      
    T1 = kron(speye(c2),U1);    %linear map of q(:) to (U1*q)(:)
    T2 = transpose(r1,c2);      %linear map of (U1q)(:) to (U1q)'(:)
    T3 = kron(speye(r1),U2);    %linear map of (U1q)'(:) to U2*(U1q)'(:)

    T = T3*T2*T1;

    %compute tranpose unless we flipped U1 and U2 
    if (tranposeMap == 1)
        T = transpose(r2,r1)*T;
    end
    
end


function test(T,U1,U2)

    [r1,c1] = size(U1);
    [r2,c2] = size(U2);

    Bmat = randn(c1,c2);
    B = Bmat(:);

    c = T*B(:);
    temp= [U1*Bmat*U2'] ;

     norm(temp(:)-c)
end

function T = transpose(r1,c2)

    indx = reshape(1:r1*c2,r1,c2);
    indx = indx';
    indx = indx(:);  
    T = sparse([1:length(indx)]',indx,ones(length(indx),1),length(indx),length(indx));
    
end