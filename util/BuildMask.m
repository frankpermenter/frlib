function [cliques,Ar,cr,Kr,indx] = BuildMask(A,b,c,K)

c = c(:)';
if nnz(c(1:K.f)) > 0
    error('Cannot eliminate free variables if they appear in objective');
end

if K.q + K.r  > 0
    error('Lorentz cone constraints not yet supported.');
end

Ain = A;
bin = b;
Kin = K;
cin = c;

Kin = coneBase.cleanK(Kin);
[A,b,K] = EliminateFreeVars(A,b,K);
c = c(Kin.f+1:end);
[A,~] = ConsolidateLinearAndPSDConstraints(A,K);
[c,K] = ConsolidateLinearAndPSDConstraints(c,K);

t = coneBase(K);
Apsd = [A(:,t.Kstart.s(1):end);c(:,t.Kstart.s(1):end)];
b(end+1,1) = 1;

Apsd = double(Apsd~=0);
tau = double(b~=0);
tauPrev = tau;

while(1)

    M = solUtil.mat( any(Apsd(tau~=0,:) ));
    [M,cliques] = BinaryPsdCompletion(M);
    Mvect = M(:)';
    [~,indx] = find(Mvect);
    [keep,~] = find(Apsd(:,indx));

    tau(keep) = 1;
    if all(tau == tauPrev)
        break;
    else
        tauPrev = tau;
    end

end

nnz(Mvect)/(Kin.l+sum(Kin.s.^2))
indx = zeros(nnz(Mvect),1); e = 0;
for i=1:length(cliques)
    [p,q] = meshgrid(cliques{i}, cliques{i});   
    s = e+1;
    e = s+length(cliques{i})^2-1;
    indx(s:e) = sub2ind(size(M),p(:),q(:));
end

Ktemp  = Kin; Ktemp.f = 0;
A = ConsolidateLinearAndPSDConstraints(Ain(:,Kin.f+1:end),Ktemp);  
c = ConsolidateLinearAndPSDConstraints(cin(:,Kin.f+1:end),Ktemp);  
Ar = [Ain(:,1:Kin.f),A(:,indx)];
cr = [cin(1:Kin.f),c(indx)];
Kr.f = Kin.f;
Kr.s = cellfun(@(x)size(x,1),cliques,'UniformOutput',1);
Kr.l = sum(Kr.s==1);
Kr.s = Kr.s(Kr.s~=1);









