function [A,C] = makeSymmetric(A,C,K)
    Z = coneHelp(A,[],C,K);
    A = Z.symmetrize(A);
    C = Z.symmetrize(C);
end
