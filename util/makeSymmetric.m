function [A,C] = makeSymmetric(A,C,K)
    Z = ConeBase(K);
    A = Z.symmetrize(A);
    C = Z.symmetrize(C);
end
