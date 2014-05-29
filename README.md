#Overview
This repo contains matlab code for pre-processing SDPs.  The code--frlib--reformulates an SDP that fails Slater's condition 
into a smaller, equivalent SDP using facial reduction and approximations of the PSD cone. The reformulation is "smaller"
in the sense that semidefinite constraints involve lower dimensional PSD cones. For instance, if the original problem had a single
n \times n semidefinite constraint, the reformulation will have a d x d semidefinite constraint with d < n.
  
Better approximations lead to more simplifications at the cost of more computation.  Approximations
currently supported include non-negative diagonal matrices and PSD diagonally-dominant matrices.

##Input formats
The code takes in a primal-dual SDP pair formatted in SeDuMi format.



##Example Usage:
To reduce (i.e. reformulate) an SDP expressed in terms of SeDuMi formatted inputs A,b,c,K, first call:

```Matlab
prg = frlibPrg(A,b,c,K);
```

The problem data A,b,c,K describe a primal-dual SDP pair.   You must choose which problem to reduce (primal or dual).

To reduce the dual using non-negative diagonal approximations ('d') or diagonally-dominant ('dd') approximations, issue a call
of the form:

```Matlab
prgR = prg.ReduceDual('d');
```
To solve the simplified SDP, call:
```Matlab
[~,y] = prgR.Solve();
```
The solution y solves the original dual SDP.  

To reduce the primal, call
```Matlab
prgR = prg.ReducePrimal('d');
```
To recover a solution x to the original primal SDP call
```Matlab
[xR] = prgR.Solve();
x = prgR.RecoverPrimal(xR);
```

