This repo contains matlab code for pre-processing SDPs.  Given an SDP that fails Slater''s condition, the code searches
for a lower dimensional face of the PSD cone containing the feasible set. If the search succeeds, the code reformulates the SDP explicitly
over this face.  This results in an SDP with "smaller" semidefinite constraints, i.e. if the original SDP contained a single nxn semidefinite
constraint the reformulation will contain a single dxd constraint with d < n.
  
To find a face efficiently, the code employs approximations of the PSD cone. Better approximations can lead to smaller reformulations at the cost of more computation.  Approximations currently supported include non-negative diagonal matrices and PSD diagonally-dominant matrices.  This approximations are polyhedral,
and hence allow one to search for a face using linear programming.

The code takes as input a primal-dual SDP pair expressed in SeDuMi formatted inputs A,b,c,K:
```Matlab
prg = frlibPrg(A,b,c,K);
```
The primal and/or the dual may fail Slater''s condition. You can choose which problem to reformulate (primal or dual).  

To reformulate the dual using non-negative diagonal approximations ('d') or diagonally-dominant ('dd') approximations, issue a call
of the form:

```Matlab
prgR = prg.ReduceDual('d');
```
To solve the reformulated SDP, call:
```Matlab
[~,y] = prgR.Solve();
```
The solution y solves the original dual SDP.  

To simplify the primal, call
```Matlab
prgR = prg.ReducePrimal('d');
```
To recover a solution x to the original primal SDP call
```Matlab
[xR] = prgR.Solve();
x = prgR.RecoverPrimal(xR);
```

