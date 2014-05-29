#Overview
This repo contains matlab code for pre-processing SDPs. Given an SDP that fails Slater's condition, the code searches for a lower dimensional face of the PSD cone containing the feasible set. If the search succeeds, the code reformulates the SDP explicitly over this face. This results in an SDP with "smaller" semidefinite constraints, i.e. if the original SDP contained a single nxn semidefinite constraint the reformulation will contain a single dxd constraint with d < n.

###Input formats
The code takes in a primal-dual SDP pair expressed using  SeDuMi formatted inputs A,b,c,K.  The code supports pre-processing of the either the primal or the dual SDP.


## Example (reduction of primal):
To perform reductions, you must specify the type of PSD approximation and the SDP (primal or dual) you wish the reduce. A typical use case is given below.
```Matlab
prg = frlibPrg(A,b,c,K);
prgR = prg.ReducePrimal('d');
xR = prgR.Solve();
x = prgR.RecoverPrimal(xR);

```
The second line 
```Matlab
prgR = prg.ReducePrimal('d');
```
tells the code to reduce the primal SDP using a diagonal approximation ('d').  The next line solves the reduced SDP (by calling SeDuMi).  The last line reconstructs a solution to the original primal SDP described by A,b,c,K.


## Example (reduction of dual):
Reduction of the dual is similarly done, but doesn't require any solution recovery (i.e. a solution to the reduced dual solves the original dual described by A,b,c,K):

```Matlab
prg = frlibPrg(A,b,c,K);
prgR = prg.ReduceDual('d');
[~,y] = prgR.Solve(); %y also solves original dual

```

## Other approximations
The above examples all used diagional approximations ('d').  It is also possible to use diagonally-dominant matrices ('dd'), e.g.

```Matlab
prgR = prg.ReduceDual('dd');
prgR = prg.ReducePrimal('dd');
```

