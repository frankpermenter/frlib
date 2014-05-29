#Overview
This repo contains matlab code for pre-processing SDPs. Given an SDP, the code searches for a lower dimensional face of the PSD cone containing the feasible set. If the search succeeds, the code reformulates the SDP explicitly over this face. This results in an SDP with "smaller" semidefinite constraints, i.e. if the original SDP contained a single nxn semidefinite constraint the reformulation will contain a single dxd constraint with d < n.

Assuming the SDP is feasible, reformulations are--in principal--possible if and only if the SDP fails Slater's condition.  This code imposes further restrictions on the type of SDPs that can be simplified by employing a user-specified approximation of the PSD cone.  Better approximations widen the class of SDPs that can be reformulated,
but increase the cost of finding a face.

###Input formats
The code takes in a primal-dual SDP pair expressed using  SeDuMi formatted inputs A,b,c,K.  The code supports pre-processing of the either the primal or the dual SDP.


## Example (reduction of primal):
To perform reductions, one specifies a  PSD approximation and the SDP (primal or dual) ones wishes to reduce. A typical use case is given below.
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
The above examples  used diagonal approximations ('d').  It is also possible to use diagonally-dominant matrices ('dd') for the PSD approximation, i.e. one can replace the relevant lines of the above examples with the following:

```Matlab
prgR = prg.ReducePrimal('dd');
```
or
```Matlab
prgR = prg.ReduceDual('dd');
```
While this approximation is more expensive than the diagonal approximation ('d'), it widens the class of SDPs that can be simplified and may lead to smaller reformulations.
