#Overview
This repo contains MATLAB code for pre-processing SDPs not strictly feasible. Given an SDP, the code searches for a lower dimensional face of the PSD cone containing the feasible set. If the search succeeds, the code reformulates the SDP explicitly over this face. This results in an SDP with "smaller" semidefinite constraints, i.e. if the original SDP contained a single nxn semidefinite constraint the reformulation will contain a single dxd constraint with d < n.

To make the search for a face 'easy' (and practical for a pre-processor), the code employs a user-specified approximation of the PSD cone.   Better approximations widen the class of SDPs that can be reformulated, but increase the cost of finding a face. 

The code takes in a primal-dual SDP pair expressed using  SeDuMi formatted inputs A,b,c,K and supports pre-processing of  either the primal or the dual SDP. Approximations currently supported are non-negative diagonal matrices and diagonally-dominant matrices.


## Example (reduction of primal):
To perform reductions, one specifies a  PSD approximation and the SDP (primal or dual) one wishes to reduce. A typical use case is given below.
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
To use diagonally-dominant matrices ('dd') for the PSD approximation, one can update the relevant lines of the above  examples as follows:

```Matlab
prgR = prg.ReducePrimal('dd');
```
```Matlab
prgR = prg.ReduceDual('dd');
```
While this approximation is more expensive than the diagonal approximation ('d'), it widens the class of SDPs that can be simplified and may lead to smaller reformulations.
