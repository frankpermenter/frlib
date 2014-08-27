#Overview
This repo contains MATLAB code for pre-processing SDPs using techniques described in the paper *Partial facial reduction: simplified, equivalent SDPs via approximations of the PSD cone* by Permenter and Parrilo.  The code is still under active development, so interface changes are possible.  Bug reports are also welcomed and appreciated!

Given an SDP, the code searches for a lower dimensional face of the PSD cone containing the feasible set. If the search succeeds, the code reformulates the SDP explicitly over this face. This results in an SDP with "smaller" semidefinite constraints, i.e. if the original SDP contained a single nxn semidefinite constraint the reformulation will contain a single dxd constraint with d < n.

The search space for a face is controlled by a specified approximation of the PSD cone.   The better approximation, the
larger the search space.  By specifying the approximation, the user trade offs pre-processing effort with potential problem simplifications.

The code takes in a primal-dual SDP pair expressed using  SeDuMi formatted inputs A,b,c,K and supports pre-processing of  either the primal or the dual SDP. Approximations currently supported are non-negative diagonal matrices and diagonally-dominant matrices.


## Example (reduction of primal):
To perform reductions, one specifies a  PSD approximation and the SDP (primal or dual) one wishes to reduce. A typical use case is given below:
```Matlab
prg = frlibPrg(A,b,c,K);
prgR = prg.ReducePrimal(`d');
[x_reduced,y_reduced] = prgR.Solve();
[x,y,dual_recov_success] = prgR.Recover(x_reduced,y_reduced);

```
The second line 
```Matlab
prgR = prg.ReducePrimal('d');
```
tells the code to reduce the primal SDP using a diagonal approximation ('d').  The next line solves the reduced SDP (by calling SeDuMi).  The last line reconstructs a solution to the original primal SDP described by A,b,c,K.  It also attempts to recover a solution y to the original dual.  Success
of this attempt is given by the flag dual_recov_success.


## Example (reduction of dual):
Reduction of the dual is similarly done:

```Matlab
prg = frlibPrg(A,b,c,K);
prgR = prg.ReduceDual(`d');
[x_reduced,y_reduced] = prgR.Solve();
[x,y,primal_recov_success] = prgR.Recover(x_reduced,y_reduced);
```
Note the call to 
```Matlab
[x,y,primal_recov_success] = prgR.Recover(x_reduced,y_reduced);
```
now reconstructs a solution to the original dual SDP and attempts recovery of a solution to the primal.


## Other approximations
To use diagonally-dominant matrices ('dd') for the PSD approximation, one can update the relevant lines of the above  examples as follows:

```Matlab
prgR = prg.ReducePrimal('dd');
```
```Matlab
prgR = prg.ReduceDual('dd');
```
While this approximation is more expensive than the diagonal approximation ('d'), it widens the class of SDPs that can be simplified and may lead to smaller reformulations.
