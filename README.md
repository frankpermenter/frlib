#Overview
This repo contains matlab code for pre-processing SDPs. Given an SDP that fails Slater's condition, the code searches for a lower dimensional face of the PSD cone containing the feasible set. If the search succeeds, the code reformulates the SDP explicitly over this face. This results in an SDP with "smaller" semidefinite constraints, i.e. if the original SDP contained a single nxn semidefinite constraint the reformulation will contain a single dxd constraint with d < n.

##Input formats
The code takes in a primal-dual SDP pair expressed using  SeDuMi formatted inputs A,b,c,K.  The code supports pre-processing of the either the primal or the dual SDP.


##Example Usage:
To reduce (i.e. reformulate) an SDP, first create an frlibPrg object from the problem data:

```Matlab
prg = frlibPrg(A,b,c,K);
```

To perform reductions, you must specify the type of PSD approximation and the SDP (primal or dual) you wish the reduce.


## Reduction of the dual SDP
To reduce the dual using non-negative diagonal approximations ('d'), issue the call:

```Matlab
prgR = prg.ReduceDual('d');
```

Similarly, to reduce the dual using diagonally-dominant ('dd') approximations, issue the call:
```Matlab
prgR = prg.ReduceDual('dd');
```

To solve the simplified SDP, call:
```Matlab
[~,y] = prgR.Solve();
```
The solution y solves the original dual SDP.  

### Reduction of the primal SDP

To reduce the primal using non-negative diagonal approximations ('d'), issue the call:
```Matlab
prgR = prg.ReducePrimal('d');
```


Similarly, to reduce the primal using diagonally-dominant ('dd') approximations, issue the call:
```Matlab
prgR = prg.ReducePrimal('dd');
```


####Solution recovery
Unlike the dual SDP, extra steps must be taken to recover a solution to the original SDP.  To do this, issue the call:

```Matlab
[xR] = prgR.Solve();
x = prgR.RecoverPrimal(xR);
```

