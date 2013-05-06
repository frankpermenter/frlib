function [act_flg,infeas,gplkError] = lpact(A_ineq,b_ineq,A_eq,b_eq)
%Identifies the always active inequalities in the set S:
% S = { x : 
%   A_eq x = b_ineq
%   A_ineq x <= b_eq }
%  
% If S is not empty:
%   act_flag(i) = 1 if ( A_ineq x )_i = (b_ineq)_i for all x in S 
%   act_flag(i) = 0  otherwise
%   infeas = 0
% 
% If S is empty (i.e. the system is infeasible):
%   act_flag = -1
%   infeas = 1

eps = 10^-12;
gplkError = 0;

NumEqIn = length(b_ineq);
if ~exist('A_eq','var')
	A_eq=[];
end


NumIneq = length(b_ineq);
NumEq = length(b_eq);
LenX = size(A_ineq,2);

A = [A_ineq,eye(NumIneq),-b_ineq;];
b = zeros(NumIneq,1);

if (~isempty(b_eq))
    A = [A;A_eq,zeros(NumEq,NumIneq),-b_eq];
    b = [b;zeros(NumEq,1)];
end


c = [zeros( LenX,1);-ones(NumIneq,1);0];
lbnd = [-Inf*ones(LenX,1);zeros(NumIneq,1);1];
ubnd = [Inf*ones(LenX,1);ones(NumIneq,1);Inf];
    
GLPK = 0;
if GLPK == 1
    
    ctype = char( ['U'*ones(NumIneq,1);'S'*ones(NumEq,1)]); 


    [xopt,copt,flag,extra] = glpk(c,A,b,lbnd,ubnd,ctype);


    GLPK_INFEAS = 110;
    if (flag == GLPK_INFEAS)
        act_flg = -1;
        infeas = 1;
        return
    end
    
else

    model.A = sparse(A);
    model.obj = c;
    model.rhs = full(b);
    model.lb = lbnd;
    model.ub = ubnd;
    model.sense = char( ['<'*ones(NumIneq,1);'='*ones(NumEq,1)]);
    params.outputflag = 0; 
    result = gurobi(model,params);
    flag = result.status;
    
    if (strcmp(flag,'INFEASIBLE'))
        act_flg = -1;
        infeas = 1;
        return
    else
        xopt = result.x;
    end
end

    %Sanity check on GPLK answer
    if any(xopt - lbnd < -eps)
        act_flg = zeros(NumEqIn,1);
        gplkError = 1;
        return
    end
    
    act_flg = abs(xopt(LenX+1: LenX + 1 + NumEqIn-1)) < eps;
    infeas = 0;
    

