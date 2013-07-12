function [act_flg,infeas,gplkError] = lpact(A_ineq,b_ineq,A_eq,b_eq,useSedumi)
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

% NumIneq = length(b_ineq);
% NumEq = length(b_eq);
% LenX = size(A_ineq,2);
% A = [A_ineq;A_eq];
% b = [b_ineq;b_eq];
% 
% lbnd = [-Inf*ones(LenX,1);];
% ubnd = [Inf*ones(LenX,1);];
% model.A = sparse(A);
% model.obj = zeros(size(A,2),1);
% model.rhs = full(b);
% model.lb = lbnd;
% model.ub = ubnd;
% model.sense = char( ['<'*ones(NumIneq,1);'='*ones(NumEq,1)]);
% params.outputflag = 0; 
% result = gurobi(model,params);
% flag = result.status;
%     
% 
% result

eps = 10^-12;
gplkError = 0;
infeas = 0;
act_flg = [];

if ~exist('A_eq','var')
	A_eq=[];
end

if ~exist('b_eq','var')
	b_eq=[];
end

if ~exist('useSedumi','var')
	useSedumi = 1;
end


if b_ineq == 0
    b_ineq = sparse(size(A_ineq,1),1);
end

NumEqIn = length(b_ineq);
NumIneq = length(b_ineq);
NumEq = length(b_eq);
LenX = size(A_ineq,2);

A = [A_ineq,speye(NumIneq),-b_ineq;];
b = sparse(NumIneq,1);

if (NumIneq == 0)
    return
end

if (~isempty(b_eq))
    A = [A;A_eq,sparse(NumEq,NumIneq),-sparse(b_eq)];
    b = [b;sparse(NumEq,1)];
end

c = [zeros( LenX,1);-ones(NumIneq,1);0];
lbnd = [-Inf*ones(LenX,1);zeros(NumIneq,1);1];
ubnd = [Inf*ones(LenX,1);ones(NumIneq,1);Inf];

if (useSedumi == 0)

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
	
else
		Aineq = A(1:NumIneq,:);
		Aeq = A(1+NumIneq:end,:);
		bineq = b(1:NumIneq,1);
		beq = b(1+NumIneq:end,1);	
	if (0)
		
		xopt = linprog(c,Aineq,bineq,Aeq,beq,lbnd,ubnd);
		
	else
	
	%c = [zeros( LenX,1);-ones(NumIneq,1);0];
	%lbnd = [-Inf*ones(LenX,1);zeros(NumIneq,1);1];
	%ubnd = [Inf*ones(LenX,1);ones(NumIneq,1);Inf];
	numSlack = size(Aineq,1) + NumIneq*2 + 1;
	
	K.l = numSlack;
	K.f = size(Aeq,2);
	
	%just the equations 
	A = [Aeq, zeros(size(Aeq,1), numSlack)];
	%just slacks for inequalities
	A = [A;[Aineq, eye(NumIneq),sparse(NumIneq,numSlack-NumIneq)]];

	%slacks for bounds on y lower bound
	A = [A;[sparse(NumIneq+1,LenX),speye(NumIneq+1),sparse(NumIneq+1,NumIneq),-speye(NumIneq+1),sparse(NumIneq+1,NumIneq)]];
	
	%upper bound
	A = [A; [sparse(NumIneq,LenX), speye(NumIneq),sparse(NumIneq,numSlack+1-NumIneq),speye(NumIneq)]];
	b = [beq;bineq;sparse(NumIneq,1);ones(NumIneq+1,1)];

	[A,b] = cleanLinear(A,b);
	c = [c;zeros(numSlack,1)];
	[xopt,~,info]=sedumi(A,b,[],K);
	
	
	end
end
    max(xopt(LenX+1: LenX + 1 + NumEqIn-1))
    act_flg = abs(xopt(LenX+1: LenX + 1 + NumEqIn-1)) < eps;
	
    infeas = 0;
    plot(xopt(LenX+1: LenX + 1 + NumEqIn-1))

