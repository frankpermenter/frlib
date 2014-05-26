classdef LPSolver

    methods(Static)

        function solverHandle = GetSolver

            solverExists = []; solverMethod = {};

            solverExists(end+1) = ~isempty(which('linprog'));
            solverMethod{end+1} = @LPSolver.SolveLPMosek;  

            solverExists(end+1) = ~isempty(which('gurobi'));
            solverMethod{end+1} = @LPSolver.SolveLPGurobi;  

            solverExists(end+1) = ~isempty(which('sedumi'));
            solverMethod{end+1} = @LPSolver.SolveLPSedumi;  

            for i=1:length(solverExists)
                if solverExists(i)
                    solverHandle = solverMethod{i};
                    return 
                end
            end

            solverHandle = []

        end

        function [x,numErr,infeas] = SolveLP(c,Aineq,bineq,Aeq,beq,lbnd,ubnd)

            solverHandle = LPSolver.GetSolver();
            if isempty(solverHandle) 
                error('No solvers found!')
            else
                [x,infeas,numErr] = solverHandle(c,Aineq,bineq,Aeq,beq,lbnd,ubnd);
            end

        end

        function [x,numErr,infeas] = SolveLPMosek(c,Aineq,bineq,Aeq,beq,lbnd,ubnd);
            %calls mosek if in path
            [x,~,flag] = linprog(c,Aineq,bineq,Aeq,beq,lbnd,ubnd);
            infeas = flag == -2 | flag == -5;
            numErr = ~infeas && flag ~= 1;
        end

        function [x,infeas,numErr] = SolveLPGurobi(c,Aineq,bineq,Aeq,beq,lbnd,ubnd);

            model.A = sparse([Aineq;Aeq]);
            model.obj = full(c);
            model.rhs = full([bineq;beq]);
            model.lb = full(lbnd);
            model.ub = full(ubnd);
            model.sense = char( ['<'*ones(length(bineq),1);'='*ones(length(beq),1)]);

            result = gurobi(model,params);
            flag = result.status;

            if (strcmp(flag,'INFEASIBLE'))
                infeas = 1;
                x = [];
                return;
            else
                x = result.x;
            end

        end

        function [x,infeas,numErr] = SolveLPSedumi(c,Aineq,bineq,Aeq,beq,lbnd,ubnd)

            Meq = size(Aeq,1);
            N = size(Aeq,2);

            AineqBnd = speye(N);
            AineqUBnd = AineqBnd( ~isinf(ubnd), :);
            bineqUBnd = ubnd(~isinf(ubnd));

            AineqBnd = speye(N);
            AineqLBnd = -AineqBnd( ~isinf(lbnd), :);
            bineqLBnd = -lbnd(~isinf(lbnd));

            AineqEq = [Aineq;AineqUBnd;AineqLBnd];
            bineqEq  = [bineq;bineqUBnd(:);bineqLBnd(:)];

            numSlack = size(AineqEq,1);

            K.l = numSlack;
            K.f = size(Aeq,2);

            A = [Aeq,sparse(Meq,numSlack); ...
                 AineqEq, speye(numSlack) ];

            b = [beq;bineqEq];
            c = [c;sparse(numSlack,1)];

            pars.fid = 0;
            try
                [xopt,~,info] = sedumi(A,b,[],K,pars);
            catch
                [A,b] = cleanLinear(A,b,1);
                [xopt,~,info] = sedumi(A,b,[],K,pars);
            end
            x = sparse( xopt(1:N,1));

            infeas = info.pinf == 1;
            numErr = info.numerr > 1;

        end

    end

end

