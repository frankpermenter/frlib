classdef LPSolver

    methods(Static)

        function solverHandle = GetSolver

            solverExists = []; solverMethod = {};

            solverExists(end+1) = LPSolver.CheckForMosek();
            solverMethod{end+1} = @LPSolver.SolveLPMosek;  

            solverExists(end+1) = ~isempty(which('gurobi'));
            solverMethod{end+1} = @LPSolver.SolveLPGurobi;  

            solverExists(end+1) = ~isempty(which('sedumi'));
            solverMethod{end+1} = @LPSolver.SolveLPSedumi;  

            solverExists(end+1) = ~isempty(which('linprog'));
            solverMethod{end+1} = @LPSolver.SolveLPlinprog;  
            
            
            for i=1:length(solverExists)
                if solverExists(i)
                    solverHandle = solverMethod{i};
                    return 
                end
            end

            solverHandle = [];

        end

        function exists = CheckForMosek
            
           pathlinProg = lower(which('linprog')); 
           exists = ~isempty(strfind(pathlinProg,'mosek'));
            
        end
        
        
        function [x,numErr,infeas] = SolveLP(c,Aineq,bineq,Aeq,beq,lbnd,ubnd)

            solverHandle = LPSolver.GetSolver();
            if isempty(solverHandle) 
                error('No solvers found!')
            else
                [x,infeas,numErr] = solverHandle(c,Aineq,bineq,Aeq,beq,lbnd,ubnd);
            end

        end

        function [x,numErr,infeas] = SolveLPlinprog(c,Aineq,bineq,Aeq,beq,lbnd,ubnd,noOptions)
            
            if ~exist('noOptions','var')
                noOptions = 0;                
            end
            
            if noOptions
                [x,~,flag] = linprog(c,Aineq,bineq,Aeq,beq,lbnd,ubnd);
            else
                options.Display='off';
                [x,~,flag] = linprog(c,Aineq,bineq,Aeq,beq,lbnd,ubnd,options);  
            end
            
            infeas = flag == -2 | flag == -5;
            numErr = ~infeas && flag ~= 1;
            
        end
        
        
        
        function [x,numErr,infeas] = SolveLPMosek(c,Aineq,bineq,Aeq,beq,lbnd,ubnd)
            
            noOptions = 1;
            [x,numErr,infeas] = LPSolver.SolveLPlinprog(c,Aineq,bineq,Aeq,beq,lbnd,ubnd,noOptions);
            
        end

        function [x,infeas,numErr] = SolveLPGurobi(c,Aineq,bineq,Aeq,beq,lbnd,ubnd)

            infeas = 1; numErr = 1;
            model.A = sparse([Aineq;Aeq]);
            model.obj = full(c);
            model.rhs = full([bineq;beq]);
            model.lb = full(lbnd);
            model.ub = full(ubnd);
            model.sense = char( ['<'*ones(length(bineq),1);'='*ones(length(beq),1)]);

            params.method = 1;
            params.crossover = 0;
            params.outputflag = 0;
            result = gurobi(model,params);
            flag = result.status;

            if (strcmp(flag,'INFEASIBLE'))
                infeas = 1;
                numErr = 0;
                x = [];
            else 
                if (strcmp(flag,'OPTIMAL'))
                    x = result.x;
                    infeas = 0;
                    numErr = 0;
                end
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

