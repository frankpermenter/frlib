classdef LPSolver
    
    methods(Static)
        
        function solverHandle = GetSolver
            
            solverExists = []; solverMethod = {};
            
            solverExists(end+1) = ~isempty(which('mosekopt'));
            solverMethod{end+1} = @LPSolver.SolveLPMosek;
            
            solverExists(end+1) = ~isempty(which('gurobi'));
            solverMethod{end+1} = @LPSolver.SolveLPGurobi;
  
            solverExists(end+1) = ~isempty(which('linprog'));
            solverMethod{end+1} = @LPSolver.SolveLPlinprog;          

            solverExists(end+1) = ~isempty(which('sedumi'));
            solverMethod{end+1} = @LPSolver.SolveLPSedumi;
                         
            for i=1:length(solverExists)
                if solverExists(i)
                    solverHandle = solverMethod{i};
                    return
                end
            end
            
            solverHandle = [];
            
        end
        
        function [x,numErr,infeas] = SolveLP(c,Aineq,bineq,Aeq,beq,lbnd,ubnd)
            
            solverHandle = LPSolver.GetSolver();
            if isempty(solverHandle)
                error('No solvers found!')
            else
                [x,numErr,infeas] = solverHandle(c,Aineq,bineq,Aeq,beq,lbnd,ubnd);
            end
            
        end
        
        function [x,numErr,infeas] = SolveLPlinprog(c,Aineq,bineq,Aeq,beq,lbnd,ubnd)
            
            options = optimset(@linprog);
            options.Display='off';
            [x,~,flag] = linprog(c,Aineq,bineq,Aeq,beq,lbnd,ubnd,[],options);
                
            infeas = flag == -2 | flag == -5;
            numErr = ~infeas && flag ~= 1;
                       
        end
        
        
        function [x,numErr,infeas] = SolveLPMosek(c,Aineq,bineq,Aeq,beq,lbnd,ubnd)
            
            param.MSK_IPAR_INTPNT_BASIS = 0;
            param.MSK_IPAR_MAX_NUM_WARNINGS = 0;
            param.MSK_DPAR_OPTIMIZER_MAX_TIME = -1;
            param.MSK_IPAR_INTPNT_MAX_ITERATIONS = 400;
            param.MSK_IPAR_MIO_MAX_NUM_BRANCHES = -1;
            param.MSK_IPAR_LOG = 0;
            param.MSK_IPAR_LOG_INTPNT = 0;
            param.MSK_IPAR_LOG_SIM = 0;
            param.MSK_IPAR_LOG_BI = 0;
            param.MSK_IPAR_LOG_PRESOLVE = 0;
            param.MSK_IPAR_INTPNT_BASIS = 0;

            f = c;
            A = Aineq;
            b = bineq;
            B = Aeq;
            c = beq;
            l = lbnd;
            u = ubnd;
            n = length(f);
            
            % Setup the problem that is feed into MOSEK.
            prob        = [];
            prob.c      = reshape(f,n,1);
            prob.a      = [A;B];
            if ( isempty(prob.a) )
                prob.a = sparse(0,length(f));
            elseif ~issparse(prob.a)
                prob.a = sparse(prob.a);
            end
            prob.blc    = [-inf*ones(size(b));c];
            prob.buc    = [b;c];
            prob.blx    = l;
            prob.bux    = u;
            
            param.MSK_IPAR_INTPNT_BASIS = 0;
            cmd = 'minimize info echo(0) statuskeys(1) symbcon';
            [rcode,res] = mosekopt(cmd,prob,param);
            
            if ( isfield(res,'sol') )
                if ( isfield(res.sol,'itr') )
                    x = res.sol.itr.xx;
                else
                    x = res.sol.bas.xx;
                end
            else
                x = [];
            end            
              
            if rcode==0 && isfield(res,'sol')
                
                if ( isfield(res.sol,'itr') )
                    
                    if (res.sol.itr.solsta==res.symbcon.MSK_SOL_STA_OPTIMAL || ...
                        res.sol.itr.solsta==res.symbcon.MSK_SOL_STA_NEAR_OPTIMAL)
                        flag = 1;
                    elseif (res.sol.itr.solsta==res.symbcon.MSK_SOL_STA_PRIM_INFEAS_CER || ...
                        res.sol.itr.solsta==res.symbcon.MSK_SOL_STA_NEAR_PRIM_INFEAS_CER)
                        flag = -2;
                    elseif (res.sol.itr.solsta==res.symbcon.MSK_SOL_STA_DUAL_INFEAS_CER || ...
                        res.sol.itr.solsta==res.symbcon.MSK_SOL_STA_NEAR_DUAL_INFEAS_CER)
                        flag = -3;
                    end
                    
                end
                
            else
                
                flag = -999;
                
            end
            
            infeas = flag == -2;
            numErr = ~infeas && flag ~= 1;
                                
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
                [A,b] = CleanLinear(A,b,1);
                [xopt,~,info] = sedumi(A,b,[],K,pars);
            end
            x = sparse( xopt(1:N,1));
            
            infeas = info.pinf == 1;
            numErr = info.numerr > 1;
            
        end
                
    end
    
end

