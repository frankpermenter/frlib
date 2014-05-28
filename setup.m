function setup()

    %Checking dependencies
    sedumiExists = ~isempty(which('sedumi'));
    display('Checking dependencies...')
    if ~sedumiExists
        error('Cannot find SeDuMi. Aborting...')
    end

    if (isempty(LPSolver.GetSolver()))
       error('No LP solver found. Add linprog, Gurobi, SeDuMi, or Mosek to path.') 
    end

    display('Updating path...')
    addpath([pwd,'/util'])
    addpath([pwd,'/test'])
    addpath(pwd)

    display('Done!')
    display('Type runTests to check installation.')
