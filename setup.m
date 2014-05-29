function setup()


    display('Updating path...')
    addpath([pwd,'/util'])
    addpath([pwd,'/test'])
    addpath(pwd)


    %Checking dependencies
    sedumiExists = ~isempty(which('sedumi'));
    display('Checking dependencies...')
    if ~sedumiExists
        error('Cannot find SeDuMi. Aborting...')
    end

    if (isempty(LPSolver.GetSolver()))
       error('No LP solver found. Add linprog, Gurobi, SeDuMi, or Mosek to path.') 
    end

    display('Done!')
    display('Type frlibTests to check installation.')
