function setup()

    display('Updating path...')
    addpath([pwd,'/util'])
    addpath([pwd,'/test'])
    addpath(pwd)

    %Checking dependencies
    display('Checking dependencies...')

    if (isempty(LPSolver.GetSolver()))
       error('No LP solver found. Add linprog, Gurobi, SeDuMi, or Mosek to path.') 
    end

    display('Done!')
    display('Type frlibTests to check installation.')
