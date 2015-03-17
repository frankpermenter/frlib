function setup()

    display('Updating path...')
    directory = fileparts(which('setup.m'));
    addpath([directory,'/util'])
    addpath([directory,'/test'])
    addpath(directory)

    %Checking dependencies
    display('Checking dependencies...')

    if (isempty(LPSolver.GetSolver()))
       error('No LP solver found. Add linprog, Gurobi, SeDuMi, or Mosek to path.') 
    end

    display('Done!')
    display('Type frlibTests to check installation.')
