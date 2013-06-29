function setup()
    %Checking dependencies
    sedumiExists = ~isempty(which('sedumi'));
    gurobiExists = ~isempty(which('gurobi'));

    display('Checking dependencies...')
    if ~sedumiExists
        error('Cannot find sedumi. Aborting...')
    end

    if ~gurobiExists
        error('Cannot find gurobi. Aborting...')
    end

    display('Updating path...')
    addpath([pwd,'/util'])
    addpath([pwd,'/test'])
    addpath(pwd)

    display('Done!')
    display('Type runTests to check installation')
