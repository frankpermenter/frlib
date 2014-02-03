function [varToSolveForValues,varToSolveFor,tFormRemoveVars,tFormAddVars,Aup,bup] = PreSolveLinearEq(Aeq,beq);

    [~,n] = size(Aeq);
    eqsInOneVar = find(sum(Aeq~=0,2) == 1);
    [eqsInOneVarPermute,varToSolveFor] = find(Aeq(eqsInOneVar,:));
    eqsInOneVar = eqsInOneVar(eqsInOneVarPermute);

    AeqEntries = sub2ind([size(Aeq,1),size(Aeq,2)],eqsInOneVar,varToSolveFor);
    varToSolveForValues = beq(eqsInOneVar)./Aeq(AeqEntries);

    varsRemove = varToSolveFor;
    varsKeep = setdiff(1:n,varsRemove);
 
    x0 = sparse(n,1);
    x0(varToSolveFor) = varToSolveForValues;
 
    bup = beq - Aeq*x0;
    Aup = Aeq(:,varsKeep);


    %remove variables: xRed = tFormRemoveVars*x
    tFormRemoveVars = sparse(1:length(varsKeep),varsKeep,1,length(varsKeep),n);

    %add back variables: x = tFormAddVars xRed + x0
    tFormAddVars = sparse(varsKeep,1:length(varsKeep),1,n,length(varsKeep));





