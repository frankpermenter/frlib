function runTests


%%%%%%%%%%%%%%%%%%%
%Test no reductions
%%%%%%%%%%%%%%%%%%%

sdpvar y1 y2 y3
Constraint = [y1 y2 0 ;
               y2  y1 0;
               0  0  y3-1]>=0;
Constraint = [Constraint;y1==3;y2==2];
Objective =  y3;
           
ops = sdpsettings('solver','frlib','frlib.solver','sdpt3');
ops = sdpsettings(ops,'frlib.removeDualEq',1);
ops = sdpsettings(ops,'frlib.useQR',1);
ops = sdpsettings(ops,'frlib.approximation','d');





%Dual
ops = sdpsettings(ops,'dualize',0);
ops = sdpsettings(ops,'frlib.reduce','primal');

stats = solvesdp(Constraint,Objective,ops);

if  double(norm([y1,y2,y3]-[3,2,1])) > 10^-3
   error('Incorrect results');
end

if (stats.problem ~= 0)
    error('Solver Error');
end


%Primal
ops = sdpsettings(ops,'dualize',1);
ops = sdpsettings(ops,'frlib.reduce','primal');

stats = solvesdp(Constraint,Objective,ops);

if  double(norm([y1,y2,y3]-[3,2,1])) > 10^-3
   error('Incorrect results');
end

if (stats.problem ~= 0)
    error('Solver Error');
end


%%%%%%%%%%%%%%%%%%%%%%
%Test basic reductions
%%%%%%%%%%%%%%%%%%%%%%

sdpvar y1 y2 y3
Constraint = [y1-pi 0 0 ;
               0 -y1+pi y2;
               0 y2 y2+y3]>=0;
Objective =  y3;
           
         
ops = sdpsettings(ops,'frlib.reduce','dual');
 

stats = solvesdp(Constraint,Objective,ops);

if  double(norm([y1,y2,y3]-[pi,0,0])) > 10^-3
   error('Incorrect results');
end

if (stats.problem ~= 0)
    error('Solver Error');
end




ops = sdpsettings(ops,'frlib.removeDualEq',1);
stats = solvesdp(Constraint,Objective,ops);

if  double(norm([y1,y2,y3]-[pi,0,0])) > 10^-3
   error('Incorrect results');
end

if (stats.problem ~= 0)
    error('Solver Error');
end




%%%%%%%%%%%%%%%%%%%%%%
%Test SOS Example
%%%%%%%%%%%%%%%%%%%%%%
J = [1 -1 1 1 -1;-1 1 -1 1 1;1 -1 1 -1 1;1 1 -1 1 -1;-1 1 1 -1 1];
m = 5;
x = sdpvar(m,1);
xx = [x;zeros(5-length(x),1)];
p = (sum(xx.^2))*(xx.^2)'*J*(xx.^2);

% Vanilla SOS
ops = sdpsettings('sos.newton',0,'sos.congruence',0,'sos.model',2);
[sol,v,Q] = solvesos(sos(p),[],ops);


% Only partial face reductions
ops = sdpsettings('sos.newton',0,'sos.congruence',0,'sos.model',2);
ops = sdpsettings(ops,'solver','frlib','frlib.reduce','dual','frlib.approximation','dd');
[sol,v,Q] = solvesos(sos(p),[],ops);



ops = sdpsettings('solver','frlib','frlib.solver','sdpt3');
ops = sdpsettings(ops,'frlib.removeDualEq',1);
ops = sdpsettings(ops,'frlib.useQR',1);
ops = sdpsettings(ops,'frlib.reduce','dual','frlib.approximation','dd');
yalmiptest(ops,1);

display('Press F5 to continue tests...')
keyboard;
ops = sdpsettings(ops,'frlib.reduce','primal','frlib.approximation','dd');
yalmiptest(ops,1);


