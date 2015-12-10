classdef reducedPrg < frlibPrg

    properties(SetAccess=protected)

        faces
        opts
        unreducedPrg
        
    end
          
    properties(Access=protected)
        
        Ty 
        y0
        primalOrDual

    end

    methods(Abstract)

         Recover(self,x,y,eps)
         RecoverPrimal(self,x)
         RecoverDual(self,yinput,eps)

    end

    methods

        function self = reducedPrg(A,b,c,K)
            self = self@frlibPrg(A,b,c,K);
        end

        function [pass,errors] = VerifyRedCert(self,eps)

            if nargin < 2
               eps = 0; 
            end
            
            pass = 1;  errors = {};
            for i=1:length(self.faces)-1
                
                if ~self.faces{i}.InDualCone(self.faces{i+1}.redCert.S,eps)
                    pass = 0;
                    min_eig = eigK(self.faces{i+1}.redCert.S, self.unreducedPrg.K);
                    min_eig = min(min_eig);
                else
                    min_eig = 0;
                end
               
                if isa(self, 'reducedPrimalPrg')
                    err1 = self.faces{i+1}.redCert.y'*self.unreducedPrg.b;
                    err2 = norm(  self.unreducedPrg.A'*self.faces{i+1}.redCert.y- ...
                            self.faces{i+1}.redCert.S');
                end
            
                if isa(self, 'reducedDualPrg')
                    err1 = abs( self.unreducedPrg.c(:)'*self.faces{i+1}.redCert.S);
                    err2 = max( abs( self.unreducedPrg.A*self.faces{i+1}.redCert.S));
                end
                

                if err1 > eps || err2 > eps 
                    pass = 0; 
                end
                
                errors{i} = [err1,err2,min_eig];
            end

        end
            
        function [datar,dataNoRed] = PrintError(self)

            if strcmp(self.primalOrDual,'primal')
                CheckFeas = @(info) info.pinf ~= 1 && info.dinf ~= 1;
                Recover = @(prg,x,y) {prg.RecoverPrimal(x),y};
                DistToFace = @(prg,x,y) norm(x- prg.faces{end}.ProjFace(x))/norm(x);
            else
                CheckFeas = @(info) info.pinf ~= 1 && info.dinf ~= 1;
                Recover = @(prg,x,y) {x,prg.RecoverDual(y)};
                DistToFace = @(prg,x,y) norm(prg.unreducedPrg.c(:)'-y'*prg.unreducedPrg.A  - ...
                     prg.faces{end}.ProjFace( prg.unreducedPrg.c(:)'-y'*prg.unreducedPrg.A)' );
            end
            
    
           % pars.fid = 0;
          
            [x,y,info] = self.unreducedPrg.Solve();
            infoNoRed = info;
            if CheckFeas(info)
                faceDist = DistToFace(self,x,y);
                
                if ~isfield(info,'dimacs')
                    [~,errors] = self.unreducedPrg.CheckSolution(x,y,1);
                else
                   errors = info.dimacs(:)';
                end
                
                
                data = [errors,faceDist];
                dataNoRed = data; infoNoRed = info;
            end
            
            [x,y,info] = self.Solve();
            if  CheckFeas(info)
                sol = Recover(self,x,y); 
                faceDist =  DistToFace(self,sol{1},sol{2});
                
                if ~isfield(info,'dimacs')
                [~,errors] = self.CheckSolution(x,y,1);
                else
                   errors = info.dimacs(:)';
                end
                data = [errors,faceDist];
         
            
        
                fprintf('Reduced \n')
                data = [data,info.time,nnz(self.A)+nnz(self.c)];
                fprintf([repmat('%2.3E &\t',1,length(data)-2), '%d &\t %d\n'],data);
                datar=data;
                fprintf('Not Reduced\n')
                data = dataNoRed; info = infoNoRed;
                data = [data,info.time,nnz(self.unreducedPrg.A)+nnz(self.unreducedPrg.c) ];
                fprintf([repmat('%2.3E &\t',1,length(data)-2), '%d & \t %d\n'],data); 
                dataNoRed=data;
            end
        end
    
        function PrintSparsity(self)

            nnz_orig = nnz(self.unreducedPrg.A)+nnz(self.unreducedPrg.c);
            nnz_red = nnz(self.A)+nnz(self.c);
            fprintf('nnz Orig\t nnz Red   \n');
            fprintf('%d \t %d ',[nnz_orig,nnz_red]);
            fprintf('\n');
        end
           


        
        function PrintStats(self)
            
            K = self.K;
            Korig = self.unreducedPrg.K;
   
            if  all(K.s == Korig.s)
                display('-------------------------------------------------------------------------')
                display('frlib: No reductions found.')
                display('-------------------------------------------------------------------------')
            else
                display('-------------------------------------------------------------------------')
                display('frlib: reductions found!')
                display('-------------------------------------------------------------------------')
                display(['  Dim PSD constraint(s) (original):  ',sprintf('%d ',Korig.s)])
                display(['  Dim PSD constraint(s) (reduced):   ',sprintf('%d ',K.s)])              
                if  self.opts.useQR == 0 && size(self.A,2) < size(self.A,1)
                   display(' ')
                   display('  Warning: Reduced problem has more equations than primal variables. This may ')
                   display('  cause solver errors.  Perform reductions with opts.useQR = 1 to remove');
                   display('  linearly dependent equations.');
                end   
            end

        end

        function y = noReductions(self)
            y = ~self.faces{end}.isProper;
        end

    end

end
