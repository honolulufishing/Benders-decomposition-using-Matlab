
function [sub_solution_flag,UB_y, u, bZX_fea, v, bZX_inf] = mgbd_sp(x0,D_block,A,B_block,b, E_block ,h_block,Q_block,g_block,l_block ,n_y, start_block)

[nRowA_block,nColA] = size(A);
        [nRowB_block,nColB] = size(B_block);
        [nRowE_block,nColE] = size(E_block);
        
        param.MSK_IPAR_LOG = 0;   % all log information is suppressed in mosek output table
        param.MSK_IPAR_OPF_WRITE_HEADER = 'MSK_OFF';
        
        % Specify the non-conic part of the problem.
        prob.c = [D_block' zeros(1,n_y+1) ];
        a_socp = zeros(n_y,n_y+1+ n_y);
        
        a_socp(1:n_y,1:n_y) = -sqrt(Q_block)   ;
        a_socp(1:n_y,n_y+1:n_y+n_y+1) =  [ zeros(n_y,1) eye(n_y,n_y) ];
        
        prob.a = sparse([B_block zeros(nRowB_block,n_y+1) ;
            E_block zeros(nRowE_block, n_y+1);
            a_socp]);
        blc_socp = 0.5.*( inv(sqrt(Q_block)) )'*l_block;
        bs = b - A*x0;
        prob.buc = [bs((start_block-1)*n_y+1:start_block*n_y,1) ;h_block; blc_socp];
        prob.blc = [-inf.*ones(2,1);h_block; blc_socp];  % 2 changed
        blx_socp = [sqrt(0.25.*l_block'*Q_block'*l_block + g_block)  -inf.*ones(1,n_y)] ;
        
        prob.blx = [zeros(n_y,1)',   blx_socp ];
        bux_socp =  [sqrt(0.25.*l_block'*Q_block'*l_block + g_block)  inf.*ones(1,n_y)];
        prob.bux = [inf*ones(n_y,1)',bux_socp];
        
        % Specify the number of cones.
        prob.cones = cell(1,1);
        prob.cones{1}.type = 'MSK_CT_QUAD';      % MSK_CT_QUAD for Rotated Qudratic Cone
        prob.cones{1}.sub = [n_y+1, n_y + 2 :n_y + n_y + 1 ];
        
        [~,res] = mosekopt('minimize echo(0)',prob,param);
        
        sub_solution_flag = 0;
        % feasible: attain the optimal solution from sub-problem
        if strcmp(res.sol.itr.solsta , 'OPTIMAL')
            
            sub_solution_flag = 1;
            % Display the primal solution.
            y0 = res.sol.itr.xx(1:n_y,1)';
            
            % Dual variables for upper bounds of linear constraints
            u = [res.sol.itr.suc(1:nRowB_block+nRowE_block,1);    %  for upper bounds of linear constraints
                res.sol.itr.slc(1+nRowB_block:nRowB_block+nRowE_block,1);    % for lower bounds of linear constraints
                ]';
            
            UB_y = D_block'*y0' ;
            
            ff =  res.sol.itr.snx(n_y+1,1)*( g_block - (y0*Q_block*y0' + l_block'*y0'));
            
            
            bZX_fea =  -(  D_block'*y0(1,:)' +  ...
                u(1,1:nRowB_block)*(B_block * y0(1,:)'-b((start_block-1)*n_y+1:start_block*n_y,1)) + u(1,1+nRowB_block:nRowB_block+nRowE_block)*( E_block *y0(1,1:n_y)'-h_block) ...
                + u(1,1+nRowB_block+nRowE_block:nRowB_block+2*nRowE_block)*(-E_block*y0' + h_block ) -  ff );
            
            v = [];
            bZX_inf = [];
            
            % infeasible: if y in sub-problem is infeasible, then solve the realxed sub-problem model
        elseif strcmp(res.sol.itr.solsta, 'PRIMAL_INFEASIBLE_CER')
            
            sub_solution_flag = 2;
            
            % convert SOCP-based model to a mosek standard representation
            param.MSK_IPAR_LOG = 0;   % all log information is suppressed in mosek output table
            
            % Specify the non-conic part of the problem.
            prob.c = [1 zeros(1,n_y) zeros(1,n_y+1) ];
            prob.a = sparse([-ones(nRowB_block,1) B_block zeros(nRowB_block,n_y+1) ;
                zeros(nRowE_block,1) E_block zeros(nRowE_block, n_y+1 );
                zeros(n_y,1) a_socp]);
            
            bs = b - A*x0;
            prob.buc = [bs((start_block-1)*n_y+1:start_block*n_y,1) ;h_block; blc_socp];
            prob.blc = [-inf.*ones(2,1);h_block; blc_socp];
            
            prob.blx = [zeros(n_y + 1,1)',   blx_socp ];
            prob.bux = [inf*ones(n_y + 1,1)',bux_socp];
            
            % Specify the number of cones.
            prob.cones = cell(1 ,1);
            prob.cones{1}.type = 'MSK_CT_QUAD';      % MSK_CT_QUAD for Rotated Qudratic Cone
            prob.cones{1}.sub = [n_y+2, n_y+1+2:n_y+n_y+2];
            
            [~,res] = mosekopt('minimize echo(0)',prob,param);
            
            yf0 = res.sol.itr.xx(2:n_y + 1,1)';
            %s = res.sol.itr.xx(1,1);
            
            % Dual variables
            v = [res.sol.itr.suc(1:nRowB_block+nRowE_block,1);     %  for upper bounds of linear constraints
                res.sol.itr.slc(1+nRowB_block:nRowB_block+nRowE_block,1);     % for lower bounds of linear constraints
                ]';
            
            dd =  res.sol.itr.snx( n_y+1+1,1)*( g_block - (yf0(1,:)*Q_block*yf0(1,:)' + l_block'*yf0(1,:)') )   ...
                + res.sol.itr.snx( n_y+1+2:n_y+1+n_y+1 ,1)'*( 0.5.*( inv(sqrt(Q_block)) )'*l_block + sqrt(Q_block)*yf0(1,:)' -  res.sol.itr.xx(n_y+1+2:n_y+1+n_y+1,1) );
            
            bZX_inf = (  -v(1,1:nRowB_block)*(B_block*yf0(1,1:n_y)'-b((start_block-1)*n_y+1:start_block*n_y,1)) - v(1,1+nRowB_block:nRowB_block+nRowE_block)*(E_block*yf0(1,1:n_y)'-h_block) - v(1,1+nRowB_block+nRowE_block:nRowB_block+2*nRowE_block)*(-E_block*yf0(1,1:n_y)'+h_block) ...
                -  dd);
            
            UB_y = inf;
            u = [];
            bZX_fea = [];
        end