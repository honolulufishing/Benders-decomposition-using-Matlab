function [OptX,OptY,OptValue,k,delta_T,delta_T_MP] = General_MBD_socp(C,D,A,B,b,E,h,F,r_le,G,r_ls,Qs,gs,ls,n_block,n_y,lb_y,ub_y)
% modified generalized benders decomposition (Modified GBD) matlab-based programming code
% Author: Chao Lei with The Hong Kong Polytechnic University
% What kind of problems can hire this Modified GBD algorithm to solve?
% There are two conditions to be satisfied:
% i) Convexity: the problem should be convex on y given the discrete variables x;
% ii) Linear separability: the problem should be linear on x given the continuous variables y.
% iii) Linear independence: with given discrete variables x, diffferent groups of continuous variables y are linearly independent.
%  This kind of problems is w.r.t. the following form:
%               min  C*x+D*y
%               s.t. A*x+B*y<=b; x in {0,1},and lb_y <= y <= ub_y
%                    E*y=h;
%                    F*x<=r_le;
%                    G*x=r_ls;
%                    y'*Q*y+l'*y<=g
% Reference:[1] Geoffrion, Arthur M, "Generalized benders decomposition," Journal of optimization theory and applications, vol. 10, no. 4: pp. 237-260, 1972


epsilon = 1e-4;                % stopping criteria of Modified GBD algorithm when abs(LB-UB) is less than epsilon
[nRowF,nColF] = size(F);
[nRowG,nColG] = size(G);
[nRowB,nColB] = size(B);
[nRowE,nColE] = size(E);
n_x = size(C,1);
l = zeros(n_y,n_block);
for j = 1 : n_block
    Q{j} = Qs{j}((j-1)*n_y + 1:j*n_y,(j-1)*n_y + 1:j*n_y);
    l(1:n_y,j) = ls((j-1)*n_y + 1:j*n_y,j);
end
g = gs;

% initialization for Modified GBD
x0 = zeros(n_x,1);            % initial points for 0-1 integer variables
LB = -1.0e18;
UB = inf;
p = 0;
q = 0;
options_cplex = cplexoptimset;  % cplex solver parameters
options_cplex.Display = 'off';
kmax = 10;
k = 1;
U = zeros(kmax,size(B,1));
V = zeros(kmax,size(B,1));
BZX_fea = zeros(kmax,1);
BZX_inf = zeros(kmax,1);
delta_T = zeros(n_block,kmax);
delta_T_MP = zeros(kmax,1);
while k<kmax
    
    % -----step 1:solve sub-problem (SOCP-based model only with y)-------%
    % generate multiple feasibility cuts of many independent sub-problems
    true_sub_solution_flag = zeros(n_block,1);
    U_temp = zeros(1,size(B,1));bZX_fea_temp = 0; UB0 = 0; y_optimal_temp = zeros(1,size(B,1));
    
    for i = 1 : n_block
        
        % convert SOCP-based model to a mosek standard representation
        D_block = D(n_y*(i-1)+1:i*n_y,1);
        B_block = B(n_y*(i-1)+1:i*n_y,n_y*(i-1)+1:i*n_y);
        E_block = E(i,n_y*(i-1)+1:i*n_y);
        h_block = h(i,1);
        
        Q_block = Q{i};
        g_block = g(i,1);
        l_block = l(:,i) ;
        
        t11=clock;
        start_block = i;

        % ------------------ mgbd_sp function ---------%
        [nRowA_block,nColA] = size(A);
        [nRowB_block,nColB] = size(B_block);
        [nRowE_block,nColE] = size(E_block);
        
        param.MSK_IPAR_LOG = 0;   % all log information is suppressed in mosek output table
        param.MSK_IPAR_OPF_WRITE_HEADER = 'MSK_OFF';
        
        % Specify the non-conic part of the problem.
        prob.c = [D_block' zeros(1,n_y+1) ];
        a_socp = zeros(n_y,n_y+1+ n_y);
        
        a_socp(1:n_y,1:n_y) = -sqrt(Q_block)   ;
        a_socp(1:n_y,n_y+1:n_y+n_y+1) =  sparse([ zeros(n_y,1) eye(n_y,n_y) ]);
        
        prob.a = sparse([B_block zeros(nRowB_block,n_y+1) ;
            E_block zeros(nRowE_block, n_y+1);
            a_socp]);
        blc_socp = 0.5.*( inv(sqrt(Q_block)) )'*l_block;
        bs = b - A*x0;
        prob.buc = [bs((start_block-1)*n_y+1:start_block*n_y,1) ;h_block; blc_socp];
        prob.blc = [-inf.*ones(2,1);h_block; blc_socp];  % 2 changed
        blx_socp = [sqrt(0.25.*l_block'*Q_block'*l_block + g_block)  -inf.*ones(1,n_y)] ;
        
        prob.blx = [lb_y',   blx_socp ];
        bux_socp =  [sqrt(0.25.*l_block'*Q_block'*l_block + g_block)  inf.*ones(1,n_y)];
        prob.bux = [ub_y',bux_socp];
        
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
            u = [res.sol.itr.suc(1:nRowB_block+nRowE_block,1);               %  for upper bounds of linear constraints
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
            
            prob.blx = [0, lb_y',   blx_socp ];
            prob.bux = [inf, ub_y',bux_socp];
            
            % Specify the number of cones.
            prob.cones = cell(1 ,1);
            prob.cones{1}.type = 'MSK_CT_QUAD';      % MSK_CT_QUAD for Rotated Qudratic Cone
            prob.cones{1}.sub = [n_y+2, n_y+1+2:n_y+n_y+2];
            
            [~,res] = mosekopt('minimize echo(0)',prob,param);
            
            yf0 = res.sol.itr.xx(2:n_y + 1,1)';
            %s = res.sol.itr.xx(1,1);
            
            % Dual variables
            v = [res.sol.itr.suc(1:nRowB_block+nRowE_block,1);                %  for upper bounds of linear constraints
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
        % ------------------ End of mgbd_sp function ---------%

        if sub_solution_flag == 2      % generate multiple feasibility cuts
            q = q+1;
            V(q,[(start_block-1)*nRowB/n_block+1:(start_block-1)*nRowB/n_block + nRowB/n_block]) = v(1,1:nRowB/n_block);
            BZX_inf(q,1) = bZX_inf;
        elseif sub_solution_flag == 1   % check if all sub-problems have optimality cuts
            true_sub_solution_flag(i,1) = 1;
            U_temp(1,[(start_block-1)*nRowB/n_block+1:(start_block-1)*nRowB/n_block + nRowB/n_block]) = u(1,1:nRowB/n_block);
            bZX_fea_temp = bZX_fea_temp + bZX_fea;
            UB0 = UB0 + UB_y;
            
           y_optimal_temp(1,[(start_block-1)*nRowB/n_block+1:(start_block-1)*nRowB/n_block + nRowB/n_block]) = y0;
        end
        
        t22 = clock;
        delta_T(i,k) = etime(t22,t11);
        
    end
    
    if  all(true_sub_solution_flag) == 1  % generate one optimality cut
        p = p+1;
        U(p,:) = U_temp;
        BZX_fea(p,1) = bZX_fea_temp;
        UB = UB0 + C'*x0;
        OptY = y_optimal_temp';
    end
    
    % -----------step2: solve relaxed master model-------------%
    CoefZ = 1;
    f = [CoefZ,zeros(1,n_x)];                                         % coefficients of objective function
    if p>0 && q > 0
        CoefMatZX = sparse([zeros(nRowF,1) F
            -repmat(CoefZ,p,1)  repmat(C',p,1)+ U(1:p,1:nRowB)*A  ;   % optimality cuts
            zeros(q,1) V(1:q,1:nRowB)*A                               % multiple feasibility cuts
            ]);
    elseif p==0
        CoefMatZX = [zeros(nRowF,1) F
            zeros(q,1)      V(1:q,1:nRowB)*A                          % multiple feasibility cuts
            ];
    elseif q==0
        CoefMatZX = sparse([zeros(nRowF,1) F
            -repmat(CoefZ,p,1)  repmat(C',p,1)+ U(1:p,1:nRowB)*A  ;   % optimality cuts
            ]);
    end
    
    bZX = [ r_le;
        BZX_fea(1:p,1);
        BZX_inf(1:q,1)
        ];
    
    vlb = zeros(1,n_x + 1 );
    vub = [ inf, ones(1,n_x) ];
    ctype = char([ repmat({'C'},1,1) repmat({'I'},1,size(C,1))])';
    
    t1 = clock;
    [OptZX,minZ,ExitflagBint] = cplexmilp(f,CoefMatZX,bZX, sparse([zeros(nRowG,1) G]), r_ls,[ ], [ ], [ ], vlb, vub, ctype, [ ],options_cplex );
    t2 = clock;
    delta_T_MP(k,1) = etime(t2,t1);
    
    
    if ExitflagBint ==1
        LB = minZ ;
        x0 = OptZX(2:end);
    else
        fprintf('Error:iterations of GBD algorithm are terminated.');
        break;
    end
    
    if  abs(LB-UB)<epsilon
        fprintf('Success:GBD algorithm is terminated in abs(LB-UB)/e-3: %f\n',abs(LB-UB)/1000);
        break;
    elseif k==kmax
        fprintf('Warning:the maximum number of iterations was reached.');
        break;
    end
    
    k = k+1;
    
end
OptX = x0;
OptValue = C'*OptX + D'*OptY;
