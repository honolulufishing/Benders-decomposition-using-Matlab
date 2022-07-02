
function [OptX,OptY,OptValue,k,delta_T,delta_T_MP] = General_BD_socp(H,C,D,A,B,b,E,h,F,r_le,G,r_ls,Qij,Q_diag,g,l,lb_y,ub_y)
% Generalized benders decomposition (GBD) matlab-based programming code
% Author: Chao Lei with The Hong Kong Polytechnic University
% What kind of problems can hire this GBD algorithm to solve?
% There are two conditions to be satisfied:
% i) Convexity: the problem should be convex on y given the discrete variables x;
% ii) Linear separability: the problem should be linear on x given the continuous variables y.
%  This kind of problems is w.r.t. the following form:
%               min  C*x+D*y + y'*H*y
%               s.t. A*x+B*y<=b; x in {0,1},and lb_y <= y <= ub_y
%                    E*y=h;
%                    F*x<=r_le;
%                    G*x=r_ls;
%                    y_i'*Q*y_i + l'*y_i + g <= y_j  (j>i)
% Reference: [1] Geoffrion, Arthur M, "Generalized benders decomposition," Journal of optimization theory and applications, vol. 10, no. 4: pp. 237-260, 1972

epsilon = 1e-4;                % stopping criteria of GBD algorithm when abs(LB-UB) is less than epsilon
[nRowA,nColA] = size(A);
[nRowB,nColB] = size(B);
[nRowE,nColE] = size(E);
[nRowF,nColF] = size(F);
[nRowG,nColG] = size(G);
n_x = size(C,1);n_y = size(D,1);
n_socp = size(Qij,2);
n_Q_element = 2;   % new added variables for each SOC constraints are n_Q_element+1

% initialization for GBD
 x0 = zeros(nColA,1);            % initial points for 0-1 integer variables
% x0=[   1.0000
%          0
%     1.0000
%     0.0000
%     0.0000];


LB = -1.0e18;
UB = inf;
p = 0;
q = 0;
options_cplex = cplexoptimset;  % cplex solver parameters
options_cplex.Display = 'off';
kmax = 10;
k = 1;
u = zeros(kmax,nRowB+2*nRowE);
v = zeros(kmax,nRowB+2*nRowE);
yf0 = zeros(kmax,size(D,1));
y0 = zeros(kmax,size(D,1));
bZX_fea = zeros(kmax,1);
bZX_inf = zeros(kmax,1);
delta_T = zeros(kmax,1);
delta_T_MP = zeros(kmax,1);
while k<kmax
    
    t11=clock;
    % -----step 1:solve sub-problem (SOCP-based model only with y)-------%
    % convert SOCP-based model to a mosek standard representation
    param.MSK_IPAR_LOG = 0;   % all log information is suppressed in mosek output table
    param.MSK_IPAR_OPF_WRITE_HEADER = 'MSK_OFF';

    % Specify the non-conic part of the problem.
    prob.c = [D'  zeros(1,(n_Q_element+1)*n_socp) ];

    a_socp_1 = zeros(n_socp,(n_Q_element+1)*n_socp + n_y );
    for j = 1 : n_socp
        a_socp_1( j , [Qij(j,1), n_y + j ]) = [ sqrt(Q_diag(j,1)) -1]   ;
    end
    
    a_socp_2 = zeros(n_socp,(n_Q_element+1)*n_socp + n_y );
    for j = 1 : n_socp
        a_socp_2( j , [Qij(j,2), n_y+n_socp+j]) = [1 -1]  ;
    end
    
    prob.a = sparse([B zeros(nRowB,(n_Q_element+1)*n_socp)  ;
        E zeros(nRowE,(n_Q_element+1)*n_socp);
        a_socp_1;  a_socp_2]);
    
    blc_socp_1 = zeros(n_socp,1);
    blc_socp_2 =zeros(n_socp,1);
    for j = 1 : n_socp
        blc_socp_1(j,1) = [ -0.5.*l(j,1).*(1./sqrt(Q_diag(j,1)))'   ];
        blc_socp_2(j,1)= -0.25*Q_diag(j,1)*l(j,1)^2 + g(j,1);
    end
    prob.blc = [-inf.*ones(size(b,1),1);h; blc_socp_1;blc_socp_2];
    prob.buc = [b - A*x0 ;h; blc_socp_1;blc_socp_2];

    prob.blx = [lb_y', -inf.*ones(1,n_socp), zeros(1,n_socp), 0.5*ones(1,n_socp) ];
    prob.bux = [ub_y',inf.*ones(1,n_socp), inf.*ones(1,n_socp), 0.5*ones(1,n_socp)  ];
    
    % Specify the number of cones.
    prob.cones = cell(n_socp ,1);
    for j = 1 : n_socp
        prob.cones{j}.type = 'MSK_CT_RQUAD';      % MSK_CT_QUAD for Rotated Qudratic Cone
        prob.cones{j}.sub = [  n_y+n_socp+j ,  n_y+2*n_socp+j,  n_y+j     ];
    end
    
    [~,res] = mosekopt('minimize echo(0)',prob,param);
    
    
    % feasible: attain the optimal solution from sub-problem
    if strcmp(res.sol.itr.solsta, 'OPTIMAL')
        
        p = p+1;
        % Display the primal solution.
        yf0(p,:) = res.sol.itr.xx(1:n_y,1)';
        J = res.sol.itr.xx(n_y+1:n_y+n_socp,1);
        M = res.sol.itr.xx(n_y+n_socp+1:n_y+2*n_socp,1);
        Tao = res.sol.itr.xx(n_y+2*n_socp+1:n_y+3*n_socp,1);
        
        % Dual variables for upper bounds of linear constraints
        u(p,:) = [res.sol.itr.suc(1:nRowB+nRowE,1);    %  for upper bounds of linear constraints
            res.sol.itr.slc(1+nRowB:nRowB+nRowE,1);    % for lower bounds of linear constraints
            ]';
        
        for j = 1: nRowB+nRowE
            if u(p,j) < 0.001
                u(p,j) = 0;
            end
        end
        UB = D'*yf0(p,:)' + C'*x0;

        % res.sol.itr.snx = [y1~y4 J1~Jn_socp m1~m_socp tao_1~tao_socp]
        ff =  res.sol.itr.snx( n_y+1:n_y+n_socp  ,1)' * J ...
            + res.sol.itr.snx(  n_y+n_socp+1:n_y+2*n_socp ,1)'* M + res.sol.itr.snx(  n_y+2*n_socp+1:n_y+3*n_socp ,1)'* Tao;
        
        bZX_fea(p,:)=  -(  D'*yf0(p,1:n_y)' +  ...
            u(p,1:nRowB)*(B*yf0(p,1:n_y)'-b) + u(p,1+nRowB:nRowB+nRowE)*(E*yf0(p,1:n_y)'-h) ...
            + u(p,1+nRowB+nRowE:nRowB+2*nRowE)*(-E*yf0(p,:)'+h) -  ff ) ;
        
        % infeasible: if y in sub-problem is infeasible, then solve the realxed sub-problem model
    elseif strcmp(res.sol.itr.solsta, 'PRIMAL_INFEASIBLE_CER')
        
        % convert SOCP-based model to a mosek standard representation
        param.MSK_IPAR_LOG = 0;   % all log information is suppressed in mosek output table
        
        % Specify the non-conic part of the problem.
        prob.c = [1 zeros(1,n_y)    zeros(1,(n_Q_element+1)*n_socp )  ];
        
        prob.a = sparse([-ones(nRowB,1) B zeros(nRowB,(n_Q_element+1)*n_socp)  ;
            zeros(nRowE,1) E zeros(nRowE,(n_Q_element+1)*n_socp);
            zeros(n_socp,1) a_socp_1;
            zeros(n_socp,1)  a_socp_2]);
        
        prob.blc = [-inf.*ones(size(b,1),1);h; blc_socp_1;blc_socp_2];
        prob.buc = [b - A*x0 ;h; blc_socp_1;blc_socp_2];
        
        
        prob.blx = [0,lb_y', -inf.*ones(1,n_socp),  zeros(1,n_socp), 0.5*ones(1,n_socp) ];
        prob.bux = [inf,ub_y',inf.*ones(1,n_socp), inf.*ones(1,n_socp), 0.5*ones(1,n_socp)  ];
        
        % Specify the number of cones.
        prob.cones = cell(n_socp ,1);
        for j = 1 : n_socp
            prob.cones{j}.type = 'MSK_CT_RQUAD';      % MSK_CT_QUAD for Rotated Qudratic Cone
            prob.cones{j}.sub = [  1+n_y+n_socp+j,  1+n_y+2*n_socp+j,  1+n_y+j     ];
        end
        
        [~,res] = mosekopt('minimize echo(0)',prob,param);
        
        q = q+1;
        y0(q,:) = res.sol.itr.xx(2:n_y + 1,1)';        % s = res.sol.itr.xx(1,1);
        J = res.sol.itr.xx(n_y+1+1:n_y+1+n_socp,1);
        M = res.sol.itr.xx(n_y+1+n_socp+1:n_y+1+2*n_socp,1);
        Tao = res.sol.itr.xx(n_y+1+2*n_socp+1:n_y+1+3*n_socp,1);
        
        % Dual variables
        v(q,:) = [res.sol.itr.suc(1:nRowB+nRowE,1);     %  for upper bounds of linear constraints
            res.sol.itr.slc(1+nRowB:nRowB+nRowE,1);     % for lower bounds of linear constraints
            ]';
        
        for j = 1: nRowB+nRowE
            if v(q,j) < 0.001
                v(q,j) = 0;
            end
        end
        
        % res.sol.itr.snx = [s y1~y4 J1~Jn_socp m1~m_socp tao_1~tao_socp]
        dd =  res.sol.itr.snx( n_y+1+1:n_y+1+n_socp  ,1)' * J ...
            + res.sol.itr.snx(  n_y+1+n_socp+1:n_y+1+2*n_socp ,1)'* M + res.sol.itr.snx(  n_y+1+2*n_socp+1:n_y+1+3*n_socp ,1)'* Tao;

        bZX_inf(q,:)= (  -v(q,1:nRowB)*(B*y0(q,1:n_y)'-b) - v(q,1+nRowB:nRowB+nRowE)*(E*y0(q,1:n_y)'-h) - v(q,1+nRowB+nRowE:nRowB+2*nRowE)*(-E*y0(q,1:size(D,1))'+h) ...
            -  dd);
    end
    t22 = clock;
    delta_T(k,1) = etime(t22,t11);
    
    % -----------step2: solve relaxed master model-------------%
    CoefZ = 1;
    f = [CoefZ,zeros(1,n_x)];                                        % coefficients of objective function
    if p>0
        CoefMatZX = full([zeros(nRowF,1) F
            -repmat(CoefZ,p,1)  repmat(C',p,1)+ u(1:p,1:nRowB)*A  ;   % optimality cuts
            zeros(q,1) v(1:q,1:nRowB)*A                               % feasibility cuts
            ]);
    else
        CoefMatZX = full([zeros(nRowF,1) F
            zeros(q,1) v(1:q,1:nRowB)*A                               % feasibility cuts
            ]);
    end
    
    bZX = [ r_le;
        bZX_fea(1:p,1);
        bZX_inf(1:q,1)
        ];
    
    vlb = zeros(1,n_x + 1 );
    vub = [ inf, ones(1,n_x) ];
    ctype = char([ repmat({'C'},1,1) repmat({'I'},1,size(C,1))])';
    t1=clock;
    [OptZX,minZ,ExitflagBint] = cplexmilp(f,CoefMatZX,bZX, sparse([zeros(nRowG,1) G]), r_ls,[ ], [ ], [ ], vlb, vub, ctype, [ ],options_cplex );
    t2=clock;
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
OptY = yf0(p,:)';
OptValue = C'*OptX + D'*OptY;