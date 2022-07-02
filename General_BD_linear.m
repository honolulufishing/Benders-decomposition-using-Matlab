
function [OptX,OptY,OptValue,UB_solution,LB_solution,k] = General_BD_linear(H,C,D,A,B,b,E,h,F,r_le,G,r_ls,lb_y,ub_y)
% Generalized benders decomposition (GBD) matlab-based programming code
% Author: Chao Lei with The Hong Kong Polytechnic University
%  This kind of problems is w.r.t. the following form:
%               min  C*x+D*y+y'*H*y
%               s.t. A*x+B*y<=b; x in {0,1},and lb_y <= y <= ub_y
%                    E*y=h;
%                    F*x<=r_le;
%                    G*x=r_ls;
% Reference:[1] Geoffrion, Arthur M, "Generalized benders decomposition," Journal of optimization theory and applications, vol. 10, no. 4: pp. 237-260, 1972

epsilon = 1e-4;                % stopping criteria of GBD algorithm when abs(LB-UB) is less than epsilon
[nRowA,nColA] = size(A);
[nRowB,nColB] = size(B);
[nRowE,nColE] = size(E);
[nRowF,nColF] = size(F);
[nRowG,nColG] = size(G);
n_x = size(C,1);n_y = size(D,1);

% initialization for GBD
x0 = zeros(nColA,1);            % initial points for 0-1 integer variables
LB = -1.0e18;
UB = 50;
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
UB0 = zeros(kmax,1);
LB0 = zeros(kmax,1);
while k<kmax
    
    % -----step 1:solve sub-problem (SOCP-based model only with y)-------%
    % convert SOCP-based model to a mosek standard representation
    param.MSK_IPAR_LOG = 0;   % all log information is suppressed in mosek output table
    param.MSK_IPAR_OPF_WRITE_HEADER = 'MSK_OFF';
    
    % Specify the non-conic part of the problem.
    prob.c = [D' ];
    
    prob.qosubi = [1:n_y]';
    prob.qosubj = [1:n_y ]';
    prob.qoval = 2.*H';
    
    prob.a = sparse([B  ;
        E ]);
    prob.blc = [-inf.*ones(size(b,1),1);h];
    prob.buc = [b - A*x0 ;h];
    prob.blx = lb_y';
    prob.bux = ub_y';
    
    [~,res] = mosekopt('minimize echo(0)',prob,param);
    
    % feasible: attain the optimal solution from sub-problem
    if strcmp(res.sol.itr.solsta, 'OPTIMAL')
        
        p = p+1;
        % Display the primal solution.
        yf0(p,:) = res.sol.itr.xx(1:n_y,1)';
        
        % Dual variables for upper bounds of linear constraints
        u(p,:) = [res.sol.itr.suc(1:nRowB+nRowE,1);    %  for upper bounds of linear constraints
            res.sol.itr.slc(1+nRowB:nRowB+nRowE,1);    % for lower bounds of linear constraints
            ]';
        
        UB = D'*yf0(p,:)' + C'*x0;
        bZX_fea(p,:)=  -(  D'*yf0(p,1:n_y)' +  ...
            u(p,1:nRowB)*(B*yf0(p,1:n_y)'-b) + u(p,1+nRowB:nRowB+nRowE)*(E*yf0(p,1:n_y)'-h) ...
            + u(p,1+nRowB+nRowE:nRowB+2*nRowE)*(-E*yf0(p,:)'+h)  ) ;
        
        % infeasible: if y in sub-problem is infeasible, then solve the realxed sub-problem model
    elseif strcmp(res.sol.itr.solsta, 'PRIMAL_INFEASIBLE_CER')
        % convert SOCP-based model to a mosek standard representation
        param.MSK_IPAR_LOG = 0;   % all log information is suppressed in mosek output table
        
        % Specify the non-conic part of the problem.
        prob.c = [1 zeros(1,n_y)  ];
        prob.a = sparse([-ones(nRowB,1) B  ;
            zeros(nRowE,1) E ]);
        prob.blc = [-inf.*ones(size(b,1),1);h];
        prob.buc = [b - A*x0 ;h];
        prob.blx = [ 0, lb_y'  ];
        prob.bux = [inf,ub_y'];
        
        [~,res] = mosekopt('minimize echo(0)',prob,param);
        
        q = q+1;
        y0(q,:) = res.sol.itr.xx(2:n_y + 1,1)';
        s = res.sol.itr.xx(1,1);
        
        % Dual variables
        v(q,:) = [res.sol.itr.suc(1:nRowB+nRowE,1);     %  for upper bounds of linear constraints
            res.sol.itr.slc(1+nRowB:nRowB+nRowE,1);     % for lower bounds of linear constraints
            ]';
        
        bZX_inf(q,:)= (  -v(q,1:nRowB)*(B*y0(q,1:n_y)'-b) - v(q,1+nRowB:nRowB+nRowE)*(E*y0(q,1:n_y)'-h) - v(q,1+nRowB+nRowE:nRowB+2*nRowE)*(-E*y0(q,1:size(D,1))'+h) ...
            );
        
    end
    
    % -----------step2: solve relaxed master model-------------%
    CoefZ = 1;
    f = [CoefZ,zeros(1,n_x)];                                        % coefficients of objective function
    if p>0
        CoefMatZX = [zeros(nRowF,1) F
            -repmat(CoefZ,p,1)  u(1:p,1:nRowB)*A  ;    % optimality cuts
            zeros(q,1) v(1:q,1:nRowB)*A                                % feasibility cuts
            ];
    else
        CoefMatZX = [zeros(nRowF,1) F
            zeros(q,1) v(1:q,1:nRowB)*A                                % feasibility cuts
            ];
    end
    
    bZX = [ r_le;
        repmat(C'*x0,p,1)+ bZX_fea(1:p,1);
        bZX_inf(1:q,1)
        ];
    
    vlb = zeros(1,n_x + 1 );
    vub = [ inf, ones(1,n_x) ];
    ctype = char([ repmat({'C'},1,1) repmat({'I'},1,size(C,1))])';
    [OptZX,minZ,ExitflagBint] = cplexmilp(f,CoefMatZX,bZX, [zeros(nRowG,1) G], r_ls,[ ], [ ], [ ], vlb, vub, ctype, [ ],options_cplex );
    
    if ExitflagBint ==1
        LB = minZ ;
        x0 = OptZX(2:end);
    else
        fprintf('Error:iterations of GBD algorithm are terminated.');
        break;
    end
    
    LB0(k,1) = LB;
    UB0(k,1) = UB;
    
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

UB_solution = [UB0(2,1)+20;UB0(2:k,1)];
LB_solution = LB0(1:k,1);
figure(1)
plot(1:k, [UB0(2,1)+20;UB0(2:k,1)],'r.-',1:k, LB0(1:k,1),'b.-')
legend('upper bound','lower bound');