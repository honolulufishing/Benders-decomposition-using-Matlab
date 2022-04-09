
function [OptX,OptY,OptValue,k] = General_BD_socp(C,D,A,B,b,E,h,F,r,Q,g,l)
% generalized benders decomposition (GBD) matlab-based programming code
% Author: Chao Lei with The Hong Kong Polytechnic University
% What kind of problems can hire this GBD algorithm to solve?
% There are two conditions to be satisfied:
% i) Convexity: the problem should be convex on y given the discrete variables x;
% ii) Linear separability: the problem should be linear on x given the continuous variables y.
%  This kind of problems is w.r.t. the following form:
%               min  C*x+D*y
%               s.t. A*x+B*y<=b; x in {0,1},and y>=0
%                    E*y=h;
%                    F*x=r;
%                    y'*Q*y+l'*y<=g
% Reference:[1] Lee, Mengyuan, et al. "Accelerating generalized benders decomposition for wireless resource allocation," IEEE Trans.on Wireless Communications, vol.20, no.2, pp.1233-1247, 2020.

epsilon = 1e-4;                % stopping criteria of GBD algorithm when abs(LB-UB) is less than epsilon
[nRowA,nColA] = size(A);
[nRowB,nColB] = size(B);
[nRowE,nColE] = size(E);
[nRowF,nColF] = size(F);
n_x = size(C,1);n_y = size(D,1);
n_socp = size(Q,2);

% initialization for GBD
x0 = zeros(nColA,1);            % initial points for 0-1 integer variables
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
while k<kmax
    
    % -----step 1:solve sub-problem (SOCP-based model only with y)-------%
    % convert SOCP-based model to a mosek standard representation
    param.MSK_IPAR_LOG = 0;   % all log information is suppressed in mosek output table
    param.MSK_IPAR_OPF_WRITE_HEADER = 'MSK_OFF';
    
    % Specify the non-conic part of the problem.
    prob.c = [D' zeros(1,(n_y+1)*n_socp) ];
    a_socp = zeros(n_y*n_socp,(n_y+1)*n_socp + n_y);
    for j = 1 : n_socp
        a_socp(n_y*(j-1)+1:n_y*j,1:n_y) = -sqrt(Q{j})   ;
        a_socp(n_y*(j-1)+1:n_y*j,n_y+(n_y+1)*(j-1)+1:n_y+(n_y+1)*j) =  [ zeros(n_y,1) eye(n_y,n_y) ];
    end
    prob.a = sparse([B zeros(nRowB,(n_y+1)*n_socp) ;
        E zeros(nRowE, (n_y+1)*n_socp  );
        a_socp]);
    blc_socp =[];
    for j = 1 : n_socp
        blc_socp = [   blc_socp ;0.5.*( inv(sqrt(Q{j})) )'*l(:,j)];
    end
    prob.blc = [-inf.*ones(size(b,1),1);h; blc_socp];
    prob.buc = [b - A*x0 ;h; blc_socp];
    blx_socp = [];
    for j = 1 : n_socp
        blx_socp = [ blx_socp  sqrt(0.25.*l(:,j)'*Q{j}'*l(:,j) + g(1,j))  -inf.*ones(1,n_y) ];
    end
    prob.blx = [zeros(n_y,1)',   blx_socp ];
    bux_socp = [];
    for j = 1 : n_socp
        bux_socp = [ bux_socp  sqrt(0.25.*l(:,j)'*Q{j}'*l(:,j) + g(1,j))  inf.*ones(1,n_y) ];
    end
    prob.bux = [inf*ones(n_y,1)',bux_socp];
    
    % Specify the number of cones.
    prob.cones = cell(n_socp ,1);
    for j = 1 : n_socp
        prob.cones{j}.type = 'MSK_CT_QUAD';      % MSK_CT_QUAD for Rotated Qudratic Cone
        prob.cones{j}.sub = [n_y + (n_y+1)*(j-1)+1, n_y + 2 + (n_y+1)*(j-1):n_y + n_y + 1 + (n_y+1)*(j-1)];
    end
    
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
        
       ff = 0;
        for j = 1 : n_socp
            ff =  ff + res.sol.itr.snx( (j-1)*(n_y+1)+n_y+1,1)*sqrt( g(1,j) - (yf0(p,:)*Q{j}*yf0(p,:)' + l(:,j)'*yf0(p,:)'));
        end
   
        bZX_fea(p,:)=  -(  D'*yf0(p,1:n_y)' +  ...
            u(p,1:nRowB)*(B*yf0(p,1:n_y)'-b) + u(p,1+nRowB:nRowB+nRowE)*(E*yf0(p,1:n_y)'-h) ...
            + u(p,1+nRowB+nRowE:nRowB+2*nRowE)*(-E*yf0(p,:)'+h) -  ff ) ;
        
        % infeasible: if y in sub-problem is infeasible, then solve the realxed sub-problem model
    elseif strcmp(res.sol.itr.solsta, 'PRIMAL_INFEASIBLE_CER')
        
        % convert SOCP-based model to a mosek standard representation
        param.MSK_IPAR_LOG = 0;   % all log information is suppressed in mosek output table
        
        % Specify the non-conic part of the problem.
        prob.c = [1 zeros(1,n_y) zeros(1,(n_y+1)*n_socp) ];
        prob.a = sparse([-ones(nRowB,1) B zeros(nRowB,(n_y+1)*n_socp) ;
            zeros(nRowE,1) E zeros(nRowE, (n_y+1)*n_socp  );
            zeros(n_y*n_socp,1) a_socp]);
        prob.blc = [-inf.*ones(size(b,1),1);h; blc_socp];
        prob.buc = [b - A*x0 ;h; blc_socp];
        prob.blx = [zeros(n_y + 1,1)',   blx_socp ];
        prob.bux = [inf*ones(n_y + 1,1)',bux_socp];
        
        % Specify the number of cones.
        prob.cones = cell(n_socp ,1);
        for j = 1 : n_socp
            prob.cones{j}.type = 'MSK_CT_QUAD';      % MSK_CT_QUAD for Rotated Qudratic Cone
            prob.cones{j}.sub = [n_y + (n_y+1)*(j-1)+2, n_y+1+2+(n_y+1)*(j-1):n_y+n_y+2+(n_y+1)*(j-1)];
        end
        
        [~,res] = mosekopt('minimize echo(0)',prob,param);
        
        q = q+1;
        y0(q,:) = res.sol.itr.xx(2:n_y + 1,1)';
        s = res.sol.itr.xx(1,1);
        
        % Dual variables
        v(q,:) = [res.sol.itr.suc(1:nRowB+nRowE,1);     %  for upper bounds of linear constraints
            res.sol.itr.slc(1+nRowB:nRowB+nRowE,1);     % for lower bounds of linear constraints
            ]';
         
        dd = 0;
        for j = 1 : n_socp
            dd =  dd + res.sol.itr.snx( j*(n_y+1)+1,1)*sqrt( g(1,j) - (y0(q,:)*Q{j}*y0(q,:)' + l(:,j)'*y0(q,:)')) ...
                + res.sol.itr.snx( j*(n_y+1)+2:j*(n_y+1)+n_y+1 ,1)'*( 0.5.*( inv(sqrt(Q{j})) )'*l(:,j) + sqrt(Q{j})*y0(q,:)' -  res.sol.itr.xx(j*(n_y+1)+2:j*(n_y+1)+n_y+1,1) );
        end


        bZX_inf(q,:)= (  -v(q,1:nRowB)*(B*y0(q,1:n_y)'-b-s) - v(q,1+nRowB:nRowB+nRowE)*(E*y0(q,1:n_y)'-h-s) - v(q,1+nRowB+nRowE:nRowB+2*nRowE)*(-E*y0(q,1:size(D,1))'+h+s) ...
            -  dd);
        
    end
    
    % -----------step2: solve relaxed master model-------------%
    CoefZ = 1;
    f = [CoefZ,zeros(1,n_x)];                                        % coefficients of objective function
    if p>0
        CoefMatZX = [zeros(nRowF,1) F
            zeros(nRowF,1) -F
            -repmat(CoefZ,p,1)  repmat(C',p,1)+ u(1:p,1:nRowB)*A  ;   % optimal cuts
            zeros(q,1) v(1:q,1:nRowB)*A                               % feasible cuts
            ];
    else
        CoefMatZX = [zeros(nRowF,1) F
            zeros(nRowF,1) -F
            zeros(q,1) v(1:q,1:nRowB)*A                                % feasible cuts
            ];
    end
    
    bZX = [ r;
        -r;
        bZX_fea(1:p,1);
        bZX_inf(1:q,1)
        ];
    
    vlb = zeros(1,n_x + 1 );
    vub = [ inf, ones(1,n_x) ];
    ctype = char([ repmat({'C'},1,1) repmat({'I'},1,size(C,1))])';
    [OptZX,minZ,ExitflagBint] = cplexmilp(f,CoefMatZX,bZX, [], [],[ ], [ ], [ ], vlb, vub, ctype, [ ],options_cplex );
    
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

% convert SOCP-based model to a mosek standard representation
param.MSK_IPAR_LOG = 0;   % all log information is suppressed in mosek output table
param.MSK_IPAR_OPF_WRITE_HEADER = 'MSK_OFF';

% Specify the non-conic part of the problem.
prob.c = [D' zeros(1,(n_y+1)*n_socp) ];
a_socp = zeros(n_y*n_socp,(n_y+1)*n_socp + n_y);
for j = 1 : n_socp
    a_socp(n_y*(j-1)+1:n_y*j,1:n_y) = -sqrt(Q{j})   ;
    a_socp(n_y*(j-1)+1:n_y*j,n_y+(n_y+1)*(j-1)+1:n_y+(n_y+1)*j) =  [ zeros(n_y,1) eye(n_y,n_y) ];
end
prob.a = sparse([B zeros(nRowB,(n_y+1)*n_socp) ;
    E zeros(nRowE, (n_y+1)*n_socp  );
    a_socp]);
blc_socp =[];
for j = 1 : n_socp
    blc_socp = [   blc_socp ;0.5.*( inv(sqrt(Q{j})) )'*l(:,j)];
end
prob.blc = [-inf.*ones(size(b,1),1);h; blc_socp];
prob.buc = [b - A*x0 ;h; blc_socp];
blx_socp = [];
for j = 1 : n_socp
    blx_socp = [ blx_socp  sqrt(0.25.*l(:,j)'*Q{j}'*l(:,j) + g(1,j))  -inf.*ones(1,n_y) ];
end
prob.blx = [zeros(n_y,1)',   blx_socp ];
bux_socp = [];
for j = 1 : n_socp
    bux_socp = [ bux_socp  sqrt(0.25.*l(:,j)'*Q{j}'*l(:,j) + g(1,j))  inf.*ones(1,n_y) ];
end
prob.bux = [inf*ones(n_y,1)',bux_socp];
% Specify the number of cones.
prob.cones = cell(n_socp ,1);
for j = 1 : n_socp
    prob.cones{j}.type = 'MSK_CT_QUAD';      % MSK_CT_QUAD for Rotated Qudratic Cone
    prob.cones{j}.sub = [n_y + (n_y+1)*(j-1)+1, n_y + 2 + (n_y+1)*(j-1):n_y + n_y + 1 + (n_y+1)*(j-1)];
end
[r1,res]=mosekopt('minimize echo(0)',prob,param);
OptY = res.sol.itr.xx(1:n_y,1);
OptValue = C'*OptX + D'*OptY;