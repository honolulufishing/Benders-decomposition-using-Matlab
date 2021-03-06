function [OptX,OptY,OptValue,k] = Classic_BD(C,D,A,B,b,E,h,F,r_le,G,r_ls,lb_y,ub_y)
% Classical benders decomposition (CBD) matlab-based programming code
% Author: Chao Lei with The Hong Kong Polytechnic University 
% What kind of problems can hire this CBD algorithm to solve?
% This kind of problems is w.r.t. the following form:
%         min  C*x+D*y
%         s.t. A*x+B*y<=b; x in {0,1},and lb_y <= y <= ub_y
%              E*y=h;
%              F*x<=r_le;
%              G*x=r_ls;

% stopping criteria of GBD algorithm when abs(LB-UB) is less than epsilon
epsilon = 1e-4;
[nRowE,nColE] = size(E);
At = [A;zeros(2*nRowE,size(A,2))];        % convert to [A;0;0]*x+[B;E;-E]*y<=[b;h;-h]
Bt = [B;E;-E];
bt = [b;h;-h];

[nRowA,nColA] = size(At);
[nRowB,nColB] = size(Bt);
% initialization for CBD
x0 = zeros(nColA,1);  % initial points for 0-1 integer variables
LB = -1e10;
UB = inf;
p = 0;
q = 0;
k = 1;
kmax = 10;
u = zeros(kmax,nRowA);
v = zeros(kmax,nRowA);
options_cplex = cplexoptimset;
options_cplex.Display = 'off';
options = optimset('LargeScale', 'off', 'Simplex', 'on');
while k<kmax

    %-----step 1:solve sub-problem (LP-based model only with y)-------%
    [UorV,fval,exitflag] = linprog((-At*x0+bt),-Bt',D',[],[], lb_y, ub_y,[],options);
    if exitflag == 1 % if sub-problem has extreme points
        p = p+1;
        u(p,:) = UorV'; % extreme points
        UB = C'*x0 + u(p,:)*(-bt+At*x0);
    elseif exitflag == -3 % if sub-problem is unbounded
        q = q+1;
        v(q,:) = UorV'/1e16;% extreme rays
        v(q,v(q,:)<0.0001) = 0;
    end

    % -----------step2: solve relaxed master model-------------%
    CoefZ = 1;
    f = [CoefZ,zeros(1,nColA)];                                  % coefficients of objective function
    if p>0
        CoefMatZX = [-ones(p,1), repmat(C',p,1) + u(1:p,:)*At;   % optimality cuts
        zeros(q,1),v(1:q,:)*At ;                                 % feasibility cuts
            ];
    else
        CoefMatZX = [
        zeros(q,1),v(1:q,:)*At ;                                % feasibility cuts     
            ];
    end
    
    bZX = [u(1:p,:)*bt;v(1:q,:)*bt];
    vlb = zeros(1,size(C,1)+1);
    vub=[inf, ones(1,size(C,1))];
    ctype = char([ repmat({'C'},1,1) repmat({'I'},1,size(C,1))])';
    
    
    [OptZX,minZ,ExitflagBint] = cplexmilp(f,CoefMatZX,bZX, [zeros(size(F,1),1),F], r_le,[zeros(size(G,1),1),G], r_ls, [ ], vlb, vub, ctype, [ ],options_cplex );
    
    if ExitflagBint ==1
        x0 = OptZX(2:end);
        LB =  minZ;  
    else
        fprintf('Error:iterations of CBD algorithm are terminated.');
        break;
    end
    
    if  abs(LB-UB)<epsilon
        fprintf('Benders Decomposition Optimization is terminated in abs(LB-UB)/e-3: %f\n',abs(LB-UB)/1000);
        break;
    elseif k==kmax
        fprintf('Warning:the maximum number of iterations was reached.');
        break;
    end
    
    k=k+1;
    
end
OptX = x0;
[OptY,OptValue] = linprog(D',B,-A*OptX+b,E,h,zeros(nColB,1),inf(nColB,1),[],options);
OptValue = C'*OptX+D'*OptY;
