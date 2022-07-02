
clc;clear all

% Example 1:a socp problem
Num = 100;
Ct = [7 7 7 7 7]';
Dt = [1 1 1 1]';
for i = 1: Num
    C(5*(i-1)+1:5*i,1) = Ct;
    D(4*(i-1)+1:4*i,1) = Dt;
end

% linear constraints
At = [-8     0     0     1     0
    0     -3     0     0     1
    0     0     -6     0     0
    0     0      0     -1     -3];
for i = 1: Num
    A(4*(i-1)+1:4*i,5*(i-1)+1:5*(i-1)+5) = At;
end

Bt = [diag(ones(1,4))];   % block matrix
for i = 1: Num
    B(4*(i-1)+1:4*i,4*(i-1)+1:4*i) = Bt;
end

bt = [0;8;3;0];
for i = 1: Num
    b(4*(i-1)+1:4*i,1) = bt;
end

Et =[1 1 0 0;
    0 0 1 1];
for i = 1: Num
    E(2*(i-1)+1:2*i,4*(i-1)+1:4*i) = Et;
end

ht = [8;5];
for i = 1: Num
    h(2*(i-1)+1:2*i,1) = ht;
end
F = zeros(1,5*Num);
r_le = 0;

Gt = [1 1 1 0 0];
for i = 1: Num
    G(i,5*(i-1)+1:5*i) = Gt;
end
r_ls = 2.*ones(Num,1);

% SCOP constraints
for i = 1 :  Num*4/2
    dd  = 1.0e-10.*ones(1,Num*4);
    dd(1,2*i-1:2*i) = 1;
    Qs{i} = diag(dd);
end
gst = [64;50];
gs = repmat(gst,Num*4/2,1);

ls = zeros(4*Num,Num*4);

% Number of blocks in matrix B
n_block = Num*4/2;                  
n_y = 2;
lb_y = zeros(1,n_y)';
ub_y = inf.*ones(1,n_y)';

% -------- Approach 1:use multi-cut generalized benders decomposition method to tackle a MISOCP-based model------- %
[OptX,OptY,OptValue,k,delta_T,delta_T_MP] = General_AMBD_socp(C,D,A,B,b,E,h,F,r_le,G,r_ls,Qs,gs,ls,n_block,n_y,lb_y,ub_y);

% time consumption of SP and MP
sum(max(delta_T)) + sum(delta_T_MP) - sum(sum(delta_T_MP)) 

% -------- Approach 2:use generalized benders decomposition method to tackle a MISOCP-based model------- %
[OptX,OptY,OptValue,k,delta_T,delta_T_MP] = General_MBD_socp(C,D,A,B,b,E,h,F,r_le,G,r_ls,Qs,gs,ls,n_block,n_y,lb_y,ub_y);

% time consumption of SP and MP
sum(max(delta_T)) + sum(delta_T_MP) - sum(sum(delta_T_MP)) 

% -------- Approach 3:use mosek tools to solve this problem------- %
t1 = clock;
% Specify the non-conic part of the problem.
n_x = size(C,1);n_y = size(D,1);
prob.c = [C' D' zeros(1,Num*4)];
prob.a = sparse([A B zeros(size(A,1),Num*4);
    zeros(size(E,1),n_x) E zeros(size(E,1),Num*4);
    F  zeros(size(F,1),n_y) zeros(size(F,1),Num*4);
    G   zeros(size(G,1),n_y) zeros(size(G,1),Num*4);
   zeros(Num*4,n_x+n_y) eye(Num*4)]);

prob.blc = [-inf.*ones(size(b,1),1);h;-inf;r_ls;gs];
prob.buc = [b ;h;r_le;r_ls;gs];

prob.blx = [zeros(n_x ,1)' zeros(n_y ,1)'   zeros(Num*4 ,1)'  ];
prob.bux = [ones(n_x,1)',  inf.*ones(1,n_y+Num*4)];

% Specify the number of cones.
prob.ints.sub= [1:n_x]';                  % x1~x5 are integer variables
prob.sol.int.xx = [zeros(1,n_x+n_y) gs']';
% Specify the number of cones.
prob.cones = cell(Num*2,1);
for j = 1 : Num
% The first cone is specified.
prob.cones{2*j-1}.type = 'MSK_CT_QUAD';     % MSK_CT_QUAD for Rotated Qudratic Cone
prob.cones{2*j-1}.sub = [n_x+n_y+2*j-1, n_x+4*j-3:n_x+4*j-2];
% The second cone is specified.
prob.cones{2*j}.type = 'MSK_CT_QUAD';     % MSK_CT_QUAD for Rotated Qudratic Cone
prob.cones{2*j}.sub = [n_x+n_y+2*j, n_x+4*j-1:n_x+4*j];
end
[r,res] = mosekopt('minimize echo(0)',prob);
% Display the primal solution.
OptY = res.sol.int.xx(n_x+1:n_x+n_y,1); % continuous solution
OptX = res.sol.int.xx(1:n_x,1);         % integer solution
OptValue = C'*OptX + D'*OptY ;
t2 = clock;
% running time is 0.1030 s
etime(t2,t1) 
% ------------------- END --------------------%


