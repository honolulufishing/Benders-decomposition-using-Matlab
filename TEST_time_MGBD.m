
clc;clear all
 
% Example 1:a socp problem
C = [7 7 7 7 7]';
D = [1 1 1 1]';
H = [1 1 1 1]';

% linear constraints
A = [-8     0     0     1     0
    0     -3     0     0     1
    0     0     -6     0     0
    0     0      0     -1     -3];
 
B = [diag(ones(1,4))];   % block matrix
b = [0;8;3;0];
 
E = [1 1 0 0;
    0 0 1 1];
h = [8;5];
 
F = zeros(1,5);
r_le = 0;
 
G = [1 1 1 0 0];
r_ls = 2;
 
% SCOP constraints
% yi'*Q*yi + l'*yi + g <=yj  (j>i)
% i.e.
% y1^2 - 9 < y2
% y3^2 - 25 < y4
Qij = [1 2;   % y1-y2
      3 4 ];  % y3-y4
Q_diag = [1; % y1'*1*y1
          1]; % y3*1*y3
l = [0 ;   % y1
    0 ;]; % y3
g = [-9;  % y1
    -26];  % y3

n_y = size(B,2);
lb_y = zeros(1,n_y)';
ub_y = inf.*ones(1,n_y)';

% Number of blocks in matrix B
n_block = 2;                  

% % -------- Approach 1:use multi-cut generalized benders decomposition method to tackle a MISOCP-based model------- %
% [OptX,OptY,OptValue,k,delta_T,delta_T_MP] = General_AMBD_socp(C,D,A,B,b,E,h,F,r_le,G,r_ls,Qs,gs,ls,n_block,n_y,lb_y,ub_y);
% 
% % time consumption of SP and MP
% sum(max(delta_T)) + sum(delta_T_MP) - sum(sum(delta_T_MP)) 
% 
% % -------- Approach 2:use modified generalized benders decomposition method to tackle a MISOCP-based model------- %
% [OptX,OptY,OptValue,k,delta_T,delta_T_MP] = General_MBD_socp(C,D,A,B,b,E,h,F,r_le,G,r_ls,Qs,gs,ls,n_block,n_y,lb_y,ub_y);
% 
% 
% % time consumption of SP and MP
% sum(max(delta_T)) + sum(delta_T_MP) - sum(sum(delta_T_MP)) 
% 
% -------- Approach 3:use generalized benders decomposition method to tackle a MISOCP-based model------- %
[OptX,OptY,OptValue,k,delta_T,delta_T_MP] = General_BD_socp(H,C,D,A,B,b,E,h,F,r_le,G,r_ls,Qij,Q_diag,g,l,lb_y,ub_y);


% time consumption of SP and MP
sum(delta_T,1) + sum(delta_T_MP,1)

% -------- Approach 4:use mosek tools to solve this problem------- %
t1 = clock;
% Specify the non-conic part of the problem.
n_x = size(C,1);n_y = size(D,1);
prob.c = [C' D' 0 0 0 0];
 
prob.qosubi = n_x + [1 2 3 4 ]';
prob.qosubj = n_x + [1 2 3 4 ]';
prob.qoval = [1.0 1 1.0 1.0]';

prob.a = sparse([A B zeros(size(A,1),4);
   zeros(size(E,1),n_x) E zeros(size(E,1),4);
   F  zeros(size(F,1),n_y) zeros(size(F,1),4);
 G   zeros(size(G,1),n_y) zeros(size(G,1),4);
   zeros(1,n_x) 0 1 0  0 -1 0 0 0;
   zeros(1,n_x) 0 0 0  1  0 -1 0 0;]);
prob.blc = [-inf.*ones(size(b,1),1);h;-inf;r_ls;g];
prob.buc = [b ;h;r_le;r_ls;g];  
prob.blx = [zeros(n_x ,1)' zeros(n_y ,1)'   zeros(2 ,1)' 0.5 0.5  ];
prob.bux = [ones(n_x,1)',  inf.*ones(1,n_y+2) 0.5 0.5];

% Specify the number of cones.
prob.ints.sub= [1:n_x]';                  % x1~x5 are integer variables
prob.sol.int.xx = [zeros(1,n_x+n_y+2) 0.5 0.5]';
% Specify the number of cones.
prob.cones = cell(2,1);
% The first cone is specified.
prob.cones{1}.type = 'MSK_CT_RQUAD';     % MSK_CT_QUAD for Rotated Qudratic Cone
prob.cones{1}.sub = [n_x+n_y+1, n_x+n_y+2+1,  n_x+1  ];
% The second cone is specified.
prob.cones{2}.type = 'MSK_CT_RQUAD';     % MSK_CT_QUAD for Rotated Qudratic Cone
prob.cones{2}.sub = [n_x+n_y+2, n_x+n_y+2+2,  n_x+3 ];
[r,res] = mosekopt('minimize echo(0)',prob);
% Display the primal solution.
OptY = res.sol.int.xx(n_x+1:n_x+n_y,1); % continuous solution
OptX = res.sol.int.xx(1:n_x,1);         % integer solution
OptValue = C'*OptX + D'*OptY ; 
t2 = clock;
% running time is 0.1030 s
etime(t2,t1) 
% ------------------- END --------------------%




