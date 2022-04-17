# Multi-cut-Generalized-Benders-decomposition-using-Matlab
This is a matlab-based package to solve a large-scale linear or quadratic constrained convex programming problem using the classical benders decomposition, generalized benders decomposition method and multi-cut generalized benders decomposition method. 

For a linear problem, we use the classical benders decomposition tools to solve this problem, which is cast as the following form:

         min  C*x+D*y
         s.t. A*x+B*y<=b; 
              E*y=h;
              F*x<=r_le;
              G*x=r_ls;
              x in {0,1},and y>=0
              
For a MISOCP problem, we use the generalized benders decomposition method to solve this problem, which is expressed as the following form:

         min  C*x+D*y
         s.t. A*x+B*y<=b; 
              E*y=h;
              F*x<=r_le;
              G*x=r_ls;
              y'*Q*y+l'*y<=g
              x in {0,1},and y>=0
              
where sub-problem is a SOCP-based model and realxed master model is a MILP-based model. Please kindly note that this MISOCP problem has two conditions to be satisfied: i) Convexity: the problem should be convex on y given the discrete variables x; ii) Linear separability: the problem should be linear on x given the continuous variables y.

For a MISOCP problem, we use the generalized benders decomposition method to solve this problem, which is 

         min  C*x+D*y
         s.t. A*x+B*y<=b; 
              E*y=h;
              F*x<=r_le;
              G*x=r_ls;
              y'*Q*y+l'*y<=g
              x in {0,1},and y>=0

Please note that this model aslo satisfy the third conditon: iii) Linear independence: with given discrete variables x, groups of continuous variables y are linearly independent. These three characteristics indicate that a SOCP-based model can be decompsoed into a relaxed master problem with respect to x and many sub problems with respect to y, which sub-problem can be solved in parallel. 

For example, 

![CodeCogsEqn](https://user-images.githubusercontent.com/102128721/163716051-8febedae-9bf9-45fe-a64f-c1665ae1f087.png)

We use our the family of benders decompostion methods to solve this problem as

clc;clear all
 
% -------- Approach 1:use multi-cut generalized benders decomposition method to tackle a MISOCP-based model------- %
% where sub-problem is a SOCP-based model and relaxed master model is a MILP-based model
% w.r.t. the following form:
%         min  C*x+D*y
%         s.t. A*x+B*y<=b; x in {0,1},and y>=0
%              E*y=h;
%              F*x<=r_le;
%              G*x=r_ls;
%              y'*Q*y+l'*y<=g
 
C = [7 7 7 7 7]';
D = [1 1 1 1]';
 
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
Qs = {diag([1;1;1.0e-10;1.0e-10]);diag([1.0e-10;1.0e-10;1;1])}; % y1-y2; y3-y4
gs = [64;50];
ls = zeros(4,2);

 % -------- Approach 2:use multi-cut generalized benders decomposition method to tackle a MISOCP-based model------- %
% use multi-cut generalized benders decomposition method
n_block = 2;                  % Number of blocks in matrix B
n_y = n_block;

% t1 = clock;
[OptX,OptY,OptValue,k,delta_T,delta_T_MP] = General_MBD_socp(C,D,A,B,b,E,h,F,r_le,G,r_ls,Qs,gs,ls,n_block,n_y);
% t2 = clock;
% % running time is 0.0980 s
% etime(t2,t1) + sum(max(delta_T_MP)) - sum(sum(delta_T_MP))  

% time consumption of SP and MP
sum(sum(delta_T)) + sum(delta_T_MP) - sum(sum(delta_T_MP)) 


% -------- Approach 2:use generalized benders decomposition method to tackle a MISOCP-based model------- %
% use generalized benders decomposition method
%t1 = clock;
[OptX,OptY,OptValue,k,delta_T,delta_T_MP] = General_BD_socp(C,D,A,B,b,E,h,F,r_le,G,r_ls,Qs,gs,ls);
%t2 = clock;
% running time is 0.1030 s
% etime(t2,t1) 

% time consumption of SP and MP
sum(delta_T,1) + sum(delta_T_MP,1)


% -------- Approach 3:use mosek tools to solve this problem------- %
% Specify the non-conic part of the problem.
n_x = size(C,1);n_y = size(D,1);
prob.c = [C' D' 0 0];
prob.a = sparse([A B zeros(size(A,1),2);
   zeros(size(E,1),n_x) E zeros(size(E,1),2);
   F  zeros(size(F,1),n_y) zeros(size(F,1),2);
 G   zeros(size(G,1),n_y) zeros(size(G,1),2);
   zeros(1,n_x+n_y) 1 0;
   zeros(1,n_x+n_y) 0 1;]);
prob.blc = [-inf.*ones(size(b,1),1);h;-inf;r_ls;gs];
prob.buc = [b ;h;r_le;r_ls;gs];  
prob.blx = [zeros(n_x ,1)' zeros(n_y ,1)'   zeros(2 ,1)'  ];
prob.bux = [ones(n_x,1)',  inf.*ones(1,n_y+2)];
 
% Specify the number of cones.
prob.ints.sub= [1:n_x]';                  % x1~x5 are integer variables
prob.sol.int.xx = [zeros(1,n_x+n_y) gs']';
% Specify the number of cones.
prob.cones = cell(2,1);
% The first cone is specified.
prob.cones{1}.type = 'MSK_CT_QUAD';     % MSK_CT_QUAD for Rotated Qudratic Cone
prob.cones{1}.sub = [n_x+n_y+1, n_x+1:n_x+2];
% The second cone is specified.
prob.cones{2}.type = 'MSK_CT_QUAD';     % MSK_CT_QUAD for Rotated Qudratic Cone
prob.cones{2}.sub = [n_x+n_y+2, n_x+3:n_x+4];
 
t1 = clock;
[r,res] = mosekopt('minimize echo(0)',prob);
t2 = clock;
% Display the primal solution.
OptY = res.sol.int.xx(n_x+1:n_x+n_y,1); % continuous solution
OptX = res.sol.int.xx(1:n_x,1);         % integer solution
OptValue = C'*OptX + D'*OptY ; 

% running time is 0.1030 s
etime(t2,t1) 
% ------------------- END --------------------%






If you have any questions, please feel free to contact me. Thank you.

Author: Chao Lei

Email: 21118924r@connect.polyu.hk 

April, 2022
