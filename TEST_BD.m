
clc;clear all

% Example 1:a linear problem
% -------- Approach 1:use classical benders decomposition method to solve this problem------- %
% This kind of problems is w.r.t. the following form:
%         min  C*x+D*y
%         s.t. A*x+B*y<=b; x in {0,1},and y>=0
%              E*y=h;
%              F*x=r;

C = [7 7 7 7 7]';
D = [1 1 1 1 1]';

% linear constraints
A = [-diag([8;3;5;5;3])];
B = [diag(ones(1,5))];
b = [zeros(5,1)];

E = [1 0 0 1 1;
    0 1 0 0  1;
    0 0 1 1 0];
h = [8;3;5];

F = zeros(1,5);
r = 0;

% use classical benders decomposition method
[OptX,OptY,OptValue,k] = Classic_BD(C,D,A,B,b,E,h,F,r);
% ------------------- END --------------------%

 
% -------- Approach 2:use CPLEX tools to solve this problem------- %
A_t = [A B ];
vlb = zeros(1,size(C,1)+size(D,1));
vub= [ ones(1,size(D,1)) inf*ones(1,size(C,1)),];
ctype = char([ repmat({'I'},1,size(C,1)) repmat({'C'},1,size(D,1))])';
% use cplex toolbox      
[OptX,minZ,ExitflagBint]=cplexmilp([C' D'], A_t ,b ,[ F zeros(size(F,1),size(D,1));zeros(size(E,1),size(C,1)) E ], [r;h],[ ], [ ], [ ], vlb, vub, ctype, [ ] );
% ------------------- END --------------------%   


%% Example 2:a socp problem
% -------- Approach 1:use generalized benders decomposition method to tackle a MISOCP-based model------- %
% where sub-problem is a SOCP-based model and relaxed master model is a MILP-based model
% w.r.t. the following form:
%         min  C*x+D*y
%         s.t. A*x+B*y<=b; x in {0,1},and y>=0
%              E*y=h;
%              F*x=r;
%              y'*Q*y+l'*y<=g

C = [7 7 7 7 7]';
D = [1 1 1 1 1]';

% linear constraints
A = [-diag([8;3;5;5;3])];
B = [diag(ones(1,5))];
b = [zeros(5,1)];

E = [1 0 0 1 1;
    0 1 0 0  1;
    0 0 1 1 0];
h = [8;3;5];

F = zeros(1,5);
r = 0;

% SCOP constraints
Q = {eye(5,5)};
g = [40];
l = zeros(5,1);

[OptX,OptY,OptValue,k] = General_BD_socp(C,D,A,B,b,E,h,F,r,Q,g,l);
% ------------------- END --------------------%

% -------- Approach 2:use mosek tools to solve this problem------- %
% Specify the non-conic part of the problem.
n_x = size(C,1);n_y = size(D,1);
prob.c = [C' D' 0];
prob.a = sparse([A B zeros(size(A,1),1);
   zeros(size(E,1),n_x) E zeros(size(E,1),1);
   zeros(1,n_x+n_y) 1]);
prob.blc = [-inf.*ones(size(b,1),1);h;g];
prob.buc = [b ;h;g];  
prob.blx = [zeros(n_x ,1)' zeros(n_y ,1)'   zeros(1 ,1)'  ];
prob.bux = [ones(n_x,1)',  inf.*ones(1,n_y+1)];
% Specify the number of cones.
prob.ints.sub= [1:5]';                  % x1~x5 are integer variables
prob.sol.int.xx = [0 0 0 0 0 0 0 0 0 0 34]';
% Specify the number of cones.
prob.cones = cell(1,1);
% The first cone is specified.
prob.cones{1}.type = 'MSK_CT_QUAD';     % MSK_CT_QUAD for Rotated Qudratic Cone
prob.cones{1}.sub = [n_x+n_y+1, n_x+1:n_x+n_y];
[r,res] = mosekopt('minimize echo(0)',prob);
% Display the primal solution.
OptY = res.sol.int.xx(n_x+1:n_x+n_y,1); % continuous solution
OptX = res.sol.int.xx(1:n_x,1);         % integer solution
OptValue = C'*OptX + D'*OptY ; 
% ------------------- END --------------------%

