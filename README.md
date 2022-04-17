# Multi-cut-Generalized-Benders-decomposition-using-Matlab
This is a matlab-based package to solve a large-scale linear or quadratic constrained convex programming problem using the classical benders decomposition (CBD), generalized benders decomposition method (GBD) and multi-cut generalized benders decomposition method (MGBD). The theory of multi-cut generalized benders decomposition method can be found in our research paper. Please note that the matlab-based codes are constructed on the MOSEK ApS and CPLEX toolbox. The uploaded matlab files are our initial codes, but the complete version will be issued very soon. 

For a linear problem, we use the classical benders decomposition method to solve this problem, which is cast as the following form:

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

For a MISOCP problem, we use the multi-cut generalized benders decomposition method to solve this problem, which is 

         min  C*x+D*y
         s.t. A*x+B*y<=b; 
              E*y=h;
              F*x<=r_le;
              G*x=r_ls;
              y'*Q*y+l'*y<=g
              x in {0,1},and y>=0

Please note that this model aslo satisfy the third conditon: iii) Linear independence: with given discrete variables x, groups of continuous variables y are linearly independent. These three characteristics indicate that a SOCP-based model can be decompsoed into a relaxed master problem with respect to x and multiple sub problems with respect to y, which multiple sub-problems can be solved in parallel. 

For example, 

![CodeCogsEqn](https://user-images.githubusercontent.com/102128721/163716575-bdaa3d1a-771c-437f-a1a4-83eb652e0024.png)

We use our the family of benders decompostion methods to solve this problem as

     clc;clear all
     C = [7 7 7 7 7]';
     D = [1 1 1 1]';
 
     % linear constraints
     A = [-8     0     0     1     0
         0     -3     0     0     1
         0     0     -6     0     0
         0     0      0     -1     -3];
 
     B = [diag(ones(1,4))];   % two block matrices
     b = [0;8;3;0];
 
     E = [1 1 0 0; 0 0 1 1];
     h = [8;5];
 
     F = zeros(1,5);
     r_le = 0;
 
     G = [1 1 1 0 0];
     r_ls = 2;
 
     % SOCP constraints
     Qs = {diag([1;1;1.0e-10;1.0e-10]);diag([1.0e-10;1.0e-10;1;1])}; 
     gs = [64;50];
     ls = zeros(4,2);

     % -------- Approach 1:use multi-cut generalized benders decomposition method to tackle a MISOCP-based model------- %
     % use multi-cut generalized benders decomposition method
     n_block = 2;                  % Number of blocks in matrix B
     n_y = n_block;
     [OptX,OptY,OptValue,k,delta_T,delta_T_MP] = General_MBD_socp(C,D,A,B,b,E,h,F,r_le,G,r_ls,Qs,gs,ls,n_block,n_y);
     % time consumption of SP and MP
     sum(sum(delta_T)) + sum(delta_T_MP) - sum(sum(delta_T_MP)) 

     % -------- Approach 2:use generalized benders decomposition method to tackle a MISOCP-based model------- %
     % use generalized benders decomposition method
     [OptX,OptY,OptValue,k,delta_T,delta_T_MP] = General_BD_socp(C,D,A,B,b,E,h,F,r_le,G,r_ls,Qs,gs,ls);
     % time consumption of SP and MP
     sum(delta_T,1) + sum(delta_T_MP,1)
     % ------------------- END --------------------%

The running time for MGBD is 0.0120 s, while the time cost using GBD is 0.0360 s. The MGBD method runs more quickly than GBD method for this model only with two block matrices.

If you have any questions, please feel free to contact me. Thank you.

Author: Chao Lei

Email: 21118924r@connect.polyu.hk 

April, 2022
