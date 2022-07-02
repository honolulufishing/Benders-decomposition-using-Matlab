# Multi-cut-Generalized-Benders-decomposition-using-Matlab
This is a matlab-based package to solve a large-scale linear or quadratic constrained convex programming problem using classical benders decomposition method, generalized and multi-cut benders decomposition methods. This matlab-based package is constructed with the MOSEK Aps 7.0 and CPLEX v12.7 using their MILP and SOCP solvers. Regarding the theory of multi-cut benders decomposition methods, please refer to our research work paper. Please kindly note that this matlab programming package is issued as our initial version, after which full version will be updated soon. 

For a linear problem, we use the classical benders decomposition method (CBD) to solve this problem, which is cast as the following form:

         min  C*x+D*y
         s.t. A*x+B*y<=b; 
              E*y=h;
              F*x<=r_le;
              G*x=r_ls;
              x in {0,1},and y>=0
              
For a MISOCP problem, we use the generalized benders decomposition method (GBD) to solve this problem with the following standard form:

         min  C*x+D*y
         s.t. A*x+B*y<=b; 
              E*y=h;
              F*x<=r_le;
              G*x=r_ls;
              y'*Q*y+l'*y<=g
              x in {0,1},and y>=0
              
where sub-problem is a SOCP-based model and realxed master model is a MILP-based model. Please kindly note that this MISOCP problem has two conditions to be satisfied: i) Convexity: the problem should be convex on y given the discrete variables x; ii) Linear separability: the problem should be linear on x given the continuous variables y.

Based on the above-mentioned MISOCP problem, we use the multi-cut generalized benders decomposition method (MGBD) to solve this problem expressed as the following standard form:

         min  C*x+D*y
         s.t. A*x+B*y<=b; 
              E*y=h;
              F*x<=r_le;
              G*x=r_ls;
              y'*Q*y+l'*y<=g
              x in {0,1},and y>=0

Please note that this model needs to satisfy the above-mentioned two conditions and also the third condition: iii) Linear independence: with given discrete variables x, different groups of continuous variables y are linearly independent.

To examine the efficiency of this family of benders decompositon methods, we compare the running time for a large-scale MISOCP-based optimization problem. This MISOCP problem contains 5 binary variables and 4 continuous variables in the basic scenario. Increasing the number of scenarios M, the number of binary variables and continuous variables are risen to 5*M and 4*M. Please refers to TEST_Time_large.m in our package.

Solvers:
- MILP solver : CPLEX v12.7
- SOCP solver: MOSEK v7

Methods for Comparison:
1. GBD - Generalized benders decomposition method
2. MGBD - Multi-cut generalized benders decomposition method
3. GSOCP - Global solver provided by MOSEK


                              Computation Time (seconds)
        No. of Scenarios  Num of variables      GBD        MGBD       GSOCP         
              1                 9              0.060       0.005      0.010
             10                90              0.104       0.006      0.012
             50               450             11.351       0.006      0.025
            100               900           out of memory  0.006      0.049
            500              4500           out of memory  0.010      0.325
            600              5400           out of memory  0.011      0.434

In summary, we can see that our MGBD runs more quickly than GBD and global solver by MOSEK. 
 
If you have any questions, please feel free to contact me. Thank you.

Author: Chao Lei

Email: 21118924r@connect.polyu.hk 

July, 2022
