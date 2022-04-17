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




If you have any questions, please feel free to contact me. Thank you.

Author: Chao Lei

Email: 21118924r@connect.polyu.hk 

April, 2022
