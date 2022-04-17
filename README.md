# Benders-decomposition-using-Matlab
This is a matlab-based package to solve a large-scale linear programming problem using classical benders decomposition and a socp problem using generalized benders decomposition method. 

For a linear problem, we use classical benders decomposition tools to solve this problem, which is cast as the following form:

         min  C*x+D*y
         s.t. A*x+B*y<=b; 
              E*y=h;
              F*x=r;
              x in {0,1},and y>=0
              
For a MISOCP problem, we use generalized benders decomposition method to solve this problem, which is cast as the following form:

         min  C*x+D*y
         s.t. A*x+B*y<=b; 
              E*y=h;
              F*x=r;
              y'*Q*y+l'*y<=g
              x in {0,1},and y>=0
where sub-problem is a SOCP-based model and realxed master model is a MILP-based model. Please kindly note that this MISOCP problem has two conditions to be satisfied: i) Convexity: the problem should be convex on y given the discrete variables x; ii) Linear separability: the problem should be linear on x given the continuous variables y.






If you have any questions, please feel free to contact me. Thank you.

Author: Chao Lei

Email: 21118924r@connect.polyu.hk 

April, 2022
