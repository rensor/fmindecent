# fmindescent
Matlab based optimizer framework based on the steepest descent method. The algorithm can handle linear equality, linear in-equality, non-linear in-equality, and non-linear equality constraints. The constraints are added to the objective function in an external penalization manner.  Here,  two different methods can be applied. The first method uses a merit function approach, where a global penalty parameters is used to make infeasible constraints expensive for the objective function. Alternatively, the Augmented Lagrange method is also available. This method updates the penalty parameters/Lagrange multipliers in each iteration, see [1] for details. The final Lagrange multipliers are available for the user to analyze. 
The descent direction can be determined using four different methods. 
•	Conjugate gradient method
•	Broyden–Fletcher–Goldfarb–Shanno (BFGS)
•	Davidon–Fletcher–Powell (DFP)
•	Method of Feasible Directions (MFD)
The current implementation of the MFD algorithm relies on the glpk linear solver used in Octave. In future releases, the user should be able to choose among other well-known linear solvers. 
In each iteration, the step size is determined by use of the golden section method.

[1]  Nocedal J, Wright SJ: Numerical Optimization, second edition, ISBN-10:0-387-303003-0, p 514.
 
Usage: see also examples.m 

myProblem = fmindecent(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon); % Default options are applied

myProblem = fmindecent(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options); % Input options structure

myProblem = fmindecent(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon ,'MaxIterations',100); % Input/modify individual options

[x,fval,exitflag,output] = myProblem.solve(); % Solve problem

Available options:

Initialize values to default

options = fmindecent.options();

Algorithm: 'CG'
Display: 'off'
ConstraintMethod: ‘AL’
LineSearch: ‘Golden’
InfeasibilityPenalization: 1000
MaxFunctionEvaluations: 1000
MaxIterations: 1000
OptimalityTolerance: 1.0000e-06
StepTolerance: 1.0000e-10
