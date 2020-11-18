# toolbox
Toolbox for "A Solution Method for Continuous-Time General Equilibrium Models"
Adrien d’Avernas, Valentin Schubert, and Quentin Vandeweyer (November 17, 2020)

List of code files in the folder ’files’
- main.m: solves the model and produces the figures
- model.m: includes all the model specific specifications. Model parameters and variables are declared in this file. Initial guesses are provided. Further, the equilibrium conditions for secondary variables as function of the endogenous variables are written first. The equations for endogenous variables are defined after. The different sections in the file are called to solve the model
- writefun.m: Matlab file setting up the model to be solved. Equations declared in model.m are called here.
- HJB.m: the model specific Hamiltonian-Jacobi Bellman equation is written separately from the other model equation in this file. The expression of the HJB after substituting for the optimal conditions has to be used.
- plotgraphs1D.m: creates the figures of a model with one state variable.
- plotgraphs2D.m: creates the figures of a model with two state variables.

Process to follow for solving a new model:

A short step-by-step description to use the algorithm and code for your own model.
1. Define variables and parameters in model.m
2. Rewrite the system of equation in the file model.m (section ’Model’ in the code)
3. The Hamiltonian-Jacobi Bellman Equations solved for the optimum solution has to be adapted in the file HJB.m. Important to also change the variables in the function file to the variables in your model.
4. After inputing the new model the set of parameters, endogenous variables, and secondary variables have to be changed in the main.m file.
5. When running the main.m file after having made changes in the 'Model' section of model.m file, it is imperative that the option ’par.write’ is set to ’on’ in the main.m file.
6. When running the main.m for the first time, you will probably not have an initial guess.mat, if this is the case, you have to set par.guess’ to ’off ’ and par.loop’ to ’on’. Further, you should provide sensible initial guesses for the endogenous variables (and secondary variables if available) in the section ’Initial Guess for X ’ in main.m
7. After having solved the model, the values of the equilibrium variables are saved as guess.mat and can be used by setting par.guess’ to ’on’.
8. If the model should be solved for different sets of parameters (without changing writefun.m), the option ’par.write’ can be set to ’off ’.
