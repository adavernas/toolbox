Toolbox for “A Solution Method for Continuous-Time General Equilibrium Models”

Adrien d’Avernas, Valentin Schubert and Quentin Vandeweyer

This repository contains the MATLAB toolbox to solve a variety of continuous-time model using the method from the paper "A Solution Method for Continuous-Time General Equilibrium Model," by Adrien d’Avernas, Valentin Schubert and Quentin Vandeweyer (2020). 

List of code files in the folder ‘files’:

- initpath.m: initialises paths required to solve a model
- mod_BruSan.m: specifies all model specific parameters, variables, and equations for a simple extension of the model provided by Brunnermeier and Sannikov (2014). A potential initial guess is also provided in this file.
- par_BruSan.m: this file collects all parameters required by the algorithm to solve the model, such as convergence criteria. The toolbox allows for different specifications for different models
- mod_DiTella.m: specifies all model specific parameters, variables, and equations for the model presented by DiTella (2016). A potential initial guess is also provided in this file.
- par_DiTella.m: collects all parameters required by the algorithm to solve the model based on Di Tella (2016), such as convergence criteria.

The folder ‘core’ contains all files needed by the algorithm. As a general rule, no file in this folder should be changed or moved. 

To solve a model using the algorithm, run the code

initpath

followed by 

solvemod(‘model name’, options)

in Matlab

Input options are:
- ‘write’: set to either ‘on’ or ‘off’. Default value is ‘off’. The option has to be set to ‘on’ whenever changes are made in the model equations.
- ‘guess’: set to either ‘on’ or ‘off’. Default value is ‘off’. If the algorithm should use a given matrix with an initial guess, the option has to be set to ‘on’.
- ‘dimensions’: takes input values ‘1D’ or 2D’. Default value is ‘2D’. If the model has only one state variable, for a faster convergence, the option can be set to ‘1D’.
- ‘method’: takes input variables ‘first_order’ and ‘forward’. Default is set for ‘first_order’. Defines which type of derivatives to be taken.
- ‘HJBupdate’: input values are ‘on’ and ‘off’. Default value is ‘off’
- ‘loop’: input values are ‘on’ and ‘off’. Default value is ‘off’
- 'parallel': input values are 'on' or 'off'. Lets you solve the model in parallel. 
- ‘search’: input values are ‘on’ and ‘off’. Default value is ‘off’
- ‘outerplot’: input values are ‘on’ and ‘off’. Default value is ‘off’. Plot figures during outer loop iteration
- ‘innerplot’: input values are ‘on’ and ‘off’. Default value is ‘off’. Plot figures during inner loop iteration
- ‘allplot’: input values are ‘on’ and ‘off’. Default value is ‘off’. Plot figures of all variables.
- ‘savegraph’: input values are ‘on’ and ‘off’. Default value is ‘off’. Save created figures
- ‘dispT’: takes integer as inputs. Default value is 10. Defines iterations to be displayed in the loop.


Example:

The toolbox provide two working examples. 
The first example is a simple extension of the Brunnermeier and Sannikov (2014) model. The derivations and set of equations required to numerically solve the model are provided in the file 'example_01_model_BruSan.pdf'. The model specific Matlab files are 'mod_BruSan.m' and 'par_BruSan.m'. The folders 'tmp_BruSan' and 'sol_BruSan' collect the function files created by Matlab to solve the model and the solution to the model respectively. If not provided, the folders are automatically created during the first iteration of the model. 
Once the model specific files are correct, the process to solve the model is relatively straightforward. 
For instance, when running the model for the first time, the following command is run:

initpath
solvemod('BruSan', 'write', 'on', 'loop', 'on')

The model 'BruSan' is solved. Matlab writes the necessary function files (saved in 'tmp_BruSan'), and when converging allows for extra looping around values. 
After having solved the model, the solution can be saved as initial guess for future use. If we have an initial guess:

initpath
solvemod('BruSan', 'write', 'on', 'guess', 'on')

The option 'write' has to be included whenever changes are made to the 'Model' part in the file 'mod_BruSan' as the function files have to be rewritten. When no changes in the 'Model' part are made, setting the option 'write' to 'off' is recommended as writing the function files takes time.  
Changes in the parameter values do not require rewriting the function files. Hence, the command can be written as:

initpath
solvemod('BruSan', 'guess', 'on', 'parallel', 'on)

Where the option 'parallel' allows to solve the model in parallel, which accelerates the process. 

The second example is based on the model by Di Tella (2016). The approach is analogue to the first example. 
