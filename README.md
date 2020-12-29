Toolbox for “A Solution Method for Continuous-Time General Equilibrium Models”

Adrien d’Avernas, Valentin Schubert and Quentin Vandeweyer

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
- ‘search’: input values are ‘on’ and ‘off’. Default value is ‘off’
- ‘outerplot’: input values are ‘on’ and ‘off’. Default value is ‘off’. Plot figures during outer loop iteration
- ‘innerplot’: input values are ‘on’ and ‘off’. Default value is ‘off’. Plot figures during inner loop iteration
- ‘allplot’: input values are ‘on’ and ‘off’. Default value is ‘off’. Plot figures of all variables.
- ‘savegraph’: input values are ‘on’ and ‘off’. Default value is ‘off’. Save created figures
- ‘dispT’: takes integer as inputs. Default value is 10. Defines iterations to be displayed in the loop.
