# BCE_PDE
A time-dependent sediment diagenesis model for blue carbon ecosystems.
The main code is named as BCE_1 and the rest of the files are the functions that needs to be saved in the same folder as the main code.
The parameter "IC_PDEs" defines the type of initial conditions to be used to solve the partial different equations in the model. 
IC_PDEs equal to zero means no blue carbon added to the system and when is 1 it would start with the steady-state solution when 
there is no blue carbon (IC_PDEs = 0). For IC_PDEs = 1, the model needs to be run at IC_PDEs=0 first so that the steady-state solution of PDEs can be
used as intial conditions for the case when there is blue carbon in the model.

