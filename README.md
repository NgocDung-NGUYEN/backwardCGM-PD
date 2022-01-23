# *backwardCGM-PD*: **a backward elimination stepwise procedure** for **C**olored **G**raphical **M**odels for **P**aired **D**ata

R code repository for scripts implements the algorithm discussed in:
> [Doctoral thesis:] *Model selection for colored graphical models for paired data*

Colored graphs for paired data are represented by (&#x1D543;, E, &#x1D53C;<sub>L</sub>). In R language, a colored graph for paired data is displayed by a list of a vector L.as,  and the two-column matrices E and E.as, respectively.

## Contents

### Main scripts
Users can find the scripts that implement the backward elimination stepwise procedure for the twin lattice, developed in Section 3.3, and the similar approach for the model inclusion lattice, described in Section 3.4, which are, respectively,

- **backward_CGM_PD_tau.R**: implements the backward stepwise procedure for the lattice of RCON models for paired data with the twin order &#x227C;<sub>&tau;</sub>. The step-by-step idea is shown in Section 3.3 of the thesis.

- **backward_CGM_PD_submod.R**: implements the backward stepwise procedure for the lattice of RCON models for paired data equipped by the model inclusion order &#x227C;<sub>&#x1D46A;</sub> .  The step-by-step idea is shown in Section 3.4 of the thesis.

---

### Supplementary script
Users can find **supplementary_functions.R** that contains R functions supporting for functions in the main scripts. In particular,

- **intersectMat()**: returns the set intersection of two edge sets;
- **matdiff()**: returns the set difference of two edge sets;
- **tau()**, **tauMat()**: returns &tau;(A) and &tau;(E) where A is a set of vertices and E is the set of edges. respectively;
- **vLabel()**: returns the labels of two homologous sets of vertices L and R where L = {"L1",..., "Lq"}, R = {"R1",..., "Rq"}, and q = p/2;
- **fullEdges**: returns F<sub>V</sub>, F<sub>L</sub>, F<sub>R</sub>, and F<sub>T</sub> when p is given;
- **outEdges()**: returns E<sub>L</sub>, E<sub>R</sub>, E<sub>T</sub>, and E<sub>L</sub> &cap; &tau;(E<sub>R</sub>) when the edge set E is given;
- **pval.function()**: returns p-values of the likelihood ratio test of a model relative to the saturated model;
- some other functions that display the visualization from the list of inputs of a colored graph for paired data.


---


### Simulation script file
**simulation.R** contains the reproducible code for recorded results in **Table 3.2**. Multiple different scenarios in the simulation settings are described.

By running the script, the needed functions, the main scripts and supplementary script are sourced.

For reproduciblity, users can also find the a saved version of the output in .RData format of **simulation_results** file and simulated data in **simulated-data** if one does not want to spend time of running the code again.

---

### Application script file
**analysis-fMRI.R** contains commands for analyzed colored graphs selected from the backward elimination stepwise procedure on the twin lattice for real data with 36 variables. The outputs can be found in .RData format of **output-fMRI** file.