# *backwardCGM-PD*: **backward** stepwise procedure for **C**olored **G**raphical **M**odels for **P**aired **D**ata

R code repository for scripts implementing the algorithm discussed in:
> [Working paper:] *Model selection for colored graphical models for paired data*

h<sub>&theta;</sub>(x) = &theta;<sub>o</sub>x + &theta;<sub>1</sub>x&#x1D543;

(&#x1D543;, E, &#x1D53C;<sub>L</sub>)

Remain that, given the number of vertices p, the colored graphs representing models for paired data are defined by (&#x1D543;, E, &#x1D53C;<sub>L</sub>). In practice, a such colored graph is translated by a list of a vector L.as, the two-column matrices E and E.as, respectively.

## Contents

### Main scripts
Users can find the scripts that implement the procedure described in the section 4.1 and the similar approach applying the submodel relation.
More specifically:

- **backward_CGM_PD_submod.R**: implements the backward stepwise procedure for the class of RCON models for paired data adapting the partial ordering &#x227C;<sub>&tau;</sub>. The step-by-step idea is shown in the section 4.1.

- **backward_CGM_PD_submod.R**: implements the backward stepwise procedure for the class of RCON models for paired data equipped by the submodel relation which is induced by the ordering &#x227C;<sub>&#x1D46A;</sub> .

---

### Supplementary script
Users can find **supplementary_functions.R** that contains R functions supporting for functions in the main scripts. In particular,

- **intersectMat()**: returns the intersection of two edge sets;
- **matdiff()**: returns the difference of two edge sets;
- **tau()**, **tauMat()**: returns &tau;(A) and &tau;(E) where A is a set of vertices and E is the set of edges. respectively;
- **vLabel()**: returns the labels of two homologous sets of vertices L and R where L = {"L1",..., "Lq"}, R = {"R1",..., "Rq"}, and q = p/2;
- **fullEdges**: returns F<sub>V</sub>, F<sub>L</sub>, F<sub>R</sub>, and F<sub>T</sub> when p is given;
- **outEdges()**: returns E<sub>L</sub>, E<sub>R</sub>, E<sub>T</sub>, and E<sub>L</sub> &cap; &tau;(E<sub>R</sub>) when the edge set E is given;
- **pval.function()**: returns p-values of the likelihood ratio test of a model relative to the saturated model;
- some other functions that display the visualization from the input of colored graph.


---


### Simulation script
**simulation()** contains the reproducible code for results reported in **Table 1** of [1].

By running the script **fit_sgl_example.R**, the input simulated data (**input_example.RData**) is automatically loaded, and the scripts in the main folder sourced. The results can be visualized by typing in R the name of the object *summary_perf* in the console.

For reproduciblity, users can also find the a saved version of the output if one wants to avoid running the code again (**output_example.RData**).

---
