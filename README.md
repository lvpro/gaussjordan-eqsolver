# gaussjordan-eqsolver
Solves a system of up to 65k algebraic equations using Gauss-Jordan elimination.  

This was a component of a multi-threaded win32 application which taught linear algebra to kids.

Module Description:
- C++ object which solves integer-based linear algebraic equations based on Gauss-Jordan Elimination. The solver algorithm is a custom implementation slightly optimized and reordered compared to a more common implementation such as the one described in "Numerical Recipies."
- 100% integer arithmetic, no FPU needed
- Detects Infinite Solutions & No Solutions
- Detects Overflow
  
Requirements:
- The number of equations and unknowns submitted to the object must be equal.
- For purposes of 100% accuracy, input to the equation solver must be 100% integer-based; floating point calculation IS NOT SUPPORTED in this module.
- The maximum number of simulatenous equations (and thus unknowns) this module supports is 65535.
