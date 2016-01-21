/*
	Author: Clint Kennedy 

	Module Description:
	- C++ object which solves integer-based linear algebraic equations
	based on Gauss-Jordan Elimination. The solver algorithm is a custom
	implementation slightly optimized and reordered compared to a more common
	implementation such as the one described in "Numerical Recipies."
	- 100% integer arithmetic, no FPU needed
	- Detects Infinite Solutions & No Solutions
	- Detects Overflow
	
	Requirements:
	- The number of equations and unknowns submitted to the object
	must be equal.
	- For purposes of 100% accuracy, input to the equation solver
	must be 100% integer-based; floating point calculation
	IS NOT SUPPORTED in this module. 
	- The maximum number of simulatenous equations (and thus unknowns)
	this module supports is 65535.
*/

/* 64-bit Integer Type Used In Overflow Checking */
/* Different Compilers Use Different Mechanisms Of Representation */

/* GCC / G++ */
#ifdef GCC_BUILD	
typedef unsigned long long int UINT64;
typedef long long int INT64;
#define UINT32MAX 42946972195LL
#define INT32MAX 2147483647LL
#define INT32MIN -2147483648LL

/* Visual C++ 6.0 */
#else
typedef signed __int64 INT64;
typedef unsigned __int64 UINT64;
#define UINT32MAX 42946972195I64
#define INT32MAX 2147483647I64
#define INT32MIN -2147483648I64
#endif

/* Definitions */
#define SOLVED 0x0001
#define NO_SOLUTIONS 0x0002
#define INFINITE_SOLUTIONS 0x0003
#define MEMORY_ERROR 0x0004
#define OVERFLOW 0x0005

/* The following datatype "fraction" is defined to offer an alternative
to floating-point datatypes. */
struct fraction
{
	unsigned int numerator;		/* Full 32-bits For 0-42946972195 */
	unsigned int denominator;
	unsigned int sign;	/* Holds Sign, Allows numerator & denominator Another Bit */
						/* 1 = Negative 0 = Positive */
};

/* eqsolver Class Defintion */
class eqsolver
{	
	/* Private Data */
	
	struct fraction **coefficient;	/* Holds N+1 x N Matrix Values On Which We May Operate */
	struct fraction **originalCoefficient;	/* Unaltered Storage For Matrix Values */
	unsigned short int eqCount;	/* # Of Simultaneous Equations In System */ 

	/* Private Methods */

	void swapRows(unsigned short int row1, unsigned short int row2, struct fraction **coeffPtr);	/* Swaps Matrix Rows */
	struct fraction reduce(struct fraction unreducedFraction);	/* Reduces Fractions */
	struct fraction divide(struct fraction dividend, struct fraction divisor);	/* Divides Fractions */
	struct fraction multiply(struct fraction fraction1, struct fraction fraction2);	/* Multiplies Fractions */
	struct fraction add(struct fraction fraction1, struct fraction fraction2);	/* Adds Fractions */
	void multiplyMatrixRow(unsigned short int row, struct fraction multiplier, struct fraction **coeffPtr);	/* Multiply Specified Matrix Row By Value "multiplier" */
	void divideMatrixRow(unsigned short int row, struct fraction divisor, struct fraction **coeffPtr);	/* Divide Specified Row By Value "divisor" */ 
	void addMatrixRows(unsigned short int row, unsigned short int rowToAdd, struct fraction **coeffPtr);	/* Add "rowToAdd" to "row" in specified matrix */

public:

	/* Public Data */

	int overFlow;	/* Overflow Flag To Be Used After Arithmetic Calculations 1 = Overflow 0 = No Overflow */
	struct fraction *solutionCoefficient;	/* Holds Solution Coefficients After "solve()" Is Called */
	
	/* Public Methods */

	/* Constructor */
	eqsolver()
	{
		coefficient = NULL;
		originalCoefficient = NULL;
		solutionCoefficient = NULL;
		eqCount = 0;
		overFlow = 0;
	}

	unsigned int setSystemEqCount(unsigned short int count);	/* Sets Dimensions */
	void setCoefficient(unsigned short int row, unsigned short int column, short int value);	/* Sets Coefficient Value */
	void setCoefficientFraction(unsigned short int row, unsigned short int column, short int numerator, short int denominator);	/* Sets Coefficient Value In Fraction Form */
	int getOriginalMatrixCoefficient(unsigned short int row, unsigned short int column);	/* Retrives Unaltered Matrix Coefficient */
	void getAlteredMatrixCoefficient(unsigned short int row, unsigned short int column, struct fraction &coefficientValue);	/* Retrieves Altered Matrix Coefficient */
	void getOriginalMatrixCoefficientFraction(unsigned int row, unsigned short int column, int *numerator, int *denominator);
	void swapRows(unsigned short int row1, unsigned short int row2);	/* Swaps Altered Matrix Rows */
	void multiplyMatrixRow(unsigned short int row, struct fraction multiplier);	/* Multiply Altered Matrix Row By Value "multiplier" */
	void divideMatrixRow(unsigned short int row, struct fraction divisor);	/* Divide Specified Row By Value "divisor" */ 
	void addMatrixRows(unsigned short int row, unsigned short int rowToAdd);	/* Add "rowToAdd" To "row" In Altered Matrix */
	unsigned int solveSystem(void);	/* Solves System Specified In originalCoefficient, Places Solution In solutionCoefficient Array */
	void cleanup(void);	/* Deallocates Memory */
};
