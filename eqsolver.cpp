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

#include <stdlib.h>
#include <memory.h>
#include "eqsolver.h"

/*	The purpose of this function is to allocate and zero-initialize
	the memory required to store the N+1xN matrix coefficients. It 
	allocates space for an "original" matrix which will be unaltered,
	as well as a working copy. The original matrix is located by the
	member variable originalCoefficient. The altered matrix is located
	by the member variable coefficient.

	Parameters: 
		count - Unsigned 16-bit value between 0-65536 which represents
				the number of simultaneous equations to solve.
	Returns:
		1 on success
		0 in case of an error (such as memory allocation errors)
*/
unsigned int eqsolver::setSystemEqCount(unsigned short int count)
{
	unsigned int i,j;	/* Loop Counter */
	
	/* No Need To Allocate Space, Empty System */
	if(count == 0)
		return 1;
	
	eqCount = count;	/* Store # Of Simulatenous Equations In System */
	
	/* Allocate Array Of Matrix Row Pointers */
	coefficient = (struct fraction **) malloc(count * sizeof(struct fraction *));
	originalCoefficient = (struct fraction **) malloc(count * sizeof(struct fraction *));

	/* Allocate Array Of Fractions For Solution Storage*/
	solutionCoefficient = (struct fraction *) malloc(count * sizeof(struct fraction));

	/* Check For Memory Allocation Error */
	if((coefficient == NULL) || (originalCoefficient == NULL) || (solutionCoefficient == NULL))
		return 0;	/* System Most Likely Crashed At This Point */
					/* Clean-up Is Not Necessary, Reboot Is In Order */
	
	/* Zero Initialize Solution Array */
	for(i=0; i<count; i++)
	{
		solutionCoefficient[i].numerator = solutionCoefficient[i].denominator = 0;
		solutionCoefficient[i].sign = 0;

	}
	/* Allocate & Zero Initialize Row Coefficients */
	for(i=0; i<count; i++)
	{
		coefficient[i] = (struct fraction *) malloc((count+1) * sizeof(struct fraction));
		originalCoefficient[i] = (struct fraction *) malloc((count+1) * sizeof(struct fraction));
		
		/* Check For Memory Allocation Error */
		if((coefficient[i] == NULL) || (originalCoefficient[i] == NULL))
			return 0;	/* Review Explanation Above */

		/* Zero-Init Memory */
		for(j=0; j<=count; j++)
		{
			coefficient[i][j].numerator = coefficient[i][j].denominator = 0;
			originalCoefficient[i][j].numerator = originalCoefficient[i][j].denominator = 0;
			coefficient[i][j].sign = originalCoefficient[i][j].sign = 0;
		}
	}

	return 1;	/* Signals No Error */
}

/*	The purpose of this function is to set the specified matrix
	coefficient (in row,column order) to the specified value. This
	value may be any number between -32768 to 32767 inclusive.

	Parameters: 
		row - matrix row # (starting at 1)
		column - matrix column # (starting at 1)
		value - value at which to set the specified coefficient

	Returns:
		In case of error (out of bound matrix positioning), the
		function simply returns without altering any memory.
*/
void eqsolver::setCoefficient(unsigned short int row, unsigned short int column, short int value)
{
	/* Verify Matrix Bounds */
	if((row > eqCount) || (column > (eqCount+1)) || (row < 1) || (column < 1))
		return;	/* Out Of Bounds, Simply Return */

	/* Set Specified Coefficient Numerator To Specified Value */
	if(value < 0)
	{
		coefficient[(row-1)][(column-1)].numerator = (unsigned int) abs((int)value);
		coefficient[(row-1)][(column-1)].sign = 1;
		originalCoefficient[(row-1)][(column-1)].numerator = (unsigned int) abs((int)value);
		originalCoefficient[(row-1)][(column-1)].sign = 1;	
	}
	else
	{
		coefficient[(row-1)][(column-1)].numerator = (unsigned int) value;	/* Cast To 32-bits */
		coefficient[(row-1)][(column-1)].sign = 0;
		originalCoefficient[(row-1)][(column-1)].numerator = (unsigned int) value;
		originalCoefficient[(row-1)][(column-1)].sign = 0;
	}

	/* Set Denominator */
	if(value == 0)	/* Value Will Be 0/0, not 0/1 */
	{
		coefficient[(row-1)][(column-1)].denominator = 0;
		originalCoefficient[(row-1)][(column-1)].denominator = 0;
	}
	else
	{
		coefficient[(row-1)][(column-1)].denominator = 1;
		originalCoefficient[(row-1)][(column-1)].denominator = 1;
	}
}

/*	The purpose of this function is to set the specified matrix
	coefficient (in row,column order) to the specified value. This
	value may be any number between -32768 to 32767 inclusive. In 
	this setCoefficient() variation, both the numerator & denominator
	can be set.

	Parameters: 
		row - matrix row # (starting at 1)
		column - matrix column # (starting at 1)
		numerator - value at which to set the specified coefficient numerator
		denominator - value at which to set the specified coefficient denominator

	Returns:
		In case of error (out of bound matrix positioning), the
		function simply returns without altering any memory.
*/
void eqsolver::setCoefficientFraction(unsigned short int row, unsigned short int column, short int numerator, short int denominator)	/* Sets Coefficient Value In Fraction Form */
{
		/* Verify Matrix Bounds */
	if((row > eqCount) || (column > (eqCount+1)) || (row < 1) || (column < 1))
		return;	/* Out Of Bounds, Simply Return */

	/* Set Specified Coefficient Numerator To Specified Value */
	coefficient[(row-1)][(column-1)].numerator = (unsigned int) abs((int)numerator);
	originalCoefficient[(row-1)][(column-1)].numerator = (unsigned int) abs((int)numerator);
	
	/* Set Denominator */
	if(denominator == 0)	/* Value Will Be 0/0 (Prevents Division By Zero */
	{
		coefficient[(row-1)][(column-1)].numerator = 0;
		originalCoefficient[(row-1)][(column-1)].numerator = 0;
		coefficient[(row-1)][(column-1)].denominator = 0;
		originalCoefficient[(row-1)][(column-1)].denominator = 0;
	}
	else
	{
		coefficient[(row-1)][(column-1)].denominator = (unsigned int) abs((int)denominator);
		originalCoefficient[(row-1)][(column-1)].denominator = (unsigned int) abs((int)denominator);
	}

	/* Set Sign */
	if(numerator < 0)
	{
		if(denominator < 0)
		{
			coefficient[(row-1)][(column-1)].sign = 0;
			originalCoefficient[(row-1)][(column-1)].sign = 0;
		}
		else
		{
			coefficient[(row-1)][(column-1)].sign = 1;
			originalCoefficient[(row-1)][(column-1)].sign = 1;
		}
	}
	else
	{
		if(denominator < 0)
		{
			coefficient[(row-1)][(column-1)].sign = 1;
			originalCoefficient[(row-1)][(column-1)].sign = 1;
		}
		else
		{
			coefficient[(row-1)][(column-1)].sign = 0;
			originalCoefficient[(row-1)][(column-1)].sign = 0;
		}
	}
}

/*	The purpose of this function is to retrieve the value of the
	specified unaltered matrix coefficient (returned in fraction form).

	Parameters: 
		row - matrix row # (starting at 1)
		column - matrix column # (starting at 1)
		coefficientValue - reference to an integer		

	Returns:
		In case of error (out of bound matrix positioning), the
		function simply returns 0. Otherwise the value is returned.
*/
int eqsolver::getOriginalMatrixCoefficient(unsigned short int row, unsigned short int column)
{
	int coefficientValue;
	
	/* Verify Matrix Bounds */
	if((row > eqCount) || (column > (eqCount+1)) || (row < 1) || (column < 1))
		return 0;	/* Out Of Bounds, Simply Return, Note coefficientValue remains unchanged */
	
	coefficientValue = originalCoefficient[(row-1)][(column-1)].numerator;
	
	if(originalCoefficient[(row-1)][(column-1)].sign == 1)
		coefficientValue = coefficientValue * -1;

	return coefficientValue;
}

/* Same As Above, Slightly Modified */
void eqsolver::getOriginalMatrixCoefficientFraction(unsigned int row, unsigned short int column,
												  int *numerator, int *denominator)
{	
	/* Verify Matrix Bounds */
	if((row > eqCount) || (column > (eqCount+1)) || (row < 1) || (column < 1))
	{
		*numerator = *denominator = 0;	
		return;	/* Out Of Bounds, Simply Return, Note coefficientValue remains unchanged */
	}

	*numerator = originalCoefficient[(row-1)][(column-1)].numerator;
	
	if(originalCoefficient[(row-1)][(column-1)].sign == 1)
		*numerator = *numerator * -1;

	*denominator = originalCoefficient[(row-1)][(column-1)].denominator;
}

/*	The purpose of this function is to retrieve the value of the
	specified altered matrix coefficient (returned in fraction form).

	Parameters: 
		row - matrix row # (starting at 1)
		column - matrix column # (starting at 1)
		coefficientValue - reference to fraction structure (stores result)		

	Returns:
		In case of error (out of bound matrix positioning), the
		function simply returns.
*/
void eqsolver::getAlteredMatrixCoefficient(unsigned short row, unsigned short column, struct fraction &coefficientValue)
{
	/* Verify Matrix Bounds */
	if((row > eqCount) || (column > (eqCount+1)) || (row < 1) || (column < 1))
		return;	/* Out Of Bounds, Simply Return */
	
	coefficientValue.numerator = coefficient[(row-1)][(column-1)].numerator;
	coefficientValue.denominator = coefficient[(row-1)][(column-1)].denominator;
	coefficientValue.sign = coefficient[(row-1)][(column-1)].sign;
}

/*	The purpose of this function is to swap the specified matrix
	rows (this is used during the Gauss-Jordan algorithm). Since
	the only memory locations altered are 2 pointers, this swap
	method is extremely efficient.

	Parameters: 
		row1, row2 - rows to swap
		coeffPtr - pointer to coefficient data structure in which you
					wish to swap rows

	Returns:
		In case of error (out of bound matrix positioning), the
		function simply returns without altering any memory.
*/
void eqsolver::swapRows(unsigned short int row1, unsigned short int row2, struct fraction **coeffPtr)
{
	struct fraction *temp;	/* Swap Address Holder */

	if((row1 < 1) || (row2 < 1) || (row1 > eqCount) || (row2 > eqCount))
		return;	/* Out Of Bounds */

	/* Swap Row Pointers */
	temp = coeffPtr[(row1-1)];
	coeffPtr[(row1-1)] = coeffPtr[(row2-1)];
	coeffPtr[(row2-1)] = temp;
}

/*	The purpose of this function is reduce the specified fraction
	to its lowest possible form.

	Parameters: 
		unreducedFraction - fraction to reduce

	Returns:
		The result of the reduction is returned.
*/
struct fraction eqsolver::reduce(struct fraction unreducedFraction)
{
	struct fraction result;
	unsigned int gcf, num1, num2, temp;	/* Greatest Common Factor Solving */

	/* First Check If the unreducedFraction = 0 */
	/* The 2nd half of this condition is not needed, but added as a 
		precautionary measure to ENSURE non-divide-by-zero. */
	if((unreducedFraction.numerator == 0) || (unreducedFraction.denominator == 0))
	{
		result.numerator = 0;
		result.denominator = 0;
		result.sign = 0;
		return result;
	}
	
	/* Next Check If The unreducedFraction = 1 */
	if(unreducedFraction.numerator == unreducedFraction.denominator)
	{
		result.numerator = 1;
		result.denominator = 1;
		result.sign = unreducedFraction.sign;
		return result;
	}
	
	num1 = unreducedFraction.numerator;
	num2 = unreducedFraction.denominator;
		
	/* Determine The Greatest Common Factor (Using Euclid's Algorithm) */
	while (num2)
	{
		temp = num1;	
		num1 = num2;
		num2 = temp % num2;
	}

	gcf = num1;

	result.numerator = unreducedFraction.numerator / gcf;
	result.denominator = unreducedFraction.denominator / gcf;
	result.sign = unreducedFraction.sign;

	return result;
}

/*	The purpose of this function is to divide two "fraction" values
	and return the result to the caller.

	Parameters: 
		dividend - fraction which will be divided
		divisor - fraction by which to divide the dividend

	Returns:
		The result of the calculation is returned. Note zero-division
		is checked but should still be handled by the caller. If a 
		zero-division occurs, the dividend is returned unaltered.
*/
struct fraction eqsolver::divide(struct fraction dividend, struct fraction divisor)
{
	UINT64 overflowCheck;
	struct fraction result;
	
	/* Accomplish Division By Reciprocal Multiplication */
	
	/* Calculate Numerator */
	overflowCheck = (UINT64)dividend.numerator * (UINT64)divisor.denominator;
	if(overflowCheck > UINT32MAX)
	{
		overFlow = 1;	/* Set Overflow Flag */
		result.numerator = 0;
		result.denominator = 0;
		result.sign = 0;
		return result;
	}

	result.numerator = dividend.numerator * divisor.denominator;

	/* Calculate Denominator */
	overflowCheck = (UINT64)dividend.denominator * (UINT64)divisor.numerator;
	if(overflowCheck > UINT32MAX)
	{
		overFlow = 1;	/* Set Overflow Flag */
		result.numerator = 0;
		result.denominator = 0;
		result.sign = 0;
		return result;
	}

	result.denominator = dividend.denominator * divisor.numerator;

	/* Prevent Divide By Zero (Should Not Happen If Called Properly) */
	if((result.numerator != 0) && (result.denominator == 0))
		return dividend;	/* Return The Dividend Unaltered */
	
	/* Set Sign */
	if(dividend.sign == 1)	/* Dividend Negative? */
	{
		if(divisor.sign == 1)	/* Divisor Also Negative? */
			result.sign = 0;	/* Negative/Negative = Positive */
		else
			result.sign = 1;	/* Negative/Positive = Negative */
	}
	else	/* Dividend Positive */
	{
		if(divisor.sign == 1)	/* Divisor Negative? */
			result.sign = 1;	/* Positive/Negative = Negative */
		else
			result.sign = 0;	/* Positive/Positive = Positive */
	}
	
	/* Reduce Result */
	result = reduce(result);

	return result;
}

/*	The purpose of this function is to multiply two "fraction" values
	and return the result to the caller.

	Parameters: 
		fraction1 - fraction to by multiplied
		fraction2 - fraction to be multiplied by

	Returns:
		The result of the calculation is returned. Note zero-division
		is checked but should still be handled by the caller. If a 
		zero-division occurs, fraction1 is returned unaltered.
*/
struct fraction eqsolver::multiply(struct fraction fraction1, struct fraction fraction2)
{
	UINT64 overflowCheck;
	struct fraction result;
	
	/* Calculate Numerator */
	overflowCheck = (UINT64)fraction1.numerator * (UINT64)fraction2.numerator;
	if(overflowCheck > UINT32MAX)
	{
		overFlow = 1;	/* Set Overflow Flag */
		result.numerator = 0;
		result.denominator = 0;
		result.sign = 0;
		return result;
	}

	result.numerator = fraction1.numerator * fraction2.numerator;

	/* Calculate Denominator */
	overflowCheck = (UINT64)fraction1.denominator * (UINT64)fraction2.denominator;
	if(overflowCheck > UINT32MAX)
	{
		overFlow = 1;	/* Set Overflow Flag */
		result.numerator = 0;
		result.denominator = 0;
		result.sign = 0;
		return result;
	}

	result.denominator = fraction1.denominator * fraction2.denominator;

	/* Prevent Divide By Zero (Should Not Happen If Called Properly) */
	if((result.numerator != 0) && (result.denominator == 0))
		return fraction1;	/* Return fraction1  Unaltered */

	/* Set Sign */
	if(fraction1.sign == 1)	/* fraction1 Negative? */
	{
		if(fraction2.sign == 1)	/* fraction2 Also Negative? */
			result.sign = 0;	/* Negative*Negative = Positive */
		else
			result.sign = 1;	/* Negative*Positive = Negative */
	}
	else	/* Dividend Positive */
	{
		if(fraction2.sign == 1)	/* fraction2 Negative? */
			result.sign = 1;	/* Positive*Negative = Negative */
		else
			result.sign = 0;	/* Positive*Positive = Positive */
	}

	/* Reduce Result */
	result = reduce(result);

	return result;
}

/*	The purpose of this function is to add two "fraction" values
	and return the result to the caller.

	Parameters: 
		fraction1, fraction 2 - fractions to be added

	Returns:
		The result of the calculation is returned. Note zero-division
		is checked but should still be handled by the caller. If a 
		zero-division occurs, fraction1 is returned unaltered.
*/
struct fraction eqsolver::add(struct fraction fraction1, struct fraction fraction2)
{
	UINT64 overflowCheck;
	INT64 overflowCheckWithSign, num1OverflowCheckWithSign, num2OverflowCheckWithSign;
	int num1, num2, resultNum;
	struct fraction result;
	
	/* Verify Niether Fraction Is 0 */
	if((fraction1.numerator == 0) && (fraction1.denominator == 0))
		return fraction2;

	if((fraction2.numerator == 0) && (fraction2.denominator == 0))
		return fraction1;

	/* Do Calculations Here */
	/* Set Numerator */

	/* Are Both Numerators Positive? */
	if((fraction1.sign == 0) && (fraction2.sign == 0))
	{
		/* Stay In True Unsigned 32-bit Representation */
		overflowCheck = ((UINT64)fraction1.numerator * (UINT64)fraction2.denominator) + ((UINT64)fraction2.numerator * (UINT64)fraction1.denominator);
		if(overflowCheck > UINT32MAX)
		{
			overFlow = 1;	/* Set Overflow Flag */
			result.numerator = 0;
			result.denominator = 0;
			result.sign = 0;
			return result;
		}
		result.numerator = (fraction1.numerator * fraction2.denominator) + (fraction2.numerator * fraction1.denominator);
	}
	else	/* One Of The Numbers Is Negative */
	{
		/* First Check For Overflow In Necessary Calculations */
		if(fraction1.sign == 1)
		{
			num1OverflowCheckWithSign = (INT64)fraction1.numerator * -1;
			if(num1OverflowCheckWithSign < INT32MIN)
			{
				
				overFlow = 1;	/* Set Overflow Flag */
				result.numerator = 0;
				result.denominator = 0;
				result.sign = 0;
				return result;
			}
		}
		else
		{
			num1OverflowCheckWithSign = (INT64)fraction1.numerator;
			if(num1OverflowCheckWithSign > INT32MAX)
			{
				
				overFlow = 1;	/* Set Overflow Flag */
				result.numerator = 0;
				result.denominator = 0;
				result.sign = 0;
				return result;
			}
		}
		if(fraction2.sign == 1)
		{
			num2OverflowCheckWithSign = (INT64)fraction2.numerator * -1;
			if(num2OverflowCheckWithSign < INT32MIN)
			{
				overFlow = 1;	/* Set Overflow Flag */
				result.numerator = 0;
				result.denominator = 0;
				result.sign = 0;
				return result;
			}
		}
		else
		{
			num2OverflowCheckWithSign = (INT64)fraction2.numerator;
			if(num1OverflowCheckWithSign > INT32MAX)
			{
				overFlow = 1;	/* Set Overflow Flag */
				result.numerator = 0;
				result.denominator = 0;
				result.sign = 0;
				return result;
			}
		}
		
		overflowCheckWithSign = num1OverflowCheckWithSign * fraction2.denominator;
		if((overflowCheckWithSign > INT32MAX) || (overflowCheckWithSign < INT32MIN))
		{
			overFlow = 1;	/* Set Overflow Flag */
			result.numerator = 0;
			result.denominator = 0;
			result.sign = 0;
			return result;
		}

		overflowCheckWithSign = num2OverflowCheckWithSign * fraction1.denominator;
		if((overflowCheckWithSign > INT32MAX) || (overflowCheckWithSign < INT32MIN))
		{
			overFlow = 1;	/* Set Overflow Flag */
			result.numerator = 0;
			result.denominator = 0;
			result.sign = 0;
			return result;
		}

		overflowCheckWithSign = (num1OverflowCheckWithSign * fraction2.denominator) + (num2OverflowCheckWithSign * fraction1.denominator);
		if((overflowCheckWithSign > INT32MAX) || (overflowCheckWithSign < INT32MIN))
		{
			overFlow = 1;	/* Set Overflow Flag */
			result.numerator = 0;
			result.denominator = 0;
			result.sign = 0;
			return result;
		}

		/* Drop Down To Signed 32-bit Representation, Lose 1 Bit Of Size */
		if(fraction1.sign == 1)
			num1 = fraction1.numerator * -1;
		else
			num1 = fraction1.numerator;
		if(fraction2.sign == 1)
			num2 = fraction2.numerator * -1;
		else
			num2 = fraction2.numerator;

		resultNum = (num1 * fraction2.denominator) + (num2 * fraction1.denominator);
		if(resultNum < 0)
			result.sign = 1;
		else
			result.sign = 0;

		result.numerator = abs(resultNum);
	}

	/* Set Denominator */
	overflowCheck = (UINT64)fraction1.denominator * (UINT64)fraction2.denominator;
	if(overflowCheck > UINT32MAX)
	{
		overFlow = 1;	/* Set Overflow Flag */
		result.numerator = 0;
		result.denominator = 0;
		result.sign = 0;
		return result;
	}

	result.denominator = fraction1.denominator * fraction2.denominator;

	/* Prevent Divide By Zero (Should Not Happen If Called Properly) */
	if((result.numerator != 0) && (result.denominator == 0))
		return fraction1;	/* Return fraction1  Unaltered */

	/* Reduce Result */
	result = reduce(result);

	return result;
}

/*	The purpose of this function is to multiply a row in the specified 
	matrix by the specified fraction value.

	Parameters: 
		row - matrix row to be multiplied
		multiplier - fraction value by which to multiply
		coeffPtr - matrix to operate on

	Returns:
		Returns nothing. If row is out of bounds the function returns 
		without altering any memory. 
*/
void eqsolver::multiplyMatrixRow(unsigned short int row, struct fraction multiplier, struct fraction **coeffPtr)
{
	unsigned short int column;
	
	/* Perform Bounds Checking */
	if((row < 1) || (row > eqCount))
		return;	/* Out Of Bounds, Return */

	/* Multiply Each Value In Specified Row By "multiplier" Fraction */
	for(column=0; column<(eqCount+1); column++)
	{
		coeffPtr[(row-1)][column] = multiply(coeffPtr[(row-1)][column], multiplier);
		if(overFlow) return;	/* Overflow Occurred, No Use To Continue */
	}
}

/*	The purpose of this function is to divide a row in the specified
	matrix by the specified fraction value.

	Parameters: 
		row - matrix row to be divided
		divisor - fraction value by which to divide
		coeffPtr - matrix to operate on

	Returns:
		Returns nothing. If row is out of bounds the function returns 
		without altering any memory. 
*/
void eqsolver::divideMatrixRow(unsigned short int row, struct fraction divisor, struct fraction **coeffPtr)
{
	unsigned short int column;
	
	/* Perform Bounds Checking */
	if((row < 1) || (row > eqCount))
		return;	/* Out Of Bounds, Return */

	/* Divide Each Value In Specified Row By "divisor" Fraction */
	for(column=0; column<(eqCount+1); column++)
	{
		coeffPtr[(row-1)][column] = divide(coeffPtr[(row-1)][column], divisor);
		if(overFlow) return;	/* Overflow Occurred, No Use To Continue */
	}
}

/*	The purpose of this function is to add a row to another row
	in the specified matrix.

	Parameters: 
		row - matrix row to be added to
		rowToAdd - matrix row to add
		coeffPtr - matrix to operate on

	Returns:
		Returns nothing. If row is out of bounds the function returns 
		without altering any memory. 
*/
void eqsolver::addMatrixRows(unsigned short int row, unsigned short int rowToAdd, struct fraction **coeffPtr)
{
	unsigned short int column;
	
	/* Perform Bounds Checking */
	if((row < 1) || (row > eqCount) || (rowToAdd < 1) || (rowToAdd > eqCount))
		return;	/* Out Of Bounds, Return */

	/* Add Each Value In Specified rowToAdd To row */
	for(column=0; column<(eqCount+1); column++)
	{
		coeffPtr[(row-1)][column] = add(coeffPtr[(row-1)][column], coeffPtr[(rowToAdd-1)][column]);
		if(overFlow) return;	/* Overflow Occurred, No Use To Continue */
	}
}

/*	The purpose of this function is to swap the altered matrix
	rows (this is used during the Gauss-Jordan algorithm). Since
	the only memory locations altered are 2 pointers, this swap
	method is extremely efficient. Note this is an overloaded version
	of the private method swapRows().

	Parameters: 
		row1, row2 - rows to swap

	Returns:
		In case of error (out of bound matrix positioning), the
		function simply returns without altering any memory.
*/
void eqsolver::swapRows(unsigned short int row1, unsigned short int row2)
{
	struct fraction *temp;	/* Swap Address Holder */

	if((row1 < 1) || (row2 < 1) || (row1 > eqCount) || (row2 > eqCount))
		return;	/* Out Of Bounds */

	/* Swap Row Pointers */
	temp = coefficient[(row1-1)];
	coefficient[(row1-1)] = coefficient[(row2-1)];
	coefficient[(row2-1)] = temp;
}

/*	The purpose of this function is to multiply a row in the altered
	matrix by the specified fraction value. Note this is an overloaded
	version of the private method multiplyMatrixRow().

	Parameters: 
		row - matrix row to be multiplied
		multiplier - fraction value by which to multiply

	Returns:
		Returns nothing. If row is out of bounds the function returns 
		without altering any memory. 
*/
void eqsolver::multiplyMatrixRow(unsigned short int row, struct fraction multiplier)
{
	unsigned short int column;
	
	/* Perform Bounds Checking */
	if((row < 1) || (row > eqCount))
		return;	/* Out Of Bounds, Return */

	/* Multiply Each Value In Specified Row By "multiplier" Fraction */
	for(column=0; column<(eqCount+1); column++)
	{
		coefficient[(row-1)][column] = multiply(coefficient[(row-1)][column], multiplier);
		if(overFlow) return;	/* Overflow Occurred, No Use To Continue */
	}
}

/*	The purpose of this function is to divide a row in the altered
	matrix by the specified fraction value. Note this is an overloaded
	version of the private method divideMatrixRow().

	Parameters: 
		row - matrix row to be divided
		divisor - fraction value by which to divide

	Returns:
		Returns nothing. If row is out of bounds the function returns 
		without altering any memory. 
*/
void eqsolver::divideMatrixRow(unsigned short int row, struct fraction divisor)
{
	unsigned short int column;
	
	/* Perform Bounds Checking */
	if((row < 1) || (row > eqCount))
		return;	/* Out Of Bounds, Return */

	/* Divide Each Value In Specified Row By "divisor" Fraction */
	for(column=0; column<(eqCount+1); column++)
	{
		coefficient[(row-1)][column] = divide(coefficient[(row-1)][column], divisor);
		if(overFlow) return;	/* Overflow Occurred, No Use To Continue */
	}
}

/*	The purpose of this function is to add a row to another row
	in the altered matrix. Note this is an overloaded version of 
	the private method addMatrixRows().

	Parameters: 
		row - matrix row to be added to
		rowToAdd - matrix row to add

	Returns:
		Returns nothing. If row is out of bounds the function returns 
		without altering any memory. 
*/
void eqsolver::addMatrixRows(unsigned short int row, unsigned short int rowToAdd)
{
	unsigned short int column;
	
	/* Perform Bounds Checking */
	if((row < 1) || (row > eqCount) || (rowToAdd < 1) || (rowToAdd > eqCount))
		return;	/* Out Of Bounds, Return */

	/* Add Each Value In Specified rowToAdd To row */
	for(column=0; column<(eqCount+1); column++)
	{
		coefficient[(row-1)][column] = add(coefficient[(row-1)][column], coefficient[(rowToAdd-1)][column]);
		if(overFlow) return;	/* Overflow Occurred, No Use To Continue */
	}
}

/*	The purpose of this function is to solve the system specified
	in the "original" matrix coefficients.

	Parameters: 
		None

	Returns:
		Returns SOLVED, NO_SOLUTIONS, or INFINITE_SOLUTIONS.
		If SOLVED, the solution coefficients can be located by
		member variable solutionCoefficient. MEMORY_ERROR is returned
		on allocation errors. Returns OVERFLOW on 32-bit Overflow.
*/
unsigned int eqsolver::solveSystem(void)
{
	unsigned short int i, j;
	unsigned short int row, column, nonZeroFound;
	short int rowCounter;
	struct fraction **coeffPtr;
	struct fraction multiplier;
	struct fraction solutionCheck;
	
	overFlow = 0;	/* Reset Overflow Flag */
	
	/* Create "Working Copy" Of Matrix To Solve */

	coeffPtr = (struct fraction **) malloc(eqCount * sizeof(struct fraction *));
	
	if(coeffPtr == NULL)
		return MEMORY_ERROR;

	/* Allocate & Zero Initialize Row Coefficients */
	for(i=0; i<eqCount; i++)
	{
		coeffPtr[i] = (struct fraction *) malloc((eqCount+1) * sizeof(struct fraction));
	
		/* Check For Memory Allocation Error */
		if((coeffPtr[i] == NULL))
			return MEMORY_ERROR;	/* Review Explanation Above */		
	}

	/* Copy Matrix Coefficients From "original" Matrix To Working Copy */
	for(i=0; i<eqCount; i++)
		for(j=0; j<(eqCount+1); j++)
			coeffPtr[i][j] = originalCoefficient[i][j];

	row = 0;
	column = 0;

	/* While Each Column Not Reduced */
	while(column < eqCount)
	{
		/* If [row, column] = 0, Perform Full Pivoting */
		if((coeffPtr[row][column].numerator == 0) &&
			(coeffPtr[row][column].denominator == 0))
		{
			rowCounter = row+1;	/* Check For Replacement Below Only */
			nonZeroFound = 0;	/* Set Flag */

			while(!nonZeroFound)	/* While Suitable Pivot Not Found */
			{
				while(rowCounter < eqCount)	/* Process All Rows Below */
				{
					/* If Suitable Pivot Found, Swap Rows */
					if((coeffPtr[rowCounter][column].numerator != 0) &&
						(coeffPtr[rowCounter][column].denominator != 0))
					{
						swapRows((row+1), (rowCounter+1), coeffPtr);
						nonZeroFound = 1;
						break;
					}
					else	/* Proceed To Check Next Possible Pivot */
						rowCounter = rowCounter+1;
				}

				if(!nonZeroFound)	/* Ultimately Will Be Infinite Or No Solutions */
				{
					rowCounter = row+1;
					column = column+1;	/* Skip Column, Goto Next */

					if(column == eqCount) /* No Pivots Available? */
					{
						if((coeffPtr[row][column].numerator == 0) &&
							(coeffPtr[row][column].denominator == 0))
								return INFINITE_SOLUTIONS;
						else	
						{
							/* First Check For A Row Of All Zeros, If One Exists, STILL INFINITE_SOLUTIONS */
							/* row, column, & nonZeroFound Won't Be Used Anymore */
							for(row=0; row<eqCount; row++)
							{
								nonZeroFound = 0;
								for(column=0; column<=eqCount; column++)
									if(coeffPtr[row][column].numerator != 0)
										nonZeroFound = 1;
								if(!nonZeroFound)
									return INFINITE_SOLUTIONS;
							}

							return NO_SOLUTIONS;
						}
					}
				}
			}
		}
		
		/* Is Pivot-Point = 1? If Not, Divide To Make It So */
		if(!((coeffPtr[row][column].numerator == 1) && (coeffPtr[row][column].denominator == 1) && (coeffPtr[row][column].sign == 0)))
		{	
			divideMatrixRow((row+1), coeffPtr[row][column], coeffPtr);
			if(overFlow) return OVERFLOW;	/* Overflow Occurred, No Reason To Continue */
		}

		/* Make Pivot -1 */
		for(i=0; i<=eqCount; i++)
		{	
			if((coeffPtr[(row)][i].numerator == 0) && (coeffPtr[(row)][i].denominator == 0))
				continue;
			
			if(coeffPtr[(row)][i].sign == 0)
				coeffPtr[(row)][i].sign = 1;
			else
				coeffPtr[(row)][i].sign = 0;
		}

		/* Clear Out Column Above Row */
		for(rowCounter=row-1; rowCounter>=0; rowCounter--)
		{
			if((coeffPtr[rowCounter][column].numerator == 0) &&
				(coeffPtr[rowCounter][column].denominator == 0))
					continue;	/* Column Already Clear */
			
			multiplier = coeffPtr[rowCounter][column];
			multiplyMatrixRow((row+1), multiplier, coeffPtr);
			if(overFlow) return OVERFLOW;	/* Overflow Occurred, No Reason To Continue */
			addMatrixRows((rowCounter+1), (row+1), coeffPtr);
			if(overFlow) return OVERFLOW;	/* Overflow Occurred, No Reason To Continue */
			divideMatrixRow((row+1), multiplier, coeffPtr);
			if(overFlow) return OVERFLOW;	/* Overflow Occurred, No Reason To Continue */
		}

		/* Clear Out Column Below Row */
		for(rowCounter=row+1; rowCounter<eqCount; rowCounter++)
		{
			if((coeffPtr[rowCounter][column].numerator == 0) &&
				(coeffPtr[rowCounter][column].denominator == 0))
					continue;	/* Column Already Clear */
			
			multiplier = coeffPtr[rowCounter][column];
			multiplyMatrixRow((row+1), multiplier, coeffPtr);
			if(overFlow) return OVERFLOW;	/* Overflow Occurred, No Reason To Continue */
			addMatrixRows((rowCounter+1), (row+1), coeffPtr);
			if(overFlow) return OVERFLOW;	/* Overflow Occurred, No Reason To Continue */
			divideMatrixRow((row+1), multiplier, coeffPtr);
			if(overFlow) return OVERFLOW;	/* Overflow Occurred, No Reason To Continue */
		}

		/* Make Pivot -1 */
		for(i=0; i<=eqCount; i++)
		{
			if((coeffPtr[(row)][i].numerator == 0) && (coeffPtr[(row)][i].denominator == 0))
				continue;
			if(coeffPtr[(row)][i].sign == 0)
				coeffPtr[(row)][i].sign = 1;
			else
				coeffPtr[(row)][i].sign = 0;
		}
		
		/* Proceed To Next Pivot Point */
		row++;	
		column++;
	}

	/* Perform Added "Checking" On Result */
	if((coeffPtr[(eqCount-1)][(eqCount-1)].numerator == 0) &&
		(coeffPtr[(eqCount-1)][(eqCount-1)].denominator == 0))
	{
		if((coeffPtr[(eqCount-1)][eqCount].numerator == 0) &&
			(coeffPtr[(eqCount-1)][eqCount].denominator == 0))
			return INFINITE_SOLUTIONS;
		else
			return NO_SOLUTIONS;
	}
	else
	{
		/* System IS In Reduced Echolon Form But The Solution
			May Be Incorrect, Thus It Needs To Be Checked 
			i = row, j = column */

		for(i=0; i<eqCount; i++)
		{
			solutionCheck.numerator = 0;
			solutionCheck.denominator = 0;
			
			/* Total Row */
			for(j=0; j<eqCount; j++)
			{
				solutionCheck = add(solutionCheck, multiply(originalCoefficient[i][j], coeffPtr[j][eqCount]));
				if(overFlow) return OVERFLOW;	/* Overflow Occurred, No Reason To Continue */
			}

			if((solutionCheck.numerator != originalCoefficient[i][eqCount].numerator) ||
				(solutionCheck.denominator != originalCoefficient[i][eqCount].denominator) ||
				(solutionCheck.sign != originalCoefficient[i][eqCount].sign))
				return NO_SOLUTIONS;
		}
	}

	for(i=0; i<eqCount; i++)
		solutionCoefficient[i] = coeffPtr[i][eqCount];

	return SOLVED;
}

/*	The purpose of this function is to deallocate and "cleanup"
	memory the equation solver has used, thus "reseting" it.

	Parameters: 
		None

	Returns:
		None
*/
void eqsolver::cleanup(void)
{
	unsigned short int i;
	
	/* Deallocate Storage For solutionCoefficient */
	if(solutionCoefficient != NULL)
		free(solutionCoefficient);

	/* Deallocate Storage For "coefficient" & "originalCoefficient"
		matrix storages */
	if(coefficient != NULL)
	{
		/* Delete Rows */
		for(i=0; i<eqCount; i++)
			free(coefficient[i]);
		/* Delete Row Pointers */
		free(coefficient);
	}
	if(originalCoefficient != NULL)
	{
		/* Delete Rows */
		for(i=0; i<eqCount; i++)
			free(originalCoefficient[i]);
		/* Delete Row Pointers */
		free(originalCoefficient);
	}

	/* Reset eqCount to Zero, Reset Pointers To NULL */
	eqCount = 0;
	solutionCoefficient = NULL;
	coefficient = NULL;
	originalCoefficient = NULL;
	overFlow = 0;

	/* Done, Return */
}