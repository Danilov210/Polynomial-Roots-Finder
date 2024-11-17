import math
import time
import warnings
import numpy as np

warnings.filterwarnings("ignore")


def getSign(n):
    """
    Check the sign of a number.

    :param n: The number to check.
    :return: 1 if positive, -1 if negative, 0 if zero.
    """
    return 1 if n > 0 else 0 if n == 0 else -1


def root_finder_newton_raphson(coefficients, x, epsilon, maxDiv):
    """
    Use the Newton-Raphson method to find a root of a polynomial.

    :param coefficients: Coefficients of the polynomial
    :param x: Initial guess
    :param epsilon: Tolerance for stopping condition
    :param maxDiv: Maximum divergence allowed
    :return: Approximate root or None if not found
    """
    iterations = 4  # Number of iterations
    derivative_coefficients = calculate_poly_derivative(coefficients)  # Calculate derivative coefficients
    xn = x  # Initialize current guess to initial guess
    for i in range(iterations): #i = Integer , range=Object that represents an interval of integers
        if abs(xn) > 1:
            # Calculate f(x)/f'(x) using Horner's method for large xn
            div = calculate_polynomial_divided_by_derivative(coefficients, derivative_coefficients, xn)
        else:
            # Calculate f(x) and f'(x) for small xn
            func_value = calculate_polynomial(coefficients, xn)
            deriv_value = calculate_polynomial(derivative_coefficients, xn)
            if deriv_value == 0:
                div = None  # If the derivative is zero, Newton's method fails
            else:
                div = func_value / deriv_value  # Calculate f(x)/f'(x)
        
        if div is None or abs(div) > maxDiv:
            # If the derivative is zero or division is too large, Newton's method fails
            return None
        
        # Calculate the next guess using Newton's method
        xn = x - div
        if abs(div) <= epsilon:
            return xn  # Return root if within tolerance
        
        # Update current guess for next iteration
        x = xn
    
    # Return None if root not found within maximum iterations
    return None


def root_finder_bisection(coefficients, a, b, epsilon):
    """
    Use the Bisection method to find a root of a polynomial.
    
    :param coefficients: Coefficients of the polynomial
    :param a: Left bound
    :param b: Right bound
    :param epsilon: Tolerance for stopping condition
    :return: Approximate root
    """
    # Calculate polynomial value at left bound a
    a_value = calculate_polynomial(coefficients, a)
    
    # Iterate until the interval [a, b] is smaller than epsilon
    while b - a > epsilon:
        # Calculate midpoint of the interval
        mid = (b + a) / 2
        # Calculate polynomial value at midpoint
        mid_value = calculate_polynomial(coefficients, mid)
        
        # If the sign of the midpoint value is the same as the sign of a_value
        if getSign(a_value) * getSign(mid_value) > 0:
            # Update left bound to be the midpoint
            a = mid
            a_value = mid_value
        else:
            # Update right bound to be the midpoint
            b = mid
    
    # Return the midpoint as the approximate root
    return (a + b) / 2


def calculate_polynomial_divided_by_derivative(func_coefficients, deriv_coefficients, x):
    """
    Evaluate f(x)/f'(x) using Horner's method to prevent overflow.

    :param func_coefficients: Coefficients of the polynomial f(x)
    :param deriv_coefficients: Coefficients of the derivative f'(x)
    :param x: Point at which to evaluate
    :return: f(x)/f'(x) or None if the derivative is zero
    """
    x_inverted = 1/x
    func_value = np.longdouble(0)
    deriv_value = np.longdouble(0)
    
    # Evaluate polynomial using Horner's method with inverted x
    for coefficient in reversed(func_coefficients):
        func_value = func_value * x_inverted + coefficient
    
    # Evaluate derivative using Horner's method with inverted x
    for coefficient in reversed(deriv_coefficients):
        deriv_value = deriv_value * x_inverted + coefficient
    
    # Calculate function value at x
    func_value = func_value * x
    
    # Check for division by zero
    if deriv_value == 0:
        return None
    
    return func_value / deriv_value


def calculate_poly_derivative(coefficients):
    """
    Calculate the derivative of a polynomial.

    :param coefficients: The coefficients of the original polynomial.
    :return: The coefficients of the derivative.
    """
    # Length of the coefficients array
    n = len(coefficients)
    
    # Initialize the derivative array excluding the last term (constant term)
    derivative = np.array(coefficients[:-1], dtype=np.longdouble)
    
    # Calculate the coefficients of the derivative
    for i in range(n - 1):
        # Multiply each coefficient by its power (n-1-i) to get the derivative coefficient
        derivative[i] *= (n - 1 - i)
    
    return derivative


def calculate_poly_derivative_normalized(coefficients):
    """
    Calculate the derivative of a polynomial and normalize it based on the highest power.

    :param coefficients: The coefficients of the original polynomial.
    :return: The coefficients of the normalized derivative.
    """
    # Length of the coefficients array
    n = len(coefficients)
    
    # Initialize the derivative array excluding the last term (constant term)
    derivative = np.array(coefficients[:-1], dtype=np.longdouble)
    
    # Calculate the coefficients of the normalized derivative
    for i in range(n - 1):
        # Multiply each coefficient by its power and normalize it by the highest power (n-1)
        derivative[i] *= (n - 1 - i) / (n - 1)
    
    return derivative


def calculate_polynomial(coefficients, x):
    """
    Calculates the polynomial coefficient at a given point x.
    param coefficients: The coefficients of the polynomial.
    param x: The point at which to calculate the coefficient.
    return: The calculated polynomial coefficient.
    """
    result = np.longdouble(0)  # Initialize the result variable to 0
    # Iterate over each coefficient in the original order
    for coefficient in coefficients:
        # Multiply the result by x and add the current coefficient
        result = result * x + coefficient
    return result


def remove_similar_roots(roots, epsilon):
    """
    Removes similar roots (their difference is less than epsilon) and returns a sorted list of roots.

    :param roots: The list of roots.
    :param epsilon: The threshold for similarity.
    :return: The sorted list of distinct roots.
    """
    # If the list of roots is empty, return an empty list
    if len(roots) == 0:
        return []
    
    # Sort the roots to ensure they are in ascending order
    roots.sort()   #O(nlog(n))

    # Initialize the result list with the first root
    result = [roots[0]]
    
    # Iterate over the roots
    for root in roots:
        # Compare each root with the last root in the result list
        # If the difference is greater than or equal to epsilon, add it to the result list
        if abs(root - result[-1]) >= epsilon:
            result.append(root)
    return result


def polynomial_roots_finder(coefficients, epsilon):
    """
    Find all roots of a polynomial.

    :param coefficients: List of coefficients of the polynomial in descending order of powers (e.g., [a_n, a_(n-1), ..., a_1, a_0])
    :param epsilon: Tolerance for root-finding methods, controls the precision of root approximations
    :return: List of roots found in the given polynomial
    """
    if len(coefficients) == 2:
        # If it's a linear equation (ax + b), return the root (-b/a)
        return [-coefficients[1] / coefficients[0]]

    result_roots = []

    # Find roots of the derivative of the polynomial
    derivative_roots = polynomial_roots_finder(calculate_poly_derivative_normalized(coefficients), epsilon)
    if len(derivative_roots) == 0:
        derivative_roots = [0]

    #Using Lagrange's Bound provides a  way to estimate the maximum possible magnitude of any root of a polynomial
    max_coefficient = max(abs(coefficient) for coefficient in coefficients[1:])
    
    # Estimate the maximum distance to search based 
    max_distance = max_coefficient / abs(coefficients[0]) + 1
    
    left_bound= max_distance*(-1)
    right_bound= max_distance

    derivative_roots = [left_bound] + derivative_roots
    derivative_roots.append(right_bound)

    for i in range(len(derivative_roots) - 1):
        a, b = derivative_roots[i], derivative_roots[i + 1]
        a_value = calculate_polynomial(coefficients, a)
        b_value = calculate_polynomial(coefficients, b)
        if getSign(a_value) != getSign(b_value):
            # Apply root-finding methods within the interval
            root = root_finder_newton_raphson(coefficients, (a + b) / 2, epsilon, (b - a) / 2)
            if root is None or root < a or root > b:
                root = root_finder_bisection(coefficients, a, b, epsilon)
            result_roots.append(root)

    # Remove similar roots and sort the results
    return remove_similar_roots(result_roots, 1.0e-6)


#main()
def main():
    # Constants and data initialization
    roots_of_numpy = []
    roots_of_the_polynomial = []
    polynomial_coefficients = []
    epsilon = 1.0e-9

    # Read the polynomial coefficients from the file
    with open('poly_coeff(997).txt') as f:
        polynomial = f.read().splitlines()
    polynomial_coefficients = np.array(polynomial, np.longdouble)

    # Using numpy method to find the roots
    start = time.time()
    roots_of_numpy = np.roots(np.array(polynomial, float))
    end = time.time()

    # Print the roots from np.roots() library method
    print("\nThe roots from np.roots() library method are:\n")
    for root1 in roots_of_numpy:
        if np.isreal(root1):
          print("(" + str(root1.real) + ")")

    # Print the runtime for np.roots() method
    print("\nnp.roots() runtime is: %s seconds\n" % (end - start))

    # Our implementation to find the roots
    start = time.time()
    roots_of_the_polynomial = polynomial_roots_finder(polynomial_coefficients, epsilon)
    end = time.time()

    # Print the roots from our implementation
    print("The roots from our implementation are:\n")
    for root1 in roots_of_the_polynomial:
        print("(" + str(root1) + ")")

    # Print the runtime for our implementation
    print("\nThe runtime is: %s seconds" % (end - start))

    #print("\nThe np.polyval F(x=-5.919278610574727)= "+str(np.polyval(polynomial_coefficients,-5.919278610574727 )))
    #print("\nThe F(x=-0.28564650800220487)= "+str(calculate_polynomial(polynomial_coefficients,-0.28564650800220487 )))
    #print("\nThe F(x=-0.8315492003888125)= "+str(calculate_polynomial(polynomial_coefficients,-0.8315492003888125 )))


if __name__ == '__main__':
     main()