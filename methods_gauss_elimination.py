# -*- coding: utf-8 -*-
"""
@date:        Created on Tue Feb 16 14:26:40 2016
@author:      Jose Ortiz Costa, Group #7
@purpose:     Math 400, Assigment #1, Problem #3
@description: This file contains useful methods to perform numerical analysis
              in matrices using the following methods:
              0. Gauss-Naive
              1. Partial Pivotiong
              2. Scaled Partial Pivoting
              3. LU factorization
              4. Gauss-Seider 


"""
# Libraries containing the methods supported for this file in order to find
# solutions for matrices

import math
import numpy as np
from scipy.integrate import quad



# Enumerator that enumerates the Gauss methods supported in this file
class Method:
    NAIVE = 0
    PARTIAL = 1
    SCALED_PARTIAL = 2 
    LU_FACTORIZATION = 3
    GAUSS_SEIDEL = 4
     

# Retursn a exact copy of the matrix given as argument
def copyMatrix(M):
    return([[M[row][col] for col in range(cols(M))]for row in
            range(rows(M))])

# returns a vector of maximum pivots per row from a matrix
# e.g: given a matrix [[4,3,2], [6,7,4], [9,8,10]]
# it returns v = [4,7,10]
def vector_of_pivots (mat):
    "returns a vector of pivots"
    return [max(row) for row in mat]

# returns number of rows from a matrix 
def rows(mat):
    "return number of rows"
    return(len(mat))

# returns number of columns from a matrix 
def cols(mat):
    "return number of cols"
    return(len(mat[0]))

# returns a column from a matrix    
def GetCol(mat, col):
    "retun a column from a matrix"
    if col >= cols(mat) or col < 0:
       return "Error: col index out of range"
    return [mat[row][col] for row in range(rows(mat))];

# Gets a row from a matrix    
def GetRow(mat, row):
    "return row row from matrix mat"
    if row >= rows(mat) or row < 0:
       return "Error: row index out of range"
    return(mat[row])

# Adds a column 
def addColumn (M,col):
    i = 0
    for row in M:
        row.append(col[i][0])
        i+=1
    return M
    
# find the position of the greater pivot
def pivot_position (mat, rowIndex):
    row = GetRow(mat, rowIndex)
    return (row.index(max(row)))

# Interchanges row in a matrix given the matrix and rows indexes
# Returns the matrix with the rows interchanged.
def swap_row (mat, row1Index, row2Index):
    M = copyMatrix(mat)
    M[row1Index], M[row2Index] = M[row2Index], M[row1Index]
    return M;
    

# Replace a column in a matrix by another given column
def replace_mat_column (mat, newCol, colIndex):
    M = copyMatrix(mat) 
    C = newCol
    if (colIndex >= len(newCol)):
        return "Error: colIndex cannot be equal or greater than size of newCol"
    for rowIndex in range(len(M)):
        for tmpColIndex in range(len(M[rowIndex])):
            M[rowIndex][colIndex] = C[rowIndex]
    return M;

# Find pivot for naive Gauss 
def findPivotrow1(mat,col):
#    Finds index of the first row with nonzero entry on or
#    below diagonal.  If there isn't one return(-1).

    epsilon = 10**(-17)
    for row in range(col, rows(mat)):
#        if mat[row][col] != 0:
        if abs(mat[row][col]) > epsilon:
            return(row)
    return(-1)

# add two vectors
def addVectors(A,B):
    "add two vectors"
    if len(A) != len(B):
        print("addVectors: different lengths")
        return()
    return([A[i]+B[i] for i in range(len(A))])

# add rows
def addrows(M, f, t, scale=1):
    "add scale times row f to row t"
    N=copyMatrix(M)
    T=addVectors(scalarMult(scale,N[f]),N[t])
    N[t]=T
    return(N)  

# Test if V is a vector
def vectorQ(V):
    "mild test to see if V is a vector"
    if type(V) != type([1]):
        return(False)
    if type(V[0]) == type([1]):
        return(False)
    return(True)

# Performs scalar multipliation   
def scalarMult(a,mat):
    "multiply a scalar times a matrix"
    if vectorQ(mat):
        return([a*m for m in mat])
    for row in range(rows(mat)):
        for col in range(cols(mat)):
            mat[row][col] = a*mat[row][col]
    return(mat)
# Calculates the dot product  
def Dot(S,T):
    "return dot product of two vectors"
    if len(S) != len(T):
        return "Error: Vectors need to be same size"
    return sum([S[index] * T[index] for index in range(len(S))])

# Multiply two matrices
def MultiplyMat(mat1,mat2):
    "multiply two matrices"
    if cols(mat1) != rows(mat2):
        return "Error: number of columns in mat1 must be equal to the number of rows in mat2"
    return [[Dot(GetRow(mat1,row),GetCol(mat2,col)) for col in range(cols(mat2))] for row in range(rows(mat1))] 

# Add two vectors
def AddVectors(S,T):
    "return S+T"
    if len(S) != len(T):
       return "Error: Vectors need to be same size"
    return [S[index] + T[index] for index in range(len(S))]
    
#Add two matrices
def AddMatrices(A,B):
    "add two matrices"
    return [AddVectors(A[vectors], B[vectors]) for vectors in range(rows(A))]
    
# reduce rows in a matrix   
def rowReduce(M):
    "returns a row reduced version of M"
    N = copyMatrix(M)
    cs = cols(M)-2   # no need to consider last two cols
    rs = rows(M)
    for col in range(cs+1):
        j = findPivotrow1(N,col)
        if j < 0:
            print("\nrowReduce: No pivot found for column index %d "%(col))
            return(N)
        else:
            if j != col:
                N = swap_row(N,col,j)
            scale = -1.0 / N[col][col]
            for row in range(col+1,rs):                
                N=addrows(N, col, row, scale * N[row][col])
    return(N)

# reduce a single row in a matrix
def getRowReduced (mat,rowToReduceIndex,pivotRowIndex):
    "returns a row reduced from a given matrix"
    pivot_row = GetRow(mat, pivotRowIndex)
    row = GetRow(mat, rowToReduceIndex)
    ratio = row[pivotRowIndex] / float(pivot_row[pivotRowIndex])
    for tmpRowIndex in range(len(row)):
        egv = (ratio * pivot_row[tmpRowIndex]) - row[tmpRowIndex]
        row[tmpRowIndex] = egv
    return row
    
#   given a row reduced augmented matrix with nonzero 
#   diagonal entries, returns a solution vector
def backSub(M, solInMatrixForm):
    "Returns the x's solution of a matrix using back sustitution method"
    cs = cols(M)-1 # cols not counting augmented col
    sol = [0 for i in range(cs)] # place for solution vector
    for i in range(1,cs+1):
        row = cs-i # work backwards
        sol[row] = ((M[row][cs] - sum([M[row][j]*sol[j] for
                    j in range(row+1,cs)])) / M[row][row])
    if (solInMatrixForm):
        solm = [[sol[nrow] for __ in range(0,1)] for nrow in range (len(sol))]
        return(solm)
    return (sol)
    
# Partial pivoting for Gauss elimination
def partial_pivoting (mat):
    "Returns a reduced matrix using the scaled partial pivoting method"
    M = copyMatrix(mat)
    for colIndex in range(len(M)-1):
        M.sort(reverse=True)        
        for partialRowIndex in range(colIndex, len(M) - 1):
            row = getRowReduced(M,partialRowIndex+1,colIndex)
            M[partialRowIndex+1] = row # add row back to matrix
    return M

# Scaled partial pivoting for Gauss elimination
def scaled_partial_pivoting (mat):
    "Returns a reduced matrix using the scaled partial pivoting method"
    M = copyMatrix(mat)
    max_pivots = vector_of_pivots(M) # gets greater value from every row in mat
    for colIndex in range(len(M)-1):
        col = GetCol(M, colIndex)
        vec = [col[indexOfPivots] / float(max_pivots[indexOfPivots]) for indexOfPivots in range(len(max_pivots))]
        pivot = vec.index(max(vec)) + colIndex
        M[colIndex], M[pivot] = M[pivot], M[colIndex]
        for partialRowIndex in range(colIndex, len(M) - 1):
            row = getRowReduced(M,partialRowIndex+1,colIndex)
            M[partialRowIndex+1] = row # add row back to matrix
    return M

# Performs LU Factorization in a matrix
def LUFactorization (M):
    "Needs to be implemented"
    return M
    

# Gauss Seider aproximation
def gauss_seidel(A, b, x=None, n=1000):
    if x==None:
        x = [0 for __ in A]
    b = sum(b, [])
    L = np.tril(A)
    U = A - L
    tmpx= []
    for i in range(n):
        x = np.dot(np.linalg.inv(L), b - np.dot(U, x))
        if set(tmpx) == set(x) :
            break;
        tmpx = x
    return [[x[row] for __ in range(0,1)] for row in range(0,len(x))]  
   

# Construct the b column of a Hilbert Matrix  
# integral of x^n(c)sin(xpi/2)dx
def Hilbert_b_integral(x,n):
    "Returns the b part of a Hilbert matrix calculated by an integral"
    return math.pow(x,n)*math.sin((math.pi*x)/2)

         
# Construct the Hilbert matrix for a given N
# [[1...1/n][....][....][1/n ..... 1/(2n-1)]]
# This matrix also include the b column
def Hilbert_matrix(n):
    "Returns a Hilbert Matrix contructed by a given n"
    H = [[0 for __ in range(0,n+1)] for __ in range (0,n+1)]
    for row in range(0,len(H)):
        for col in range(0,len(H[row])):
            if (col != len(H)-1):
               H[row][col] = 1/float(col+row+1)
            else:
               H[row][col] = quad(Hilbert_b_integral, 0, 1, args=(row))[0]
    #return [[ 1/float(col+row+1) for col in range(0,n)] for row in range(0,n)]
    return H

# Construct b solution of a Hilbert matrix 
def Hilbert_matrix_b(n):
    b = [[0 for __ in range(0,1)] for __ in range(0,n)]
    for row in range(len(b)):
        for col in range(len(b[row])):
            b[row][col] = quad(Hilbert_b_integral, 0, 1, args=(row))[0]
    return b
    
# Construct a Hilbert Matrix without b
def Hilbert_matrix_without_b(n):
    "Returns a Hilbert Matrix contructed by a given n"
    H = [[0 for __ in range(0,n)] for __ in range (0,n)]
    for row in range(0,n):
        for col in range(0,n):
            H[row][col] = 1/float(col+row+1)
    return H

# Solves a matrix by a gauss method   
def solve_matrix (Hb, method, b, H=None, x=None, n=100):
    if (method == 0):
        gauss_naive_sol = backSub(rowReduce(Hb),1)
        r = residual_error(H,gauss_naive_sol,b)
        rnorm = norm_residual_error(r)
        print ("Method: Gauss-Naive\n")
        printSolutions(Hb,H,gauss_naive_sol,b,r,rnorm) # naive
    elif (method == 1):
        gauss_partial_sol = backSub(partial_pivoting(Hb),1)
        r = residual_error(H,gauss_partial_sol,b)
        rnorm = norm_residual_error(r)
        print ("Method: Gauss-Partial-Pivoting\n")
        printSolutions(Hb,H,gauss_partial_sol,b,r,rnorm) # partial
    elif (method == 2):
        gauss_scaled_partial_pivoting = backSub(scaled_partial_pivoting(Hb),1)
        r = residual_error(H,gauss_scaled_partial_pivoting,b)
        rnorm = norm_residual_error(r)
        print("Method: Gauss-Scaled-Partial-Pivoting\n")
        printSolutions(Hb, H,gauss_scaled_partial_pivoting,b,r,rnorm) # scaled partial
    elif (method == 3):
        print "\nLU FACTORIZATION NOT IMPLEMENTED YET\n"
        L = LU_factorization(H)[0] # gets L
        U = LU_factorization(H)[1] # gets U
        lu_sol = solveLU(L,U,b)
        r = residual_error(H,lu_sol,b)
        rnorm = norm_residual_error(r)
        printSolutions(Hb, H,lu_sol,b,r,rnorm) # scaled partial
    else:
        gs = gauss_seidel(H, b)
        r = residual_error(H,gs,b)
        rnorm = norm_residual_error(r)
        print("Method: Gauss-Seidel\n")
        printSolutions(Hb, H,gs,b,r,rnorm) 

# calculates the residual error vector where r = Ax - b
# and returns ||r||. That's it the norm of the r vector
# which represents the residual vector of the solutions
# A = matrix A
# xap = matrix of aproximate solutions of x1,x2...
# b = solution to the matrix
def residual_error(A,x,b):
    Ax = MultiplyMat(A,x) # Ax
    r = []
    for row in range(0,len(Ax)):
        for col in range(len(Ax[row])):
            r.append(b[row][col] - Ax[row][col])
    return r

#Given a residual error vector calculates its norm   
def norm_residual_error (r):
    return math.sqrt(sum(abs(row) for row in r)) # return r     


# Prints solutions in a nicely format
def printSolutions(Hb,H,sol,b,r,rnorm):
    x = MultiplyMat(H,sol)
    print '{:<15} {:>35} {:>35} {:>53}'.format("Hilbert X solutions","Hilbert b solutions", "Compute Ax", "Residual Error Vector (r = b - Ax)")
    for nrow in range(len(sol)):
        for ncol in range(len(sol[nrow])):
            tmpnrow = nrow + 1
            print '{:>1} {:<30} {:<3} {:>12} {:>24} {:>12} {:>24} {:>1}'.format("x"+`tmpnrow`+" =", sol[nrow][ncol], 
                  "b"+`tmpnrow`+" =",b[nrow][ncol], "Ax"+`tmpnrow`+" =",x[nrow],"r"+`tmpnrow`+" =",r[nrow])
    print ("\nNorm of the residual error vector ||r|| = %e\n\n" %rnorm)

# multipies matrix from a given vector   
def matrixMul(A, B):
    TB = zip(*B)
    return [[sum(ea*eb for ea,eb in zip(a,b)) for b in TB] for a in A]

# pivotizes for LU factorization
def pivotize(m):
    """Creates the pivoting matrix for m."""
    n = len(m)
    ID = [[float(i == j) for i in xrange(n)] for j in xrange(n)]
    for j in xrange(n):
        row = max(xrange(j, n), key=lambda i: abs(m[i][j]))
        if j != row:
            ID[j], ID[row] = ID[row], ID[j]
    return ID

# performs LU factorization and returns L,U and P
def LU_factorization(A):
    """Decomposes a nxn matrix A by PA=LU and returns L, U and P."""
    n = len(A)
    L = [[0.0] * n for i in xrange(n)]
    U = [[0.0] * n for i in xrange(n)]
    P = pivotize(A)
    A2 = matrixMul(P, A)
    for j in xrange(n):
        L[j][j] = 1.0
        for i in xrange(j+1):
            s1 = sum(U[k][j] * L[i][k] for k in xrange(i))
            U[i][j] = A2[i][j] - s1
        for i in xrange(j, n):
            s2 = sum(U[k][j] * L[i][k] for k in xrange(j))
            L[i][j] = (A2[i][j] - s2) / U[j][j]
    return (L, U, P)


# Solves x's solutions from a matrix previosly descomposed by LU factorization        
def solveLU(L,U,b):
    A = MultiplyMat(L,U)
    Ax = addColumn(A,b)
    return backSub(rowReduce(Ax),1)
       

"----------- TEST FUNCTIONS ------------------------------------------------"

for it in range(3,25):
    print ("\n\n************ HILBERT MATRIX TEST ITERATION = %d *******************\n" %it)
    Hb = Hilbert_matrix(it) # creates H with b
    H = Hilbert_matrix_without_b(it) # creates H without b
    b = Hilbert_matrix_b(it) # creates b
    solve_matrix(Hb,Method.NAIVE, b,H) # solve H by Gauss-Naive
    solve_matrix(Hb,Method.SCALED_PARTIAL, b,H) # solve H by Scaled Partial
    solve_matrix(Hb,Method.LU_FACTORIZATION, b,H) # solve H by LU Factorization
    solve_matrix(Hb,Method.GAUSS_SEIDEL, b,H) # Solve H by Gauss-Seidel 
    
    
"------------------------------- END OF TESTING ----------------------------"
  
    
    
    
    
   






      
    



  
 

      
        
        

    



