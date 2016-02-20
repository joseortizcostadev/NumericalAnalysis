# -*- coding: utf-8 -*-
"""
@date:        Created on Tue Feb 16 14:26:40 2016
@author:      Jose Ortiz Costa, Group #7
@purpose:     Math 400, Assigment #1, Problem #3
@description: This file contains useful methods to perform numerical analysis
              in matrices using the following methods:
              1. Partial Pivotiong
              2. Scaled Partial Pivoting
              3. LU factorization
              4. Gauss-Seider 


"""
# Libraries containing the methods supported for this file in order to find
# solutions for matrices
import numpy as np
import math
from scipy.integrate import quad

# Enumerator that enumerates the Gauss methods supported in this file
class Method:
    NAIVE, PARTIAL, SCALED_PARTIAL, LU_FACTORIZATION, GAUSS_SEIDER = range(0, 5) 

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
def backSub(M):
    "Returns the x's solution of a matrix using back sustitution method"
    cs = cols(M)-1 # cols not counting augmented col
    sol = [0 for i in range(cs)] # place for solution vector
    for i in range(1,cs+1):
        row = cs-i # work backwards
        sol[row] = ((M[row][cs] - sum([M[row][j]*sol[j] for
                    j in range(row+1,cs)])) / M[row][row]) 
    return(sol)
    
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
def Gauss_Seider(M, b, x, n):
    L = np.tril(M)
    U = A - L
    for k in range(n):
        x = np.dot(np.linalg.inv(L), b - np.dot(U, x))
    return x

# Construct the b column of a Hilbert Matrix   
def Hilbert_b(x,n):
    "Returns the b part of a Hilbert matrix calculated by an integral"
    return math.pow(x,n)*math.sin((math.pi*x)/2)

# Construct the Hilbert matrix for a given N
def Hilbert_matrix(n):
    "Returns a Hilbert Matrix contructed by a given n"
    H = [[0 for col in range(0,n+1)] for row in range(0,n)]
    for row in range(0,n):
        for col in range(0,n+1):
            if (col != n):
                H[row][col] = (1/float(col+row+1))
            else:
                H[row][col] = quad(Hilbert_b, 0, 1, args=(row))[0]
    return H

   

A = np.array([[25,5,1,106.8],[64,8,1,177.2],[144,12,1,279.2]])
b = [1,2,3]
x = [1, 1, 1]
B = [[25,5,1,106.8], [64,8,1,177.2],[144,12,1,279.2]]
n = 20

def solve_matrix (M, method):
    if (method == 0):
        print "\nNaive Gauss Elimination\n"
        print backSub(rowReduce(B)) # naive
    elif (method == 1):
        print "\nPartial Pivoting Gauss Elimination\n"
        print backSub(partial_pivoting(B)) # partial
    elif (method == 2):
        print "\nScaled Partial Pivoting Gauss Elimination\n"
        print backSub(scaled_partial_pivoting(B)) # scaled partial
    elif (method == 3):
        print "\nLU FACTORIZATION NOT IMPLEMENTED YET\n"
        #print backSub(scaled_partial_pivoting(B)) # scaled partial
    else:
        print "\nError aproximation in Gauss-Seider\n"
        print Gauss_Seider(A, b, x, n)
        
        

"----------- TEST FUNCTIONS ------------"

solve_matrix(B, Method.NAIVE)
solve_matrix(B, Method.PARTIAL)
solve_matrix(B, Method.SCALED_PARTIAL)
#solve_matrix(B, Method.GAUSS_SEIDER)
print "\n Hilbert Matrix with n=3\n"
print Hilbert_matrix(3) # construct a Hilber matrix of H(n)




      
    



  
 

      
        
        

    



