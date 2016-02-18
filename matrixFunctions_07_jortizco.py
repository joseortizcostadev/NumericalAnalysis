"""
matrixFunctions_07_jortizco.py

Student: Jose Ortiz Costa

2/1/16 bic M400

This is a skeleton for the matrix part of first Math 400 assignment

"""

#######  START Administrivia 
m400group = 07   # change this to your group number

m400names = ['Emily Conway', 'Samuel Gluss', '*Jose Ortiz Costa', 'Abdelmajid Samir','Edward Yao'] # change this for your names

def printNames():
    print("matrixFunctions_07_jortizco.py for group %s:"%(m400group)),
    for name in m400names:
        print("%s, "%(name)),
    print

printNames()

#######  END Administrivia


"""
Vector Functions

copy these three functions from your finished vectorFunctions.py file

"""

def ScalarMult(s,V):
    "return vector sV"
    return [s*item for item in V]

def AddVectors(S,T):
    "return S+T"
    if len(S) != len(T):
       return "Error: Vectors need to be same size"
    return [S[index] + T[index] for index in range(len(S))]
        
           
def Dot(S,T):
    "return dot product of two vectors"
    if len(S) != len(T):
        return "Error: Vectors need to be same size"
    return sum([S[index] * T[index] for index in range(len(S))])

"""

Matrix Functions

"""

def showMatrix(mat):
    "Print out matrix"
    for row in mat:
        print(row)

def rows(mat):
    "return number of rows"
    return(len(mat))

def cols(mat):
    "return number of cols"
    return(len(mat[0]))
 

#### Functions for you to finish

def GetCol(mat, col):
    "Retun a column from a matrix"
    if col >= cols(mat) or col < 0:
       return "Error: col index out of range"
    return [mat[row][col] for row in range(rows(mat))]

def Transpose(mat):
    "return transpose of mat"
    if not mat: return []
    return [[mat[row][column] for row in range(rows(mat))] for column in range(cols(mat))]
         
def GetRow(mat, row):
    "return row row from matrix mat"
    if row >= rows(mat) or row < 0:
       return "Error: row index out of range"
    return(mat[row])

def ScalarMultMatrix(a,mat):
    "return scalar multiplication of a matrix"
    return[ScalarMult(a,mat[vector]) for vector in range(rows(mat))]

def AddMatrices(A,B):
    "add two matrices"
    return [AddVectors(A[vectors], B[vectors]) for vectors in range(rows(A))]
       
def MultiplyMat(mat1,mat2):
    "multiply two matrices"
    if cols(mat1) != rows(mat2):
        return "Error: number of columns in mat1 must be equal to the number of rows in mat2"
    return [[Dot(GetRow(mat1,row),GetCol(mat2,col)) for col in range(cols(mat2))] for row in range(rows(mat1))]        
             
                

######  Initial tests

A= [[4,-2,1,11],
    [-2,4,-2,-16],
    [1,-2,4,17]]

Ae= [[4,-2,1],
    [-2,4,-2],
    [1,-2,4]]


Bv=[11,-16,17]

Bm=[[11,-16,17]]

C=[2,3,5]

print("running matrixFunction_07_jortizco.py file")

def testMatrix():
    print("A")
    showMatrix(A)
    print("Bm")
    showMatrix(Bm)
    print("Ae")
    showMatrix(Ae)
    print("multiplyMat(Ae,A)")
    showMatrix(MultiplyMat(Ae,A))
    print("scalarMultMatrix(2,A))")
    showMatrix(ScalarMultMatrix(2,A))
    print("addMatrices(A,A)")
    showMatrix(AddMatrices(A,A))
    print("transpose(A)")
    showMatrix(Transpose(A))

###  uncomment next line to run initial tests 
testMatrix()