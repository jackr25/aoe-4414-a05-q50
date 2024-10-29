# ray_ellipsoid_intersection.py
# Usage: python3 ray_ellipsoid_intersection.py d_l_x d_l_y d_l_z c_l_x c_l_y c_l_z
#  Text explaining script usage
# Parameters:
# d_l_x: x-component of origin-referenced ray direction
# d_l_y: y-component of origin-referenced ray direction
# d_l_z: z-component of origin-referenced ray direction
# c_l_x: x-component offset of ray origin
# c_l_y: y-component offset of ray origin
# c_l_z: z-component offset of ray origin
# Output:
#  A description of the script output
# l_d[0]: x-component of intersection point
# l_d[1]: y-component of intersection point
# l_d[2]: z-component of intersection point
# Written by Jack Rathert
# Other contributors: None
#
## imports---------------------------------------------------------------------------
import sys
import math as m

## constants variable input----------------------------------------------------------
E_E = 0.081819221456
R_E = 6378.1363 # km

## Classes---------------------------------------------------------------------------
class numpy_lite:
    def __init__(self):
        from math import sqrt
        self.sqrt = sqrt 
    def matrix_mult(self, matrix1, matrix2):
        # Wrap flat lists into 1xN matrices
        if isinstance(matrix1[0], (int, float)):
            matrix1 = [[elem] for elem in matrix1]  # Convert to Nx1 matrix
        if isinstance(matrix2[0], (int, float)):
            matrix2 = [[elem] for elem in matrix2]  # Convert to Nx1 matrix
        # Get matrix dimensions
        rowsA, colsA = len(matrix1), len(matrix1[0])
        rowsB, colsB = len(matrix2), len(matrix2[0])
        # Nx1 * 1xN -> NxN (outer product)
        if colsA == 1 and rowsB == 1:
            return [[matrix1[i][0] * matrix2[0][j] for j in range(colsB)] for i in range(rowsA)]
        # 1xN * Nx1 -> scalar (dot product)
        if rowsA == 1 and colsB == 1 and colsA == rowsB:
            return sum(matrix1[0][i] * matrix2[i][0] for i in range(colsA))
        # General matrix multiplication
        if colsA != rowsB:
            raise ValueError(f"Incompatible dimensions: {colsA} != {rowsB}")
        # Perform matrix multiplication
        return [[sum(matrix1[i][k] * matrix2[k][j] for k in range(colsA)) for j in range(colsB)] for i in range(rowsA)]
    def matrix_add(self,list1,list2):
      return [x+y for x,y in zip(list1, list2)]
    def matrix_sub(self,list1,list2):
      return [x-y for x,y in zip(list1,list2)]
    def mag(self,list1):
      return sqrt(sum(x**2 for x in list1)) # type: ignore
    def smul(self,scalar, vector):
      return [scalar * x for x in vector]
    def vecadd(self,vector1,vector2):
      if len(vector1) != len(vector2):
         raise ValueError("Vectors must be of the same length.")
      return [x + y for x, y in zip(vector1, vector2)]
    def vecsub(self,vector1,vector2):
      if len(vector1) != len(vector2):
         raise ValueError("Vectors must be of the same length.")
      return [x - y for x, y in zip(vector1, vector2)]  
    def dotprod(self,vector1,vector2):
      if len(vector1) != len(vector2):
         raise ValueError("Vectors must be of the same length.")
      return sum(x * y for x, y in zip(vector1, vector2))
nplite = numpy_lite()

## Functions--------------------------------------------------------------------------
def d_rayvec(d_x:float,d_y:float, d_z:float, c_x:float, c_y:float, c_z:float):
  a = d_x**2 + d_y**2+(d_z**2)/(1-E_E**2)
  b = 2*(d_x*c_x+d_y*c_y+(d_z*c_z)/(1-E_E**2))
  c = c_x**2+c_y**2+(c_z**2)/(1-E_E**2) - R_E**2
  discrim = b**2 - 4*a*c
  if discrim<0:
   raise ValueError("Discriminate is Negative")
  d1 = (-b + m.sqrt(discrim)) / (2 * a)
  d2 = (-b - m.sqrt(discrim)) / (2 * a)
  if d1 > 0 and d2 > 0: # Choose the smallest positive root
   d = min(d1, d2)  
  elif d1 > 0:
   d = d1  # Only d1 is positive 
  elif d2 > 0:
   d = d2  # Only d2 is positive
  else:
   raise ValueError("Both intersections are behind the ray origin.")
  return d

def rayepintersect(c_l:float, d_l:float, d:float):
  l_d = nplite.vecadd(nplite.smul(d,d_l),c_l)
  return l_d
# Arguments---------------------------------------------------------------------------
# parsing
if len(sys.argv)==7:
   d_l_x = float(sys.argv[1])
   d_l_y = float(sys.argv[2])
   d_l_z = float(sys.argv[3])
   c_l_x = float(sys.argv[4])
   c_l_y = float(sys.argv[5])
   c_l_z = float(sys.argv[6])
else:
  print(\
   'Usage: '\
   'python3 ray_ellipsoid_intersection.py d_l_x d_l_y d_l_z c_l_x c_l_y c_l_z'\
  )
  exit()

## Function Calls ----------------------------------------------------------------------
d = d_rayvec(d_l_x,d_l_y,d_l_z,c_l_x,c_l_y,c_l_z)
c_l = [c_l_x,c_l_y,c_l_z]
d_l = [d_l_x,d_l_y,d_l_z]
l_d = rayepintersect(c_l,d_l,d)

print(l_d[0]) # x-component of intersection point
print(l_d[1]) # y-component of intersection point
print(l_d[2]) # z-component of intersection point