import sys
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.pyplot as plt

def add(a, b):
    res = [0.0]*3
    for i in range(0, 3):
       res[i] = a[i] + b[i]
    return res

def subtract(a, b):
    res = [0.0]*3
    for i in range(0, 3):
       res[i] = a[i] - b[i]
    return res

def minus(a):
    res = [0.0]*3
    for i in range (0, 3):
        res[i] = -a[i]
    return res

def vector_product(a, b):
   res = [a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0]]
   return res

def scalar_product(a, b):
    res = 0.0
    for i in range(len(a)):
        res += a[i]*b[i]
    return res

def eq_signs(a):
    if (a[0]<0):
        for el in a:
            if (el>0 or el == 0):
                return False
    elif (a[0] > 0):
        for el in a:
            if (el<0 or el == 0):
                return False
    else:
        return False
    return True

def zeroDistance(d, eps):
    for el in d:
        if (abs(el) >= eps):
            return False
    return True

def sign(a):
    if (a>=0):
        return 1
    elif (a<0):
        return -1

def ZeroSign3D(a):
    npos = sum(x > 0 for x in a)
    nneg = sum(x < 0 for x in a)
    if (npos == 2):
        return False
    if (nneg == 2):
        return True

def signZero(a, pos):
    if (pos):
        if (a>=0):
            return 1
        elif (a<0):
            return -1
    else:
        if (a>0):
            return 1
        elif (a<=0):
            return -1
       
def DiffSigns(a, eps):
    numnegative = 0
    numpositive = 0
    for el in a:
        if (el > eps):
            numpositive += 1
        elif (el < -eps):
            numnegative += 1
    return numpositive != 0 and numnegative != 0

def RightOrder2D(Int):
    if (Int[0]>Int[1]):
        temp = Int[0]
        Int[0]=Int[1]
        Int[1]=temp

def swap(a, b, i, j):
    tempa = a[i]
    tempb = b[i]
    a[i] = a[j]
    b[i] = b[j]
    a[j] = tempa
    b[j] = tempb

#rightform: the first and the last vertice make d of the same sign to another plane
def rightform(a, d):
    pos = ZeroSign3D(d)
    if (signZero(d[0], pos) != signZero(d[2], pos)):
        if (signZero(d[0], pos) == signZero(d[1], pos)):
            swap(a, d, 1, 2)
        else:
            swap(a, d, 0, 1)

def IntIntersects(T1, T2):
    if (T2[0] <= T1[1] and T1[0] <= T2[1]):
        return True
    else:
        return False

def Sign2D(p1, p2, p3):
    return (p1[0] - p3[0]) * (p2[1] - p3[1]) - (p2[0] - p3[0]) * (p1[1] - p3[1])

def PointInTriangle2D(pt, p1, p2, p3, eps):

    b1 = Sign2D(pt, p1, p2) < eps;
    b2 = Sign2D(pt, p2, p3) < eps;
    b3 = Sign2D(pt, p3, p1) < eps;
    a = (b1 and b2) or (not b1 and not b2)
    return (b1 == b2) and (b2 == b3)
    #return ((b1 and b2) or (not b1 and not b2)) and ((b3 and b2) or (not b3 and not b2))

def PointInTriangle2DBar(pt, p1, p2, p3, eps):
     denominator = ((p3[1] - p1[1])*(p2[0] - p1[0]) - (p3[0] - p1[0])*(p2[1] - p1[1]))
     if (abs(denominator) < eps):
         return False
     else:
        a = ((p3[1] - p1[1])*(pt[0] - p1[0]) - (p3[0] - p1[0])*(pt[1] - p1[1])) / denominator
        b = ((pt[1] - p1[1])*(p2[0] - p1[0]) - (pt[0] - p1[0])*(p2[1] - p1[1])) / denominator
        c = 1 - a - b
        return 0 - eps <= a and a <= 1 + eps and 0 - eps <= b and b <= 1 + eps and 0 - eps <= c and c <= 1 + eps

def EdgesOfTriangle(p1, p2, p3):
    edge1 = [p1[1] - p2[1], p2[0]-p1[0], p1[0]*p2[1] - p2[0]*p1[1]]
    edge2 = [p1[1] - p3[1], p3[0]-p1[0], p1[0]*p3[1] - p3[0]*p1[1]]
    edge3 = [p3[1] - p2[1], p2[0]-p3[0], p3[0]*p2[1] - p2[0]*p3[1]]
    return [edge1, edge2, edge3]

def TrianglesIntersect2D(first, second, eps):
    #if (PointInTriangle2D(first[0], second[0], second[1], second[2], eps) or PointInTriangle2D(first[1], second[0], second[1], second[2], eps) or PointInTriangle2D(first[2], second[0], second[1], second[2], eps) or PointInTriangle2D(second[0], first[0], first[1], first[2], eps) or PointInTriangle2D(second[1], first[0], first[1], first[2], eps) or PointInTriangle2D(second[2], first[0], first[1], first[2], eps)):
    if (PointInTriangle2DBar(first[0], second[0], second[1], second[2], eps) or PointInTriangle2DBar(first[1], second[0], second[1], second[2], eps) or PointInTriangle2DBar(first[2], second[0], second[1], second[2], eps) or PointInTriangle2DBar(second[0], first[0], first[1], first[2], eps) or PointInTriangle2DBar(second[1], first[0], first[1], first[2], eps) or PointInTriangle2DBar(second[2], first[0], first[1], first[2], eps)):
        return True
    else:
        dist = [[0.0 for j in range(3)] for i in range(3)]
        edges = EdgesOfTriangle(second[0], second[1], second[2])
        for i in range(3):       
            for j in range(3):
                print(i, j)
                dist[i][j] = edges[i][0]*first[j][0] + edges[i][1]*first[j][1] + edges[i][2]
        if (DiffSigns(dist[0], eps) or DiffSigns(dist[1], eps) or DiffSigns(dist[2], eps)):
            return True
        else:
             return False



def Coplanar(first_3D, second_3D, N1, N2, eps):
    n_x = abs(N1[0])
    n_y = abs(N1[1])
    n_z = abs(N1[2])
    first_0 = [0.0]*2
    second_0 = [.0]*2
    first_1 = [0.0]*2
    second_1 = [.0]*2
    first_2 = [0.0]*2
    second_2 = [.0]*2
    if (( n_x > n_z ) and ( n_x >= n_y )):
        first_0[0] = first_3D[0][2]
        second_0[0] = second_3D[0][2]
        first_0[1] = first_3D[0][1]
        second_0[1] = second_3D[0][1]

        first_1[0] = first_3D[1][2]
        second_1[0] = second_3D[1][2]
        first_1[1] = first_3D[1][1]
        second_1[1] = second_3D[1][1]

        first_2[0] = first_3D[2][2]
        second_2[0] = second_3D[2][2]
        first_2[1] = first_3D[2][1]
        second_2[1] = second_3D[2][1]
    elif (( n_y > n_z ) and ( n_y >= n_x )):
        first_0[0] = first_3D[0][0]
        second_0[0] = second_3D[0][0]
        first_0[1] = first_3D[0][2]
        second_0[1] = second_3D[0][2]

        first_1[0] = first_3D[1][0]
        second_1[0] = second_3D[1][0]
        first_1[1] = first_3D[1][2]
        second_1[1] = second_3D[1][2]

        first_2[0] = first_3D[2][0]
        second_2[0] = second_3D[2][0]
        first_2[1] = first_3D[2][2]
        second_2[1] = second_3D[2][2]
    else:
        first_0[0] = first_3D[0][0]
        second_0[0] = second_3D[0][0]
        first_0[1] = first_3D[0][1]
        second_0[1] = second_3D[0][1]

        first_1[0] = first_3D[1][0]
        second_1[0] = second_3D[1][0]
        first_1[1] = first_3D[1][1]
        second_1[1] = second_3D[1][1]

        first_2[0] = first_3D[2][0]
        second_2[0] = second_3D[2][0]
        first_2[1] = first_3D[2][1]
        second_2[1] = second_3D[2][1]

    first_2D = [first_0, first_1, first_2]
    second_2D = [second_0, second_1, second_2]
    return TrianglesIntersect2D(first_2D, second_2D, eps)

def TrianglesIntersect3D(first, second):
    eps = 0.000001
    N1 = vector_product(subtract(first[1], first[0]), subtract(first[2], first[0]))
    d1 = scalar_product(minus(N1), first[0])
    N2 = vector_product(subtract(second[1], second[0]), subtract(second[2], second[0]))
    d2 = scalar_product(minus(N2), second[0])
    first_plane = [N1, d1]
    second_plane = [N2, d2]
    d1to2 = [0.0]*3
    d2to1 = [0.0]*3
    for i in range(3):
        d1to2[i] = scalar_product(second_plane[0], first[i])+ second_plane[1]
    mayIntersect = True
    if (not zeroDistance(d1to2, eps) and eq_signs(d1to2)):
        return False
    for i in range(3):
        d2to1[i] = scalar_product(first_plane[0], second[i])+ first_plane[1]
    if (not zeroDistance(d2to1, eps) and eq_signs(d2to1)):
        return False
    if zeroDistance(d1to2, eps):
        #co-planar
        return Coplanar(first, second, N1, N2,eps)
    else:
        D = vector_product(first_plane[0], second_plane[0])
        rightform(first, d1to2)
        rightform(second, d2to1)
        p1 = [0.0]*3
        p2 = [0.0]*3
        Int1 = [0.0]*2
        Int2 = [0.0]*2
        #interval on the first triangle
        for i in range(3):
            p1[i] = scalar_product(D, first[i])
        Int1[0] = p1[0] + (p1[1]-p1[0])*(d1to2[0]/(d1to2[0]-d1to2[1]))
        Int1[1] = p1[2] + (p1[1]-p1[2])*(d1to2[2]/(d1to2[2]-d1to2[1]))
        #interval on the second triangle
        for i in range(3):
            p2[i] = scalar_product(D, second[i])
        Int2[0] = p2[0] + (p2[1]-p2[0])*(d2to1[0]/(d2to1[0]-d2to1[1]))
        Int2[1] = p2[2] + (p2[1]-p2[2])*(d2to1[2]/(d2to1[2]-d2to1[1]))
        RightOrder2D(Int1)
        RightOrder2D(Int2)
        if (IntIntersects(Int1, Int2)):
            return True
        else:
            return False

def DrawTriangle(first, second, ax):
    x1 = [first[0][0], first[1][0], first[2][0]]
    y1 = [first[0][1], first[1][1], first[2][1]]
    z1 = [first[0][2], first[1][2], first[2][2]]
    x2 = [second[0][0], second[1][0], second[2][0]]
    y2 = [second[0][1], second[1][1], second[2][1]]
    z2 = [second[0][2], second[1][2], second[2][2]]
    verts1 = [list(zip(x1, y1, z1))]
    verts2 = [list(zip(x2, y2, z2))]
    tri1 = Poly3DCollection(verts1)
    tri1.set_color('r')
    tri2 = Poly3DCollection(verts2)
    tri2.set_color('g')
    ax.add_collection3d(tri1)
    ax.add_collection3d(tri2)


second = [[[0, 0, 0], [0, 5, 0], [6, 5, 0]], [[-1, 0, 0], [0, 2, 0], [0, 0, 0]], [[-1, 0, 0], [0, -1, 0], [0, 0, 0]], [[0, 0, 0], [0, 2, 0], [1, 0, 0]], [[0, 0, 0], [0, 4, 0], [4, 0, 0]], [[0, 0, 0], [0, 4, 0], [4, 0, 0]], [[0, 0, 0], [0, 4, 0], [4, 0, 0]], [[0, 0, 0], [0, 4, 0], [4, 0, 0]], [[0, 0, 0], [0, 4, 0], [5, 0, 0]]]
first = [[[1, 4, 0], [2, 4, 0], [2, 3, 0]], [[0, 0, 0], [0, 3, 0], [5, 0, 0]], [[0, 0, 0], [5, 0, 0], [0, 4, 0]], [[0, -1, 0], [0, 3, 0], [7, -1, 0]], [[1, 2, 0], [1, 1, -3], [0.5, 2, -2]], [[4, 0, 0], [1, 1, -3], [0.5, 2, -2]], [[1, 2, 0], [2, 1, 0], [0.5, 2, -2]], [[-1, 2, 0], [0, 2, 2], [0, 0, -2]], [[1, 1, 2], [5, 6, -2], [3, -4, -1]]]
for i in range(9):
    fig = plt.figure()
    ax = Axes3D(fig)
    ax.set_xlim(-7, 7)
    ax.set_ylim(-7, 7)
    ax.set_zlim(-7, 7)
    DrawTriangle(first[i], second[i], ax)
    
    print(TrianglesIntersect3D(first[i], second[i]))
    plt.show()
    print(TrianglesIntersect3D(first[i], second[i]))
#print(EdgesOfTriangle([0, 0], [2, 0], [0, 2]))