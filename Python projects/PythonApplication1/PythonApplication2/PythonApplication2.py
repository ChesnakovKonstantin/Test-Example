import sys

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
    res = a[0]*b[0] + a[1]*b[1]+a[2]*b[2]
    return res

def eq_signs(a):
    if (a[0]<0):
        for el in a:
            if (el>0):
                return False
    elif (a[0] > 0):
        for el in a:
            if (el<0):
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

def swap(a, b, i, j):
    tempa = a[i]
    tempb = b[i]
    a[i] = a[j]
    b[i] = b[j]
    a[j] = tempa
    b[j] = tempb

#rightform: the first and the last vertice make d of the same sign to another plane
def rightform(a, d):
    if (sign(d[0]) != sign(d[2])):
        if (sign(d[0]) == sign(d[1])):
            swap(a, d, 1, 2)
        else:
            swap(a, d, 0, 1)

def IntIntersects(T1, T2):
    if (T2[0] <= T1[1] and T0[0] <= T1[1]):
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
    #return (b1 == b2) and (b2 == b3)
    return ((b1 and b2) or (not b1 and not b2)) and ((b3 and b2) or (not b3 and not b2))

def TrianglesIntersect2D(first, second, eps):
    if (PointInTriangle2D(first[0], second[0], second[1], second[2], eps) or PointInTriangle2D(first[1], second[0], second[1], second[2], eps) or PointInTriangle2D(first[2], second[0], second[1], second[2], eps) or PointInTriangle2D(second[0], first[0], first[1], first[2], eps) or PointInTriangle2D(second[1], first[0], first[1], first[2], eps) or PointInTriangle2D(second[2], first[0], first[1], first[2], eps)):
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
        if (IntIntersects(Int1, Int2)):
            return True
        else:
            return False

#first = [[0, 1.5, 2.71828], [1.17, 9.01, 9.08], [4.9, 0.4, 3.4]]
first = [[-1, 3.6, 3.1415], [9.02, 7.15, 1.7], [6.08, 2.99, 1.56]]
second = [[-1, 3.6, 3.1415], [9.02, 7.15, 1.7], [6.08, 2.99, 1.56]]
print(TrianglesIntersect3D(first, second))