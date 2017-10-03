def LCS(x, y):
   #if (len(x1) < len(y1)):
   #    x = y1
   #    y = x1
   #else:
   #    x = x1
   #    y = y1
   m = len(x)
   n = len(y)
   lcs = [[0 for j in range(n)] for i in range(m)]
   for j in range(n):
      lcs[0][j] = 0
      lcs[1][j] = 0
   for i in range(0,m):
       lcs[1][0] = 0
       for j in range(0,n):
           lcs[0][j] = lcs[1][j]
           if (x[i] == y[j]):
               lcs[1][j] = lcs[0][j - 1] + 1
           else:
               if (lcs[0][j] >= lcs[1][j - 1]):
                   lcs[1][j] = lcs[0][j]
               else:
                   lcs[1][j] = lcs[1][j - 1]
   return lcs[1][n-1]

a = [1, 4, 6, 2, 1, 9, 0]
b= [4, 66, 54 , 34, 2, 9, 0]
print(LCS(a,b))
