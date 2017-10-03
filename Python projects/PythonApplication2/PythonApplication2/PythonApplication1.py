def BinarySearch(d, element):
    l = -1               
    r = len(d)    
    while l < r - 1:
        lr = (l + r) / 2
        m = int(lr)  + (not lr.is_integer()) 
        if d[m] <= element:
            l = m
        else: 
            r = m        
    return r

def LeastIncreasingSequence(a):
   n = len(a)
   d = [0]*(n)
   length = 0
   d[0] = -float('Inf')
   for i in range(1, n):
       d[i] = float('Inf')
   for i in range(n):
       if (not isinstance(a[i], int)):
           raise Exception("List consists non-integer value!")
       j = BinarySearch(d, a[i])
       a[i]
       if (d[j - 1] < a[i] and a[i] < d[j]):
           d[j] = a[i]
           length = max(length, j)
   return length
   
a = [1,2,4,5,2,7, 7.09]
print(LeastIncreasingSequence(a))