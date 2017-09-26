import sys
print(sys.platform)
K = 2048
mylist = [66, 333, 333, 1, 1234, 1234, 78, 78, 9, 78]
numlist = [0]*K
nDifferent = 0
for i in mylist:
    numlist[i]=1;
for i in range(0, len(numlist)):
    nDifferent += numlist[i]
print(nDifferent)
