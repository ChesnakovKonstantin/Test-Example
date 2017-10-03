def power(x, y):
    if (y == 0):
       return 1
    if (y == -1):
        return 1.0 / x
    if (y == 1):
        return x
    y2 = y/2
    if (y > 0):
        p = power(x, int(y2) + (not y2.is_integer()))
    else:
        p = power(x, int(y2) - (not y2.is_integer()))
    p *= p
    if y % 2:
        p *= x
    return p

print(power(2, -5))