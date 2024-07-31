import math

def Wigner3j(j1, j2, j3, m1, m2, m3):
    # error checking
    if not (2*j1 == int(2*j1) and 2*j2 == int(2*j2) and 2*j3 == int(2*j3) and
            2*m1 == int(2*m1) and 2*m2 == int(2*m2) and 2*m3 == int(2*m3)):
        raise ValueError('All arguments must be integers or half-integers.')

    if not (j1 - m1 == int(j1 - m1)):
        raise ValueError('2*j1 and 2*m1 must have the same parity')

    if not (j2 - m2 == int(j2 - m2)):
        raise ValueError('2*j2 and 2*m2 must have the same parity')

    if not (j3 - m3 == int(j3 - m3)):
        raise ValueError('2*j3 and 2*m3 must have the same parity')

    if j3 > j1 + j2 or j3 < abs(j1 - j2):
        raise ValueError('j3 is out of bounds.')

    if abs(m1) > j1:
        raise ValueError('m1 is out of bounds.')

    if abs(m2) > j2:
        raise ValueError('m2 is out of bounds.')

    if abs(m3) > j3:
        raise ValueError('m3 is out of bounds.')

    t1 = j2 - m1 - j3
    t2 = j1 + m2 - j3
    t3 = j1 + j2 - j3
    t4 = j1 - m1
    t5 = j2 + m2

    tmin = max(0, max(t1, t2))
    tmax = min(t3, min(t4, t5))

    wigner = 0

    for t in range(tmin, tmax + 1):
        wigner += (-1)**t / (math.factorial(t) * math.factorial(t-t1) * math.factorial(t-t2) *
                             math.factorial(t3-t) * math.factorial(t4-t) * math.factorial(t5-t))

    wigner *= (-1)**(j1 - j2 - m3)
    wigner *= math.sqrt(math.factorial(j1 + j2 - j3) * math.factorial(j1 - j2 + j3) *
                        math.factorial(-j1 + j2 + j3) / math.factorial(j1 + j2 + j3 + 1) *
                        math.factorial(j1 + m1) * math.factorial(j1 - m1) *
                        math.factorial(j2 + m2) * math.factorial(j2 - m2) *
                        math.factorial(j3 + m3) * math.factorial(j3 - m3))

    return wigner