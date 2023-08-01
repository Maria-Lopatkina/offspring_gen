def equation_trio_wo_mut(f1, f2, m1, m2, k1, k2):
    if k1 == k2:
        if f1 == f2:
            return k1, 1
        else:
            return k1, 2
    else:
        if k1 == m1 and k2 == m2:
            for i in f1, f2:
                if i not in [k1, k2]:
                    return k1, k2, 4
            return k1, k2, 3
        else:
            if f1 == f2:
                return f1, 1
            else:
                for i in k1, k2:
                    if i not in [m1, m2] and i in [f1, f2]:
                        return i, 2


print(equation_trio_wo_mut(12, 14, 12, 13, 12, 13))
