import random


def inx_mutation(n_p, k, m_r):
    count = 0
    list_of_mut = []
    mutation = int(n_p * m_r * k * 22)
    for m in range(mutation):
        a = random.choice(range(0, n_p * k))
        if a not in list_of_mut:
            list_of_mut.append(a)
            count += 1
    return list_of_mut, mutation

number_of_pairs = 10
kids = 3
mut_rate = 0.0045
print(inx_mutation(number_of_pairs, kids, mut_rate))
