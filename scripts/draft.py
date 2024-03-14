import copy
import pandas as pd
import random
import time
import yaml


def empty_freq_table(path):  # Create empty dictionaries for calculation of allele frequencies
    freq_df = pd.read_excel(path)
    f_columns = freq_df.columns.values.tolist()
    key_dict = []
    for el in range(1, len(f_columns) - 2, 2):
        key_dict.append(list(f_columns[el].split("_"))[0])
    freq_dict = dict.fromkeys(key_dict)
    return freq_dict


def calculate_frequencies(p, d, path):  # Calculate alleles frequencies
    freq_df = pd.read_excel(path)
    f_columns = freq_df.columns.values.tolist()
    temp_freq_df = freq_df.loc[(freq_df['population'] == p)]
    for al in range(1, len(f_columns) - 2, 2):
        allele_name = f_columns[al].split('_')[0]
        allele_df_1 = temp_freq_df.iloc[:, [al]]
        allele_df_1 = allele_df_1.set_axis([allele_name], axis=1)
        allele_df_2 = temp_freq_df.iloc[:, [al + 1]]
        allele_df_2 = allele_df_2.set_axis([allele_name], axis=1)
        allele_df = pd.concat([allele_df_1, allele_df_2], axis=0)
        a = allele_df.value_counts() / len(allele_df.index)
        aa = a.to_dict()
        keys_list = list(aa.keys())
        new_keys_list = []
        for h in keys_list:
            new_keys_list.append(h[0])
        new_aa = d.fromkeys(new_keys_list)
        for h in range(0, len(new_keys_list)):
            new_aa[new_keys_list[h]] = aa[keys_list[h]]
        d[allele_name] = new_aa
    k, v = "pmin", 0.0003
    for ii in d.keys():
        d[ii][k] = v
    return d


# functions for LR calculation
def hyp_sibl(s1, s2, k1, k2, our_dict, rus_dict, str_elem):
    #  1 - for our frequencies, 2 - for rus
    print()
    print(str_elem, ":")
    print("k1=", k1, "k2=", k2, "s1=", s1, "s2=", s2)
    print("p1, p2:")
    p1, p2 = calculate_p(s1, s2, k1, k2, our_dict, str_elem)
    print("p3, p4:")
    p3, p4 = calculate_p(s1, s2, k1, k2, rus_dict, str_elem)
    return p1, p2, p3, p4


def calculate_p(s1, s2, k1, k2, dic, allele):
    p1 = p2 = None
    if k1 == k2:
        if s1 == s2 and s1 == k1:
            if k1 in dic[allele]:
                p1 = 1
                p2 = (dic[allele][k1] * (2 - dic[allele][k1])) ** 2
            else:
                p1 = 1
                p2 = (dic[allele]["pmin"] * (2 - dic[allele]["pmin"])) ** 2
            print("s1 == s2 and s1 == k1", allele, p1, p2)
        elif k1 == s1 or k1 == s2:
            if k1 == s1:
                a = s1
                b = s2
            elif k1 == s2:
                a = s2
                b = s1
            if a in dic[allele] and b in dic[allele]:
                p1 = (2 * dic[allele][a] * (2 - dic[allele][a]) * 2 * dic[allele][a] * dic[allele][b] -
                      2 * dic[allele][a] * dic[allele][b] * 2 * dic[allele][a] * dic[allele][b]) / \
                     (2 * dic[allele][a] * (2 - dic[allele][a]) * dic[allele][b] * (2 - dic[allele][b]) -
                      2 * dic[allele][a] * dic[allele][b] * 2 * dic[allele][a] * dic[allele][b])
                p2 = (dic[allele][a] * (2 - dic[allele][a])) ** 2
            elif a not in dic[allele] and b in dic[allele]:
                p1 = (2 * dic[allele]["pmin"] * (2 - dic[allele]["pmin"]) * 2 * dic[allele]["pmin"] * dic[allele][b] -
                      2 * dic[allele]["pmin"] * dic[allele][b] * 2 * dic[allele]["pmin"] * dic[allele][b]) / \
                     (2 * dic[allele]["pmin"] * (2 - dic[allele]["pmin"]) * dic[allele][b] * (2 - dic[allele][b]) -
                      2 * dic[allele]["pmin"] * dic[allele][b] * 2 * dic[allele]["pmin"] * dic[allele][b])
                p2 = (dic[allele]["pmin"] * (2 - dic[allele]["pmin"])) ** 2
            elif a in dic[allele] and b not in dic[allele]:
                p1 = (2 * dic[allele][a] * (2 - dic[allele][a]) * 2 * dic[allele][a] * dic[allele]["pmin"] -
                      2 * dic[allele][a] * dic[allele]["pmin"] * 2 * dic[allele][a] * dic[allele]["pmin"]) / \
                     (2 * dic[allele][a] * (2 - dic[allele][a]) * dic[allele]["pmin"] * (2 - dic[allele]["pmin"]) -
                      2 * dic[allele][a] * dic[allele]["pmin"] * 2 * dic[allele][a] * dic[allele]["pmin"])
                p2 = (dic[allele][a] * (2 - dic[allele][a])) ** 2
            elif a not in dic[allele] and a not in dic[allele]:
                p1 = (2 * dic[allele]["pmin"] * (2 - dic[allele]["pmin"]) * 2 * dic[allele]["pmin"] * dic[allele]["pmin"] -
                      2 * dic[allele]["pmin"] * dic[allele]["pmin"] * 2 * dic[allele]["pmin"] * dic[allele]["pmin"]) / \
                     (2 * dic[allele]["pmin"] * (2 - dic[allele]["pmin"]) * dic[allele]["pmin"] * (2 - dic[allele]["pmin"]) -
                      2 * dic[allele]["pmin"] * dic[allele]["pmin"] * 2 * dic[allele]["pmin"] * dic[allele]["pmin"])
                p2 = (dic[allele]["pmin"] * (2 - dic[allele]["pmin"])) ** 2
            print("k1 == s1 or k1 == s2", allele, p1, p2)
        elif s1 == s2 and k1 != s1:
            a = k1
            b = s1
            if a in dic[allele] and b in dic[allele]:
                p1 = (2 * dic[allele][a] * dic[allele][b] * 2 * dic[allele][a] * dic[allele][b]) / \
                     (dic[allele][b] * (2 - dic[allele][b])) ** 2
                p2 = (dic[allele][a] * (2 - dic[allele][a])) ** 2
            elif a not in dic[allele] and b in dic[allele]:
                p1 = (2 * dic[allele]["pmin"] * dic[allele][b] * 2 * dic[allele]["pmin"] * dic[allele][b]) / \
                     (dic[allele][b] * (2 - dic[allele][b])) ** 2
                p2 = (dic[allele]["pmin"] * (2 - dic[allele]["pmin"])) ** 2
            elif a in dic[allele] and b not in dic[allele]:
                p1 = (2 * dic[allele][a] * dic[allele]["pmin"] * 2 * dic[allele][a] * dic[allele]["pmin"]) / \
                     (dic[allele]["pmin"] * (2 - dic[allele]["pmin"])) ** 2
                p2 = (dic[allele][a] * (2 - dic[allele][a])) ** 2
            elif a not in dic[allele] and a not in dic[allele]:
                p1 = (2 * dic[allele]["pmin"] * dic[allele]["pmin"] * 2 * dic[allele]["pmin"] * dic[allele]["pmin"]) / \
                     (dic[allele]["pmin"] * (2 - dic[allele]["pmin"])) ** 2
                p2 = (dic[allele]["pmin"] * (2 - dic[allele]["pmin"])) ** 2
            print("s1 == s2 and k1 != s1", allele, p1, p2)
        elif k1 != s1 and k1 != s2 and s1 != s2:
            a = k1
            b = s1
            c = s2
            if a in dic[allele] and b in dic[allele] and c in dic[allele]:
                p1 = (2 * 2 * dic[allele][a] * dic[allele][b] * 2 * dic[allele][a] * dic[allele][c]) / \
                     (2 * dic[allele][b] * (2 - dic[allele][b]) * dic[allele][c] * (2 - dic[allele][c]) -
                      2 * dic[allele][b] * dic[allele][c] * 2 * dic[allele][b] * dic[allele][c])
                p2 = (dic[allele][a] * (2 - dic[allele][a])) ** 2
            elif a not in dic[allele] and b in dic[allele] and c in dic[allele]:
                p1 = (2 * 2 * dic[allele]["pmin"] * dic[allele][b] * 2 * dic[allele]["pmin"] * dic[allele][c]) / \
                     (2 * dic[allele][b] * (2 - dic[allele][b]) * dic[allele][c] * (2 - dic[allele][c]) -
                      2 * dic[allele][b] * dic[allele][c] * 2 * dic[allele][b] * dic[allele][c])
                p2 = (dic[allele]["pmin"] * (2 - dic[allele]["pmin"])) ** 2
            elif a in dic[allele] and b not in dic[allele] and c in dic[allele]:
                p1 = (2 * 2 * dic[allele][a] * dic[allele]["pmin"] * 2 * dic[allele][a] * dic[allele][c]) / \
                     (2 * dic[allele]["pmin"] * (2 - dic[allele]["pmin"]) * dic[allele][c] * (2 - dic[allele][c]) -
                      2 * dic[allele]["pmin"] * dic[allele][c] * 2 * dic[allele]["pmin"] * dic[allele][c])
                p2 = (dic[allele][a] * (2 - dic[allele][a])) ** 2
            elif a in dic[allele] and b in dic[allele] and c not in dic[allele]:
                p1 = (2 * 2 * dic[allele][a] * dic[allele][b] * 2 * dic[allele][a] * dic[allele]["pmin"]) / \
                     (2 * dic[allele][b] * (2 - dic[allele][b]) * dic[allele]["pmin"] * (2 - dic[allele]["pmin"]) -
                      2 * dic[allele][b] * dic[allele]["pmin"] * 2 * dic[allele][b] * dic[allele]["pmin"])
                p2 = (dic[allele][a] * (2 - dic[allele][a])) ** 2
            elif a not in dic[allele] and b not in dic[allele] and c in dic[allele]:
                p1 = (2 * 2 * dic[allele]["pmin"] * dic[allele]["pmin"] * 2 * dic[allele]["pmin"] * dic[allele][c]) / \
                     (2 * dic[allele]["pmin"] * (2 - dic[allele]["pmin"]) * dic[allele][c] * (2 - dic[allele][c]) -
                      2 * dic[allele]["pmin"] * dic[allele][c] * 2 * dic[allele]["pmin"] * dic[allele][c])
                p2 = (dic[allele]["pmin"] * (2 - dic[allele]["pmin"])) ** 2
            elif a not in dic[allele] and b in dic[allele] and c not in dic[allele]:
                p1 = (2 * 2 * dic[allele]["pmin"] * dic[allele][b] * 2 * dic[allele]["pmin"] * dic[allele]["pmin"]) / \
                     (2 * dic[allele][b] * (2 - dic[allele][b]) * dic[allele]["pmin"] * (2 - dic[allele]["pmin"]) -
                      2 * dic[allele][b] * dic[allele]["pmin"] * 2 * dic[allele][b] * dic[allele]["pmin"])
                p2 = (dic[allele]["pmin"] * (2 - dic[allele]["pmin"])) ** 2
            elif a in dic[allele] and b not in dic[allele] and c not in dic[allele]:
                p1 = (2 * 2 * dic[allele][a] * dic[allele]["pmin"] * 2 * dic[allele][a] * dic[allele]["pmin"]) / \
                     (2 * dic[allele]["pmin"] * (2 - dic[allele]["pmin"]) * dic[allele]["pmin"] * (2 - dic[allele]["pmin"]) -
                      2 * dic[allele]["pmin"] * dic[allele]["pmin"] * 2 * dic[allele]["pmin"] * dic[allele]["pmin"])
                p2 = (dic[allele][a] * (2 - dic[allele][a])) ** 2
            elif a not in dic[allele] and b not in dic[allele] and c not in dic[allele]:
                p1 = (2 * 2 * dic[allele]["pmin"] * dic[allele]["pmin"] * 2 * dic[allele]["pmin"] * dic[allele]["pmin"]) / \
                     (2 * dic[allele]["pmin"] * (2 - dic[allele]["pmin"]) * dic[allele]["pmin"] * (2 - dic[allele]["pmin"]) -
                      2 * dic[allele]["pmin"] * dic[allele]["pmin"] * 2 * dic[allele]["pmin"] * dic[allele]["pmin"])
                p2 = (dic[allele]["pmin"] * (2 - dic[allele]["pmin"])) ** 2
            print("k1 != s1 and k1 != s2 and s1 != s2", allele, p1, p2)
    else:
        if s1 == s2 and (s1 == k1 or s1 == k2):
            a = s1
            b = None
            if s1 == k1:
                b = k2
            elif s1 == k2:
                b = k1
            if a in dic[allele] and b in dic[allele]:
                p1 = (2 * dic[allele][a] * (2 - dic[allele][a]) * 2 * dic[allele][a] * dic[allele][b] -
                      2 * dic[allele][a] * dic[allele][b] * 2 * dic[allele][a] * dic[allele][b]) / \
                     (dic[allele][a] * (2 - dic[allele][a])) ** 2
                p2 = 2 * dic[allele][a] * (2 - dic[allele][a]) * dic[allele][b] * (2 - dic[allele][b]) - \
                     2 * dic[allele][a] * dic[allele][b] * 2 * dic[allele][a] * dic[allele][b]
            if a not in dic[allele] and b in dic[allele]:
                p1 = (2 * dic[allele]["pmin"] * (2 - dic[allele]["pmin"]) * 2 * dic[allele]["pmin"] * dic[allele][b] -
                      2 * dic[allele]["pmin"] * dic[allele][b] * 2 * dic[allele]["pmin"] * dic[allele][b]) / \
                     (dic[allele]["pmin"] * (2 - dic[allele]["pmin"])) ** 2
                p2 = 2 * dic[allele]["pmin"] * (2 - dic[allele]["pmin"]) * dic[allele][b] * (2 - dic[allele][b]) - \
                     2 * dic[allele]["pmin"] * dic[allele][b] * 2 * dic[allele]["pmin"] * dic[allele][b]
            if a in dic[allele] and b not in dic[allele]:
                p1 = (2 * dic[allele][a] * (2 - dic[allele][a]) * 2 * dic[allele][a] * dic[allele]["pmin"] -
                      2 * dic[allele][a] * dic[allele]["pmin"] * 2 * dic[allele][a] * dic[allele]["pmin"]) / \
                     (dic[allele][a] * (2 - dic[allele][a])) ** 2
                p2 = 2 * dic[allele][a] * (2 - dic[allele][a]) * dic[allele]["pmin"] * (2 - dic[allele]["pmin"]) - \
                     2 * dic[allele][a] * dic[allele]["pmin"] * 2 * dic[allele][a] * dic[allele]["pmin"]
            if a not in dic[allele] and b not in dic[allele]:
                p1 = (2 * dic[allele]["pmin"] * (2 - dic[allele]["pmin"]) * 2 * dic[allele]["pmin"] * dic[allele]["pmin"] -
                      2 * dic[allele]["pmin"] * dic[allele]["pmin"] * 2 * dic[allele]["pmin"] * dic[allele]["pmin"]) / \
                     (dic[allele]["pmin"] * (2 - dic[allele]["pmin"])) ** 2
                p2 = 2 * dic[allele]["pmin"] * (2 - dic[allele]["pmin"]) * dic[allele]["pmin"] * (2 - dic[allele]["pmin"]) - \
                     2 * dic[allele]["pmin"] * dic[allele]["pmin"] * 2 * dic[allele]["pmin"] * dic[allele]["pmin"]
            print("s1 == s2 and (s1 == k1 or s1 == k2)", allele, p1, p2)
        if k1 == s1 and k2 == s2:
            a = k1
            b = k2
            if a in dic[allele] and b in dic[allele]:
                p1 = 1
                p2 = 2 * dic[allele][a] * (2 - dic[allele][a]) * dic[allele][b] * (2 - dic[allele][b]) - \
                     2 * dic[allele][a] * dic[allele][b] * 2 * dic[allele][a] * dic[allele][b]
            if a not in dic[allele] and b in dic[allele]:
                p1 = 1
                p2 = 2 * dic[allele]["pmin"] * (2 - dic[allele]["pmin"]) * dic[allele][b] * (2 - dic[allele][b]) - \
                     2 * dic[allele]["pmin"] * dic[allele][b] * 2 * dic[allele]["pmin"] * dic[allele][b]
            if a in dic[allele] and b not in dic[allele]:
                p1 = 1
                p2 = 2 * dic[allele][a] * (2 - dic[allele][a]) * dic[allele]["pmin"] * (2 - dic[allele]["pmin"]) - \
                     2 * dic[allele][a] * dic[allele]["pmin"] * 2 * dic[allele][a] * dic[allele]["pmin"]
            if a not in dic[allele] and b not in dic[allele]:
                p1 = 1
                p2 = 2 * dic[allele]["pmin"] * (2 - dic[allele]["pmin"]) * dic[allele]["pmin"] * (2 - dic[allele]["pmin"]) - \
                     2 * dic[allele]["pmin"] * dic[allele]["pmin"] * 2 * dic[allele]["pmin"] * dic[allele]["pmin"]
            print("k1 == s1 and k2 == s2", allele, p1, p2)
        if s1 != s2 and (k1 == s1 and k2 != s2 or k1 == s2 and k2 != s1 or k1 != s1 and k2 == s2 or k1 != s2 and k2 == s1):
            a = None
            b = None
            c = None
            if k1 == s1 and k2 != s2:
                a = k1
                b = k2
                c = s2
            if k1 == s2 and k2 != s1:
                a = k1
                b = k2
                c = s1
            if k1 != s1 and k2 == s2:
                a = k2
                b = k1
                c = s1
            if k1 != s2 and k2 == s1:
                a = k2
                b = k1
                c = s2
            if a in dic[allele] and b in dic[allele] and c in dic[allele]:
                p1 = (2 * dic[allele][a] * (2 - dic[allele][a]) * 2 * dic[allele][b] * dic[allele][c] +
                      2 * 2 * dic[allele][a] * dic[allele][b] * 2 * dic[allele][a] * dic[allele][c]) / \
                     (2 * dic[allele][a] * (2 - dic[allele][a]) * dic[allele][c] * (2 - dic[allele][c]) -
                      2 * dic[allele][a] * dic[allele][c] * 2 * dic[allele][a] * dic[allele][c])
                p2 = 2 * dic[allele][a] * (2 - dic[allele][a]) * dic[allele][b] * (2 - dic[allele][b]) - \
                     2 * dic[allele][a] * dic[allele][b] * 2 * dic[allele][a] * dic[allele][b]
            elif a not in dic[allele] and b in dic[allele] and c in dic[allele]:
                p1 = (2 * dic[allele]["pmin"] * (2 - dic[allele]["pmin"]) * 2 * dic[allele][b] * dic[allele][c] +
                      2 * 2 * dic[allele]["pmin"] * dic[allele][b] * 2 * dic[allele]["pmin"] * dic[allele][c]) / \
                     (2 * dic[allele]["pmin"] * (2 - dic[allele]["pmin"]) * dic[allele][c] * (2 - dic[allele][c]) -
                      2 * dic[allele]["pmin"] * dic[allele][c] * 2 * dic[allele]["pmin"] * dic[allele][c])
                p2 = 2 * dic[allele]["pmin"] * (2 - dic[allele]["pmin"]) * dic[allele][b] * (2 - dic[allele][b]) - \
                     2 * dic[allele]["pmin"] * dic[allele][b] * 2 * dic[allele]["pmin"] * dic[allele][b]
            elif a in dic[allele] and b not in dic[allele] and c in dic[allele]:
                p1 = (2 * dic[allele][a] * (2 - dic[allele][a]) * 2 * dic[allele]["pmin"] * dic[allele][c] +
                      2 * 2 * dic[allele][a] * dic[allele]["pmin"] * 2 * dic[allele][a] * dic[allele][c]) / \
                     (2 * dic[allele][a] * (2 - dic[allele][a]) * dic[allele][c] * (2 - dic[allele][c]) -
                      2 * dic[allele][a] * dic[allele][c] * 2 * dic[allele][a] * dic[allele][c])
                p2 = 2 * dic[allele][a] * (2 - dic[allele][a]) * dic[allele]["pmin"] * (2 - dic[allele]["pmin"]) - \
                     2 * dic[allele][a] * dic[allele]["pmin"] * 2 * dic[allele][a] * dic[allele]["pmin"]
            elif a in dic[allele] and b in dic[allele] and c not in dic[allele]:
                p1 = (2 * dic[allele][a] * (2 - dic[allele][a]) * 2 * dic[allele][b] * dic[allele]["pmin"] +
                      2 * 2 * dic[allele][a] * dic[allele][b] * 2 * dic[allele][a] * dic[allele]["pmin"]) / \
                     (2 * dic[allele][a] * (2 - dic[allele][a]) * dic[allele]["pmin"] * (2 - dic[allele]["pmin"]) -
                      2 * dic[allele][a] * dic[allele]["pmin"] * 2 * dic[allele][a] * dic[allele]["pmin"])
                p2 = 2 * dic[allele][a] * (2 - dic[allele][a]) * dic[allele][b] * (2 - dic[allele][b]) - \
                     2 * dic[allele][a] * dic[allele][b] * 2 * dic[allele][a] * dic[allele][b]
            elif a not in dic[allele] and b not in dic[allele] and c in dic[allele]:
                p1 = (2 * dic[allele]["pmin"] * (2 - dic[allele]["pmin"]) * 2 * dic[allele]["pmin"] * dic[allele][c] +
                      2 * 2 * dic[allele]["pmin"] * dic[allele]["pmin"] * 2 * dic[allele]["pmin"] * dic[allele][c]) / \
                     (2 * dic[allele]["pmin"] * (2 - dic[allele]["pmin"]) * dic[allele][c] * (2 - dic[allele][c]) -
                      2 * dic[allele]["pmin"] * dic[allele][c] * 2 * dic[allele]["pmin"] * dic[allele][c])
                p2 = 2 * dic[allele]["pmin"] * (2 - dic[allele]["pmin"]) * dic[allele]["pmin"] * (2 - dic[allele]["pmin"]) - \
                     2 * dic[allele]["pmin"] * dic[allele]["pmin"] * 2 * dic[allele]["pmin"] * dic[allele]["pmin"]
            elif a not in dic[allele] and b in dic[allele] and c not in dic[allele]:
                p1 = (2 * dic[allele]["pmin"] * (2 - dic[allele]["pmin"]) * 2 * dic[allele][b] * dic[allele]["pmin"] +
                      2 * 2 * dic[allele]["pmin"] * dic[allele][b] * 2 * dic[allele]["pmin"] * dic[allele]["pmin"]) / \
                     (2 * dic[allele]["pmin"] * (2 - dic[allele]["pmin"]) * dic[allele]["pmin"] * (2 - dic[allele]["pmin"]) -
                      2 * dic[allele]["pmin"] * dic[allele]["pmin"] * 2 * dic[allele]["pmin"] * dic[allele]["pmin"])
                p2 = 2 * dic[allele]["pmin"] * (2 - dic[allele]["pmin"]) * dic[allele][b] * (2 - dic[allele][b]) - \
                     2 * dic[allele]["pmin"] * dic[allele][b] * 2 * dic[allele]["pmin"] * dic[allele][b]
            elif a in dic[allele] and b not in dic[allele] and c not in dic[allele]:
                p1 = (2 * dic[allele][a] * (2 - dic[allele][a]) * 2 * dic[allele]["pmin"] * dic[allele]["pmin"] +
                      2 * 2 * dic[allele][a] * dic[allele]["pmin"] * 2 * dic[allele][a] * dic[allele]["pmin"]) / \
                     (2 * dic[allele][a] * (2 - dic[allele][a]) * dic[allele]["pmin"] * (2 - dic[allele]["pmin"]) -
                      2 * dic[allele][a] * dic[allele]["pmin"] * 2 * dic[allele][a] * dic[allele]["pmin"])
                p2 = 2 * dic[allele][a] * (2 - dic[allele][a]) * dic[allele]["pmin"] * (2 - dic[allele]["pmin"]) - \
                     2 * dic[allele][a] * dic[allele]["pmin"] * 2 * dic[allele][a] * dic[allele]["pmin"]
            elif a not in dic[allele] and b not in dic[allele] and c not in dic[allele]:
                p1 = (2 * dic[allele]["pmin"] * (2 - dic[allele]["pmin"]) * 2 * dic[allele]["pmin"] * dic[allele]["pmin"] +
                      2 * 2 * dic[allele]["pmin"] * dic[allele]["pmin"] * 2 * dic[allele]["pmin"] * dic[allele]["pmin"]) / \
                     (2 * dic[allele]["pmin"] * (2 - dic[allele]["pmin"]) * dic[allele]["pmin"] * (2 - dic[allele]["pmin"]) -
                      2 * dic[allele]["pmin"] * dic[allele]["pmin"] * 2 * dic[allele]["pmin"] * dic[allele]["pmin"])
                p2 = 2 * dic[allele]["pmin"] * (2 - dic[allele]["pmin"]) * dic[allele]["pmin"] * (2 - dic[allele]["pmin"]) - \
                     2 * dic[allele]["pmin"] * dic[allele]["pmin"] * 2 * dic[allele]["pmin"] * dic[allele]["pmin"]
            print("s1 != s2 and (k1 == s1 and k2 != k2 or k1 == s2 and k2 != s1 or k1 != s1 and k2 == s2 or k1 != s2 and k2 == s1)", allele, p1, p2)
        if s1 == s2 and k1 != s1 and k2 != s1:
            a = k1
            b = k2
            c = s1
            if a in dic[allele] and b in dic[allele] and c in dic[allele]:
                p1 = (2 * 2 * dic[allele][a] * dic[allele][c] * 2 * dic[allele][b] * dic[allele][c]) / \
                     (dic[allele][c] * (2 - dic[allele][c])) ** 2
                p2 = 2 * dic[allele][a] * (2 - dic[allele][a]) * dic[allele][b] * (2 - dic[allele][b]) - \
                     2 * dic[allele][a] * dic[allele][b] * 2 * dic[allele][a] * dic[allele][b]
            elif a not in dic[allele] and b in dic[allele] and c in dic[allele]:
                p1 = (2 * 2 * dic[allele]["pmin"] * dic[allele][c] * 2 * dic[allele][b] * dic[allele][c]) / \
                     (dic[allele][c] * (2 - dic[allele][c])) ** 2
                p2 = 2 * dic[allele]["pmin"] * (2 - dic[allele]["pmin"]) * dic[allele][b] * (2 - dic[allele][b]) - \
                     2 * dic[allele]["pmin"] * dic[allele][b] * 2 * dic[allele]["pmin"] * dic[allele][b]
            elif a in dic[allele] and b not in dic[allele] and c in dic[allele]:
                p1 = (2 * 2 * dic[allele][a] * dic[allele][c] * 2 * dic[allele]["pmin"] * dic[allele][c]) / \
                     (dic[allele][c] * (2 - dic[allele][c])) ** 2
                p2 = 2 * dic[allele][a] * (2 - dic[allele][a]) * dic[allele]["pmin"] * (2 - dic[allele]["pmin"]) - \
                     2 * dic[allele][a] * dic[allele]["pmin"] * 2 * dic[allele][a] * dic[allele]["pmin"]
            elif a in dic[allele] and b in dic[allele] and c not in dic[allele]:
                p1 = (2 * 2 * dic[allele][a] * dic[allele]["pmin"] * 2 * dic[allele][b] * dic[allele]["pmin"]) / \
                     (dic[allele]["pmin"] * (2 - dic[allele]["pmin"])) ** 2
                p2 = 2 * dic[allele][a] * (2 - dic[allele][a]) * dic[allele][b] * (2 - dic[allele][b]) - \
                     2 * dic[allele][a] * dic[allele][b] * 2 * dic[allele][a] * dic[allele][b]
            elif a not in dic[allele] and b not in dic[allele] and c in dic[allele]:
                p1 = (2 * 2 * dic[allele]["pmin"] * dic[allele][c] * 2 * dic[allele]["pmin"] * dic[allele][c]) / \
                     (dic[allele][c] * (2 - dic[allele][c])) ** 2
                p2 = 2 * dic[allele]["pmin"] * (2 - dic[allele]["pmin"]) * dic[allele]["pmin"] * (2 - dic[allele]["pmin"]) - \
                     2 * dic[allele]["pmin"] * dic[allele]["pmin"] * 2 * dic[allele]["pmin"] * dic[allele]["pmin"]
            elif a not in dic[allele] and b in dic[allele] and c not in dic[allele]:
                p1 = (2 * 2 * dic[allele]["pmin"] * dic[allele]["pmin"] * 2 * dic[allele][b] * dic[allele]["pmin"]) / \
                     (dic[allele]["pmin"] * (2 - dic[allele]["pmin"])) ** 2
                p2 = 2 * dic[allele]["pmin"] * (2 - dic[allele]["pmin"]) * dic[allele][b] * (2 - dic[allele][b]) - \
                     2 * dic[allele]["pmin"] * dic[allele][b] * 2 * dic[allele]["pmin"] * dic[allele][b]
            elif a in dic[allele] and b not in dic[allele] and c not in dic[allele]:
                p1 = (2 * 2 * dic[allele][a] * dic[allele]["pmin"] * 2 * dic[allele]["pmin"] * dic[allele]["pmin"]) / \
                     (dic[allele]["pmin"] * (2 - dic[allele]["pmin"])) ** 2
                p2 = 2 * dic[allele][a] * (2 - dic[allele][a]) * dic[allele]["pmin"] * (2 - dic[allele]["pmin"]) - \
                     2 * dic[allele][a] * dic[allele]["pmin"] * 2 * dic[allele][a] * dic[allele]["pmin"]
            elif a not in dic[allele] and b not in dic[allele] and c not in dic[allele]:
                p1 = (2 * 2 * dic[allele]["pmin"] * dic[allele]["pmin"] * 2 * dic[allele]["pmin"] * dic[allele]["pmin"]) / \
                     (dic[allele]["pmin"] * (2 - dic[allele]["pmin"])) ** 2
                p2 = 2 * dic[allele]["pmin"] * (2 - dic[allele]["pmin"]) * dic[allele]["pmin"] * (2 - dic[allele]["pmin"]) - \
                     2 * dic[allele]["pmin"] * dic[allele]["pmin"] * 2 * dic[allele]["pmin"] * dic[allele]["pmin"]
            print("s1 == s2 and k1 != s1 and k2 != s1", allele, p1, p2)
        if s1 != s2 and k1 != s1 and k1 != s2 and k2 != s1 and k2 != s2:
            a = k1
            b = k2
            c = s1
            d = s2
            if a in dic[allele] and b in dic[allele] and c in dic[allele] and d in dic[allele]:
                p1 = (2 * 2 * dic[allele][a] * dic[allele][c] * 2 * dic[allele][b] * dic[allele][d] + 2 *
                      2 * dic[allele][a] * dic[allele][d] * 2 * dic[allele][b] * dic[allele][c]) / \
                     (2 * dic[allele][c] * (2 - dic[allele][c]) * dic[allele][d] * (2 - dic[allele][d]) -
                      2 * dic[allele][c] * dic[allele][d] * 2 * dic[allele][c] * dic[allele][d])
                p2 = 2 * dic[allele][a] * (2 - dic[allele][a]) * dic[allele][b] * (2 - dic[allele][b]) - \
                     2 * dic[allele][a] * dic[allele][b] * 2 * dic[allele][a] * dic[allele][b]
            elif a not in dic[allele] and b in dic[allele] and c in dic[allele] and d in dic[allele]:
                p1 = (2 * 2 * dic[allele]["pmin"] * dic[allele][c] * 2 * dic[allele][b] * dic[allele][d] + 2 *
                      2 * dic[allele]["pmin"] * dic[allele][d] * 2 * dic[allele][b] * dic[allele][c]) / \
                     (2 * dic[allele][c] * (2 - dic[allele][c]) * dic[allele][d] * (2 - dic[allele][d]) -
                      2 * dic[allele][c] * dic[allele][d] * 2 * dic[allele][c] * dic[allele][d])
                p2 = 2 * dic[allele]["pmin"] * (2 - dic[allele]["pmin"]) * dic[allele][b] * (2 - dic[allele][b]) - \
                     2 * dic[allele]["pmin"] * dic[allele][b] * 2 * dic[allele]["pmin"] * dic[allele][b]
            elif a in dic[allele] and b not in dic[allele] and c in dic[allele] and d in dic[allele]:
                p1 = (2 * 2 * dic[allele][a] * dic[allele][c] * 2 * dic[allele]["pmin"] * dic[allele][d] + 2 *
                      2 * dic[allele][a] * dic[allele][d] * 2 * dic[allele]["pmin"] * dic[allele][c]) / \
                     (2 * dic[allele][c] * (2 - dic[allele][c]) * dic[allele][d] * (2 - dic[allele][d]) -
                      2 * dic[allele][c] * dic[allele][d] * 2 * dic[allele][c] * dic[allele][d])
                p2 = 2 * dic[allele][a] * (2 - dic[allele][a]) * dic[allele]["pmin"] * (2 - dic[allele]["pmin"]) - \
                     2 * dic[allele][a] * dic[allele]["pmin"] * 2 * dic[allele][a] * dic[allele]["pmin"]
            elif a in dic[allele] and b in dic[allele] and c not in dic[allele] and d in dic[allele]:
                p1 = (2 * 2 * dic[allele][a] * dic[allele]["pmin"] * 2 * dic[allele][b] * dic[allele][d] + 2 *
                      2 * dic[allele][a] * dic[allele][d] * 2 * dic[allele][b] * dic[allele]["pmin"]) / \
                     (2 * dic[allele]["pmin"] * (2 - dic[allele]["pmin"]) * dic[allele][d] * (2 - dic[allele][d]) -
                      2 * dic[allele]["pmin"] * dic[allele][d] * 2 * dic[allele]["pmin"] * dic[allele][d])
                p2 = 2 * dic[allele][a] * (2 - dic[allele][a]) * dic[allele][b] * (2 - dic[allele][b]) - \
                     2 * dic[allele][a] * dic[allele][b] * 2 * dic[allele][a] * dic[allele][b]
            elif a in dic[allele] and b in dic[allele] and c in dic[allele] and d not in dic[allele]:
                p1 = (2 * 2 * dic[allele][a] * dic[allele][c] * 2 * dic[allele][b] * dic[allele]["pmin"] + 2 *
                      2 * dic[allele][a] * dic[allele]["pmin"] * 2 * dic[allele][b] * dic[allele][c]) / \
                     (2 * dic[allele][c] * (2 - dic[allele][c]) * dic[allele]["pmin"] * (2 - dic[allele]["pmin"]) -
                      2 * dic[allele][c] * dic[allele]["pmin"] * 2 * dic[allele][c] * dic[allele]["pmin"])
                p2 = 2 * dic[allele][a] * (2 - dic[allele][a]) * dic[allele][b] * (2 - dic[allele][b]) - \
                     2 * dic[allele][a] * dic[allele][b] * 2 * dic[allele][a] * dic[allele][b]
            elif a not in dic[allele] and b not in dic[allele] and c in dic[allele] and d in dic[allele]:
                p1 = (2 * 2 * dic[allele]["pmin"] * dic[allele][c] * 2 * dic[allele]["pmin"] * dic[allele][d] + 2 *
                      2 * dic[allele]["pmin"] * dic[allele][d] * 2 * dic[allele]["pmin"] * dic[allele][c]) / \
                     (2 * dic[allele][c] * (2 - dic[allele][c]) * dic[allele][d] * (2 - dic[allele][d]) -
                      2 * dic[allele][c] * dic[allele][d] * 2 * dic[allele][c] * dic[allele][d])
                p2 = 2 * dic[allele]["pmin"] * (2 - dic[allele]["pmin"]) * dic[allele]["pmin"] * (2 - dic[allele]["pmin"]) - \
                     2 * dic[allele]["pmin"] * dic[allele]["pmin"] * 2 * dic[allele]["pmin"] * dic[allele]["pmin"]
            elif a not in dic[allele] and b in dic[allele] and c not in dic[allele] and d in dic[allele]:
                p1 = (2 * 2 * dic[allele]["pmin"] * dic[allele]["pmin"] * 2 * dic[allele][b] * dic[allele][d] + 2 *
                      2 * dic[allele]["pmin"] * dic[allele][d] * 2 * dic[allele][b] * dic[allele]["pmin"]) / \
                     (2 * dic[allele]["pmin"] * (2 - dic[allele]["pmin"]) * dic[allele][d] * (2 - dic[allele][d]) -
                      2 * dic[allele]["pmin"] * dic[allele][d] * 2 * dic[allele]["pmin"] * dic[allele][d])
                p2 = 2 * dic[allele]["pmin"] * (2 - dic[allele]["pmin"]) * dic[allele][b] * (2 - dic[allele][b]) - \
                     2 * dic[allele]["pmin"] * dic[allele][b] * 2 * dic[allele]["pmin"] * dic[allele][b]
            elif a not in dic[allele] and b in dic[allele] and c in dic[allele] and d not in dic[allele]:
                p1 = (2 * 2 * dic[allele]["pmin"] * dic[allele][c] * 2 * dic[allele][b] * dic[allele]["pmin"] + 2 *
                      2 * dic[allele]["pmin"] * dic[allele]["pmin"] * 2 * dic[allele][b] * dic[allele][c]) / \
                     (2 * dic[allele][c] * (2 - dic[allele][c]) * dic[allele]["pmin"] * (2 - dic[allele]["pmin"]) -
                      2 * dic[allele][c] * dic[allele]["pmin"] * 2 * dic[allele][c] * dic[allele]["pmin"])
                p2 = 2 * dic[allele]["pmin"] * (2 - dic[allele]["pmin"]) * dic[allele][b] * (2 - dic[allele][b]) - \
                     2 * dic[allele]["pmin"] * dic[allele][b] * 2 * dic[allele]["pmin"] * dic[allele][b]
            elif a in dic[allele] and b not in dic[allele] and c not in dic[allele] and d in dic[allele]:
                p1 = (2 * 2 * dic[allele][a] * dic[allele]["pmin"] * 2 * dic[allele]["pmin"] * dic[allele][d] + 2 *
                      2 * dic[allele][a] * dic[allele][d] * 2 * dic[allele]["pmin"] * dic[allele]["pmin"]) / \
                     (2 * dic[allele]["pmin"] * (2 - dic[allele]["pmin"]) * dic[allele][d] * (2 - dic[allele][d]) -
                      2 * dic[allele]["pmin"] * dic[allele][d] * 2 * dic[allele]["pmin"] * dic[allele][d])
                p2 = 2 * dic[allele][a] * (2 - dic[allele][a]) * dic[allele]["pmin"] * (2 - dic[allele]["pmin"]) - \
                     2 * dic[allele][a] * dic[allele]["pmin"] * 2 * dic[allele][a] * dic[allele]["pmin"]
            elif a in dic[allele] and b not in dic[allele] and c in dic[allele] and d not in dic[allele]:
                p1 = (2 * 2 * dic[allele][a] * dic[allele][c] * 2 * dic[allele]["pmin"] * dic[allele]["pmin"] + 2 *
                      2 * dic[allele][a] * dic[allele]["pmin"] * 2 * dic[allele]["pmin"] * dic[allele][c]) / \
                     (2 * dic[allele][c] * (2 - dic[allele][c]) * dic[allele]["pmin"] * (2 - dic[allele]["pmin"]) -
                      2 * dic[allele][c] * dic[allele]["pmin"] * 2 * dic[allele][c] * dic[allele]["pmin"])
                p2 = 2 * dic[allele][a] * (2 - dic[allele][a]) * dic[allele]["pmin"] * (2 - dic[allele]["pmin"]) - \
                     2 * dic[allele][a] * dic[allele]["pmin"] * 2 * dic[allele][a] * dic[allele]["pmin"]
            elif a in dic[allele] and b in dic[allele] and c not in dic[allele] and d not in dic[allele]:
                p1 = (2 * 2 * dic[allele][a] * dic[allele]["pmin"] * 2 * dic[allele][b] * dic[allele]["pmin"] + 2 *
                      2 * dic[allele][a] * dic[allele]["pmin"] * 2 * dic[allele][b] * dic[allele]["pmin"]) / \
                     (2 * dic[allele]["pmin"] * (2 - dic[allele]["pmin"]) * dic[allele]["pmin"] * (2 - dic[allele]["pmin"]) -
                      2 * dic[allele]["pmin"] * dic[allele]["pmin"] * 2 * dic[allele]["pmin"] * dic[allele]["pmin"])
                p2 = 2 * dic[allele][a] * (2 - dic[allele][a]) * dic[allele][b] * (2 - dic[allele][b]) - \
                     2 * dic[allele][a] * dic[allele][b] * 2 * dic[allele][a] * dic[allele][b]
            elif a not in dic[allele] and b not in dic[allele] and c not in dic[allele] and d in dic[allele]:
                p1 = (2 * 2 * dic[allele]["pmin"] * dic[allele]["pmin"] * 2 * dic[allele]["pmin"] * dic[allele][d] + 2 *
                      2 * dic[allele]["pmin"] * dic[allele][d] * 2 * dic[allele]["pmin"] * dic[allele]["pmin"]) / \
                     (2 * dic[allele]["pmin"] * (2 - dic[allele]["pmin"]) * dic[allele][d] * (2 - dic[allele][d]) -
                      2 * dic[allele]["pmin"] * dic[allele][d] * 2 * dic[allele]["pmin"] * dic[allele][d])
                p2 = 2 * dic[allele]["pmin"] * (2 - dic[allele]["pmin"]) * dic[allele]["pmin"] * (2 - dic[allele]["pmin"]) - \
                     2 * dic[allele]["pmin"] * dic[allele]["pmin"] * 2 * dic[allele]["pmin"] * dic[allele]["pmin"]
            elif a not in dic[allele] and b not in dic[allele] and c in dic[allele] and d not in dic[allele]:
                p1 = (2 * 2 * dic[allele]["pmin"] * dic[allele][c] * 2 * dic[allele]["pmin"] * dic[allele]["pmin"] + 2 *
                      2 * dic[allele]["pmin"] * dic[allele]["pmin"] * 2 * dic[allele]["pmin"] * dic[allele][c]) / \
                     (2 * dic[allele][c] * (2 - dic[allele][c]) * dic[allele]["pmin"] * (2 - dic[allele]["pmin"]) -
                      2 * dic[allele][c] * dic[allele]["pmin"] * 2 * dic[allele][c] * dic[allele]["pmin"])
                p2 = 2 * dic[allele]["pmin"] * (2 - dic[allele]["pmin"]) * dic[allele]["pmin"] * (2 - dic[allele]["pmin"]) - \
                     2 * dic[allele]["pmin"] * dic[allele]["pmin"] * 2 * dic[allele]["pmin"] * dic[allele]["pmin"]
            elif a not in dic[allele] and b in dic[allele] and c not in dic[allele] and d not in dic[allele]:
                p1 = (2 * 2 * dic[allele]["pmin"] * dic[allele]["pmin"] * 2 * dic[allele][b] * dic[allele]["pmin"] + 2 *
                      2 * dic[allele]["pmin"] * dic[allele]["pmin"] * 2 * dic[allele][b] * dic[allele]["pmin"]) / \
                     (2 * dic[allele]["pmin"] * (2 - dic[allele]["pmin"]) * dic[allele]["pmin"] * (2 - dic[allele]["pmin"]) -
                      2 * dic[allele]["pmin"] * dic[allele]["pmin"] * 2 * dic[allele]["pmin"] * dic[allele]["pmin"])
                p2 = 2 * dic[allele]["pmin"] * (2 - dic[allele]["pmin"]) * dic[allele][b] * (2 - dic[allele][b]) - \
                     2 * dic[allele]["pmin"] * dic[allele][b] * 2 * dic[allele]["pmin"] * dic[allele][b]
            elif a in dic[allele] and b not in dic[allele] and c not in dic[allele] and d not in dic[allele]:
                p1 = (2 * 2 * dic[allele][a] * dic[allele]["pmin"] * 2 * dic[allele]["pmin"] * dic[allele]["pmin"] + 2 *
                      2 * dic[allele][a] * dic[allele]["pmin"] * 2 * dic[allele]["pmin"] * dic[allele]["pmin"]) / \
                     (2 * dic[allele]["pmin"] * (2 - dic[allele]["pmin"]) * dic[allele]["pmin"] * (2 - dic[allele]["pmin"]) -
                      2 * dic[allele]["pmin"] * dic[allele]["pmin"] * 2 * dic[allele]["pmin"] * dic[allele]["pmin"])
                p2 = 2 * dic[allele][a] * (2 - dic[allele][a]) * dic[allele]["pmin"] * (2 - dic[allele]["pmin"]) - \
                     2 * dic[allele][a] * dic[allele]["pmin"] * 2 * dic[allele][a] * dic[allele]["pmin"]
            elif a not in dic[allele] and b not in dic[allele] and c not in dic[allele] and d not in dic[allele]:
                p1 = (2 * 2 * dic[allele]["pmin"] * dic[allele]["pmin"] * 2 * dic[allele]["pmin"] * dic[allele]["pmin"] + 2 *
                      2 * dic[allele]["pmin"] * dic[allele]["pmin"] * 2 * dic[allele]["pmin"] * dic[allele]["pmin"]) / \
                     (2 * dic[allele]["pmin"] * (2 - dic[allele]["pmin"]) * dic[allele]["pmin"] * (2 - dic[allele]["pmin"]) -
                      2 * dic[allele]["pmin"] * dic[allele]["pmin"] * 2 * dic[allele]["pmin"] * dic[allele]["pmin"])
                p2 = 2 * dic[allele]["pmin"] * (2 - dic[allele]["pmin"]) * dic[allele]["pmin"] * (2 - dic[allele]["pmin"]) - \
                     2 * dic[allele]["pmin"] * dic[allele]["pmin"] * 2 * dic[allele]["pmin"] * dic[allele]["pmin"]
            print("s1 != s2 and k1 != s1 and k1 != s2 and k2 != s1 and k2 != s2", allele, p1, p2)
    return p1, p2


def main():
    start = time.time()
    with open("./config_file_siblings.yaml", "r") as yaml_file:
        config = yaml.load(yaml_file, Loader=yaml.FullLoader)
    # Create dataframe with offsprings
    parents_df = pd.read_excel(config["main_table_path"])
    columns_list = parents_df.columns.values.tolist()
    columns_list = columns_list[:-1]
    columns_list.append("Family_ID")
    columns_list.append("Status")
    # Add columns for LR calculations (for potential brothers and sisters)
    columns_list.append("LR_ref_old_codis")
    columns_list.append("LR_rus_old_codis")
    columns_list.append("LR_ref_new_codis")
    columns_list.append("LR_rus_new_codis")

    offsprings_df = pd.DataFrame(columns=columns_list)
    for pop in config["populations"]:
        ref_dict = empty_freq_table(config["main_table_path"])
        ref_dict = calculate_frequencies(pop, ref_dict, config["main_table_path"])
        new_pair_df = pd.DataFrame(columns=columns_list, index=[0, 1])
        temp_parents_df = parents_df.loc[(parents_df['population'] == pop)]
        index_list = temp_parents_df.index.values.tolist()
        index_list_f = index_list[:1000]
        index_list_m = index_list[1000:]
        list_for_pairs = []
        n = copy.deepcopy(config["number_of_pairs"])
        while n != 0:
            father = random.choice(index_list_f)
            mother = random.choice(index_list_m)
            new_pair = [father, mother]
            pair = set(new_pair)
            if pair not in list_for_pairs and father != mother:
                list_for_pairs.append(pair)
                pair_df = temp_parents_df.loc[[father, mother]]    # create df for one random pair
                kids_df = pd.DataFrame(columns=columns_list, index=[i for i in range(0, config["kids"] + 1)])
                for k in range(config["kids"] + 1):
                    if k != 0:
                        kids_df.iloc[k]["Status"] = "real_sibling"
                    else:
                        kids_df.iloc[k]["Status"] = "kid"
                    kids_df.iloc[k]["№"] = (str(pair_df.iloc[0]['№']) + "-" + str(pair_df.iloc[1]['№']) + "-" +
                                            str(k + 1))
                    for elem in range(1, len(columns_list) - 7, 2):
                        f_alleles = [pair_df.iloc[0][columns_list[elem]], pair_df.iloc[0][columns_list[elem + 1]]]
                        f_allele_random = random.choice(f_alleles)
                        m_alleles = [pair_df.iloc[1][columns_list[elem]], pair_df.iloc[1][columns_list[elem + 1]]]
                        m_allele_random = random.choice(m_alleles)
                        possible_alleles = [f_allele_random, m_allele_random]
                        kid_alleles = random.sample(possible_alleles, 2)
                        kid_alleles.sort()
                        kids_df.iloc[k][columns_list[elem]] = kid_alleles[0]
                        kids_df.iloc[k][columns_list[elem + 1]] = kid_alleles[1]
                        f_alleles.sort()
                        m_alleles.sort()
                    if k == 0:
                        for elem in range(len(columns_list) - 6):
                            new_pair_df.iloc[0][columns_list[elem]] = pair_df.iloc[0][columns_list[elem]]
                            new_pair_df.iloc[1][columns_list[elem]] = pair_df.iloc[1][columns_list[elem]]
                        new_pair_df.iloc[0]["Status"] = "real_father"
                        new_pair_df.iloc[1]["Status"] = "real_mother"
                        new_pair_df.iloc[0]["Family_ID"] = (str(pair_df.iloc[0]['№']) + "-" + str(pair_df.iloc[1]['№']))
                        new_pair_df.iloc[1]["Family_ID"] = (str(pair_df.iloc[0]['№']) + "-" + str(pair_df.iloc[1]['№']))
                        offsprings_df = pd.concat([offsprings_df, new_pair_df], ignore_index=True)
                    kids_df.iloc[k]["Family_ID"] = (str(pair_df.iloc[0]['№']) + "-" + str(pair_df.iloc[1]['№']))
                    kids_df.iloc[k]["population"] = pop
                    offsprings_df = pd.concat([offsprings_df, kids_df.loc[[k]]], ignore_index=True)

                # Add data for potential siblings with mutation:
                siblings_mut_df = pd.DataFrame(columns=columns_list, index=[i for i in range(0, config["sibl_mut"])])
                for k in range(config["sibl_mut"]):
                    siblings_mut_df.iloc[k]["Status"] = "sibling_with_mutation"
                    siblings_mut_df.iloc[k]["№"] = (str(pair_df.iloc[0]['№']) + "-" + str(pair_df.iloc[1]['№']) + "-"
                                                    + str(config["kids"] + k + 2))
                    rand_n = 0
                    randstr = []
                    while rand_n != 7:
                        rand = random.randint(1, len(columns_list) - 7)
                        if rand % 2 == 1 and rand not in randstr:
                            randstr.append(rand)
                            rand_n += 1
                    for elem in range(1, len(columns_list) - 7, 2):
                        f_alleles = [pair_df.iloc[0][columns_list[elem]], pair_df.iloc[0][columns_list[elem + 1]]]
                        f_allele_random = random.choice(f_alleles)
                        m_alleles = [pair_df.iloc[1][columns_list[elem]], pair_df.iloc[1][columns_list[elem + 1]]]
                        m_allele_random = random.choice(m_alleles)
                        possible_alleles = [f_allele_random, m_allele_random]
                        kid_alleles = random.sample(possible_alleles, 2)
                        kid_alleles.sort()
                        if elem in randstr:
                            a = random.choice([0, 1])
                            if a == 0:
                                siblings_mut_df.iloc[k][columns_list[elem]] = kid_alleles[0] - 1
                                siblings_mut_df.iloc[k][columns_list[elem + 1]] = kid_alleles[1] + 1
                            elif a == 1:
                                siblings_mut_df.iloc[k][columns_list[elem]] = kid_alleles[0]
                                siblings_mut_df.iloc[k][columns_list[elem + 1]] = kid_alleles[1] + 1
                        else:
                            siblings_mut_df.iloc[k][columns_list[elem]] = kid_alleles[0]
                            siblings_mut_df.iloc[k][columns_list[elem + 1]] = kid_alleles[1]
                        f_alleles.sort()
                        m_alleles.sort()
                    siblings_mut_df.iloc[k]["Family_ID"] = (str(pair_df.iloc[0]['№']) + "-" + str(pair_df.iloc[1]['№']))
                    siblings_mut_df.iloc[k]["population"] = pop
                    offsprings_df = pd.concat([offsprings_df, siblings_mut_df.loc[[k]]], ignore_index=True)

                # Add data for nonsiblings:
                nonsiblings_df = pd.DataFrame(columns=columns_list, index=[i for i in range(0, config["nonsibling"])])
                for k in range(config["nonsibling"]):
                    nonsiblings_df.iloc[k]["Status"] = "nonsibling"
                    nonsiblings_df.iloc[k]["№"] = (str(pair_df.iloc[0]['№']) + "-" + str(pair_df.iloc[1]['№']) + "-"
                                                   + str(config["kids"] + config["sibl_mut"] + k + 2))
                    for elem in range(1, len(columns_list) - 7, 2):
                        f_alleles = [pair_df.iloc[0][columns_list[elem]] + 1, pair_df.iloc[0][columns_list[elem + 1]] + 1]
                        f_allele_random = random.choice(f_alleles)
                        m_alleles = [pair_df.iloc[1][columns_list[elem]] + 1, pair_df.iloc[1][columns_list[elem + 1]] + 1]
                        m_allele_random = random.choice(m_alleles)
                        possible_alleles = [f_allele_random, m_allele_random]
                        kid_alleles = random.sample(possible_alleles, 2)
                        kid_alleles.sort()
                        nonsiblings_df.iloc[k][columns_list[elem]] = kid_alleles[0]
                        nonsiblings_df.iloc[k][columns_list[elem + 1]] = kid_alleles[1]
                        f_alleles.sort()
                        m_alleles.sort()
                    nonsiblings_df.iloc[k]["Family_ID"] = (str(pair_df.iloc[0]['№']) + "-" + str(pair_df.iloc[1]['№']))
                    nonsiblings_df.iloc[k]["population"] = pop
                    offsprings_df = pd.concat([offsprings_df, nonsiblings_df.loc[[k]]], ignore_index=True)
            n -= 1

            #######################
            # Count LR
            df_for_lr = offsprings_df.loc[(offsprings_df['Family_ID'] == new_pair_df.iloc[0]["Family_ID"])]
            proband_df = df_for_lr.iloc[[2]]  # df with 1st, 2nd or 3rd child's data
            for ii in range(3, len(df_for_lr.index)):
                hyp1_ref_old = []
                hyp2_ref_old = []
                hyp1_rus_old = []
                hyp2_rus_old = []
                hyp1_ref_new = []
                hyp2_ref_new = []
                hyp1_rus_new = []
                hyp2_rus_new = []
                sibl_df = df_for_lr.iloc[[ii]]
                for loc in range(1, len(columns_list) - 7, 2):
                    el = columns_list[loc].split('_')[0]
                    kid = [proband_df.iloc[0][columns_list[loc]], proband_df.iloc[0][columns_list[loc + 1]]]
                    sibl = [sibl_df.iloc[0][columns_list[loc]], sibl_df.iloc[0][columns_list[loc + 1]]]
                    kid.sort()
                    sibl.sort()
                    p1, p2, p3, p4 = hyp_sibl(sibl[0], sibl[1], kid[0], kid[1], ref_dict, config["rus_dict"], el)
                    if el in config["codis_old"]:
                        hyp1_ref_old.append(p1)
                        hyp2_ref_old.append(p2)
                        hyp1_rus_old.append(p3)
                        hyp2_rus_old.append(p4)
                    if el in config["codis_new"]:
                        hyp1_ref_new.append(p1)
                        hyp2_ref_new.append(p2)
                        hyp1_rus_new.append(p3)
                        hyp2_rus_new.append(p4)
                print(hyp1_ref_old)
                print(hyp2_ref_old)
                print(hyp1_rus_old)
                print(hyp1_rus_old)
                print(hyp1_ref_new)
                print(hyp2_ref_new)
                print(hyp1_rus_new)
                print(hyp1_rus_new)
                multiplication_1 = 1
                multiplication_2 = 1
                multiplication_3 = 1
                multiplication_4 = 1
                multiplication_5 = 1
                multiplication_6 = 1
                multiplication_7 = 1
                multiplication_8 = 1
                for m in hyp1_ref_old:
                    multiplication_1 *= m
                for m in hyp2_ref_old:
                    multiplication_2 *= m
                for m in hyp1_rus_old:
                    multiplication_3 *= m
                for m in hyp2_rus_old:
                    multiplication_4 *= m
                for m in hyp1_ref_new:
                    multiplication_5 *= m
                for m in hyp2_ref_new:
                    multiplication_6 *= m
                for m in hyp1_rus_new:
                    multiplication_7 *= m
                for m in hyp2_rus_new:
                    multiplication_8 *= m
                offsprings_df.iloc[df_for_lr.index[ii]]["LR_ref_old_codis"] = multiplication_1 / multiplication_2
                offsprings_df.iloc[df_for_lr.index[ii]]["LR_rus_old_codis"] = multiplication_3 / multiplication_4
                offsprings_df.iloc[df_for_lr.index[ii]]["LR_ref_new_codis"] = multiplication_5 / multiplication_6
                offsprings_df.iloc[df_for_lr.index[ii]]["LR_rus_new_codis"] = multiplication_7 / multiplication_8
    offsprings_df.to_excel("SIBLINGS_OUTPUT_TEST.xlsx", index=True)
    print(round(time.time() - start, 2), 's')


if __name__ == "__main__":
    main()
