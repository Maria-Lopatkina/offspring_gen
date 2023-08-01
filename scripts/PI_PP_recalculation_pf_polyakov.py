import pandas as pd
import time


def empty_freq_table():  # Create empty dictionaries for allele frequencies
    freq_df = pd.read_excel("ASTR22_main.xlsx")
    f_columns = freq_df.columns.values.tolist()
    key_dict = []
    for el in range(1, len(f_columns) - 2, 2):
        key_dict.append(list(f_columns[el].split("_"))[0])
    freq_dict = dict.fromkeys(key_dict)
    return freq_dict


def calculate_frequencies(p, d):  # Calculate alleles frequencies for three populations
    freq_df = pd.read_excel("ASTR22_main.xlsx")
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
    return d


def father_trio(f1, f2, m1, m2, k1, k2):
    if (f1 == k1 and m1 == k2) or (f2 == k1 and m2 == k2) or (f1 == k1 and m2 == k2) or (f2 == k1 and m1 == k2) or \
            (f1 == k2 and m1 == k1) or (f2 == k2 and m2 == k1) or (f1 == k2 and m2 == k1) or (f2 == k2 and m1 == k1):
        return True
    return False


def father_duo(f1, f2, k1, k2):
    for i in f1, f2:
        for j in k1, k2:
            if i == j:
                return True
    return False


def equations(a1, a2, n, dic, allele):
    p = None
    if n == 1:
        if a1 in dic[allele]:
            p = 1 / (dic[allele][a1])
        else:
            p = 1 / (dic[allele]["pmin"])
    if n == 2:
        if a1 in dic[allele]:
            p = 1 / (2 * dic[allele][a1])
        else:
            p = 1 / (2 * dic[allele]["pmin"])
    if n == 3:
        if a1 in dic[allele] and a2 in dic[allele]:
            p = 1 / (dic[allele][a1] + dic[allele][a2])
        elif a1 not in dic[allele] and a2 in dic[allele]:
            p = 1 / (dic[allele]["pmin"] + dic[allele][a2])
        elif a1 in dic[allele] and a2 not in dic[allele]:
            p = 1 / (dic[allele][a1] + dic[allele]["pmin"])
        elif a1 in dic[allele] and a2 in dic[allele]:
            p = 1 / (dic[allele]["pmin"] + dic[allele]["pmin"])
    if n == 4:
        if a1 in dic[allele] and a2 in dic[allele]:
            p = 1 / (2 * dic[allele][a1] + 2 * dic[allele][a2])
        elif a1 not in dic[allele] and a2 in dic[allele]:
            p = 1 / (2 * dic[allele]["pmin"] + 2 * dic[allele][a2])
        elif a1 in dic[allele] and a2 not in dic[allele]:
            p = 1 / (2 * dic[allele][a1] + 2 * dic[allele]["pmin"])
        elif a1 in dic[allele] and a2 in dic[allele]:
            p = 1 / (2 * dic[allele]["pmin"] + 2 * dic[allele]["pmin"])
    if n == 5:
        if a1 in dic[allele] and a2 in dic[allele]:
            p = (dic[allele][a1] + dic[allele][a2]) / (4 * dic[allele][a1] * dic[allele][a2])
        elif a1 not in dic[allele] and a2 in dic[allele]:
            p = (dic[allele]["pmin"] + dic[allele][a2]) / (4 * dic[allele]["pmin"] * dic[allele][a2])
        elif a1 in dic[allele] and a2 not in dic[allele]:
            p = (dic[allele][a1] + dic[allele]["pmin"]) / (4 * dic[allele][a1] * dic[allele]["pmin"])
        elif a1 in dic[allele] and a2 in dic[allele]:
            p = (dic[allele]["pmin"] + dic[allele]["pmin"]) / (4 * dic[allele]["pmin"] * dic[allele]["pmin"])
    if n == 6:
        if a1 in dic[allele]:
            p = 1 / (4 * dic[allele][a1])
        else:
            p = 1 / (4 * dic[allele]["pmin"])
    return p


def p_calc_trio_wo_mut(f1, f2, m1, m2, k1, k2, dic, allele):
    if k1 == k2:
        if f1 == f2:
            return equations(k1, None, 1, dic, allele)
        else:
            return equations(k1, None, 2, dic, allele)
    else:
        if k1 == m1 and k2 == m2:
            for i in f1, f2:
                if i not in [k1, k2]:
                    return equations(k1, k2, 4, dic, allele)
            return equations(k1, k2, 3, dic, allele)
        else:
            if f1 == f2:
                return equations(f1, None, 1, dic, allele)
            else:
                for i in k1, k2:
                    if i not in [m1, m2] and i in [f1, f2]:
                        return equations(i, None, 2, dic, allele)


def p_calc_duo_wo_mut(f1, f2, k1, k2, dic, allele):
    if k1 == k2:
        if f1 == f2:
            return equations(k1, None, 1, dic, allele)
        else:
            return equations(k1, None, 2, dic, allele)
    else:
        if f1 == f2:
            return equations(f1, None, 2, dic, allele)
        else:
            if k1 == f1 and k2 == f2:
                return equations(k1, k2, 5, dic, allele)
            else:
                return equations(k1, None, 6, dic, allele)


def p_calc_trio_mut(f1, f2, m1, m2, k1, k2, dic, allele, mut):
    if k1 == k2:
        if f1 == f2:
            return equations(k1, None, 1, dic, allele)
        else:
            return equations(k1, None, 2, dic, allele)
    else:
        if k1 == m1 and k2 == m2:
            for i in f1, f2:
                if i not in [k1, k2]:
                    return equations(k1, k2, 4, dic, allele)
            return equations(k1, k2, 3, dic, allele)
        else:
            if f1 == f2:
                return equations(f1, None, 1, dic, allele)
            else:
                for i in k1, k2:
                    if i not in [m1, m2] and i in [f1, f2]:
                        return equations(i, None, 2, dic, allele)


def p_calc_duo_mut(f1, f2, k1, k2, dic, allele, mut):
    if k1 == k2:
        if f1 == f2:
            return equations(k1, None, 1, dic, allele)
        else:
            return equations(k1, None, 2, dic, allele)
    else:
        if f1 == f2:
            return equations(f1, None, 2, dic, allele)
        else:
            if k1 == f1 and k2 == f2:
                return equations(k1, k2, 5, dic, allele)
            else:
                return equations(k1, None, 6, dic, allele)




def main():
    start = time.time()
    mut_freq = float(input("Enter the frequency of STR-locus mutation: "))
    population = input("Enter the population name: ")
    # Create dictionary with alleles frequencies
    ref_dict = empty_freq_table()
    ref_dict = calculate_frequencies(population, ref_dict)
    rus_dict = {"CSF1PO": {7: 0.0004, 8: 0.0012, 9: 0.0438, 10: 0.2824, 11: 0.2882, 12: 0.3062, 13: 0.0617, 14: 0.0127,
                           15: 0.0032, 16: 0.0003, "pmin": 0.0003},
                "D3S1358": {10: 0.0003, 11: 0.0008, 12: 0.0005, 13: 0.0014, 14: 0.1073, 15: 0.2721, 16: 0.268,
                            17: 0.2148, 18: 0.1244, 19: 0.0096, 20: 0.0009, 21: 0.0003, "pmin": 0.0003},
                "D7S820": {6: 0.0003, 7: 0.012, 7.3: 0.0003, 8: 0.1701, 8.1: 0.0003, 9: 0.1261, 10: 0.2763,
                           10.1: 0.0003, 10.3: 0.0003, 11: 0.2223, 11.1: 0.0033, 11.2: 0.0003, 11.3: 0.0003,
                           12: 0.1549, 12.1: 0.0004, 13: 0.0296, 13.1: 0.0003, 14: 0.0037, 15: 0.0003, "pmin": 0.0003},
                "D8S1179": {6: 0.0003, 8: 0.0082, 9: 0.0083, 10: 0.0702, 11: 0.0725, 11.3: 0.0003, 12: 0.1596,
                            12.3: 0.0003, 13: 0.3284, 14: 0.2182, 15: 0.1019, 16: 0.0277, 17: 0.0046, 18: 0.0003,
                            "pmin": 0.0003},
                "D16S539": {7: 0.0003, 8: 0.0118, 9: 0.119, 10: 0.0709, 11: 0.2732, 12: 0.3105, 13: 0.1823, 14: 0.0295,
                            15: 0.0024, 16: 0.0003, "pmin": 0.0003},
                "D18S51": {9: 0.0003, 10: 0.0081, 11: 0.0148, 12: 0.0906, 12.3: 0.0003, 13: 0.1142, 14: 0.1647,
                           14.2: 0.0003, 15: 0.1643, 16: 0.1576, 16.2: 0.0003, 17: 0.1251, 18: 0.0747, 19: 0.0422,
                           20: 0.0207, 21: 0.0131, 21.1: 0.0003, 22: 0.0058, 23: 0.0021, 24: 0.0012, 25: 0.0003,
                           30: 0.0003, "pmin": 0.0003},
                "D21S11": {25: 0.0007, 25.2: 0.0003, 26: 0.0026, 27: 0.0212, 27.2: 0.0003, 28: 0.1557, 28.2: 0.0008,
                           28.3: 0.0003, 29: 0.2048, 29.1: 0.0003, 29.2: 0.0014, 29.3: 0.0003, 30: 0.2362, 30.2: 0.056,
                           31: 0.0698, 31.2: 0.0899, 32: 0.0103, 32.1: 0.0003, 32.2: 0.1043, 33: 0.0012, 33.1: 0.0003,
                           33.2: 0.0392, 34: 0.0003, 34.2: 0.0045, 35.2: 0.0003, "pmin": 0.0003},
                "FGA": {16: 0.0006, 17: 0.0007, 18: 0.0152, 18.2: 0.0003, 19: 0.0726, 19.2: 0.0003, 20: 0.121,
                        20.2: 0.0015, 21: 0.1675, 21.2: 0.0024, 21.3: 0.0003, 22: 0.1966, 22.1: 0.0003, 22.2: 0.0107,
                        22.3: 0.0003, 23: 0.1421, 23.2: 0.0066, 24: 0.1424, 24.1: 0.0003, 24.2: 0.0017, 25: 0.0823,
                        25.1: 0.0003, 25.2: 0.0004, 26: 0.0276, 27: 0.0061, 28: 0.0009, 29: 0.0003, 33: 0.0003,
                        "pmin": 0.0003},
                "TH01": {5: 0.0005, 6: 0.2295, 7: 0.1626, 7.3: 0.0003, 8: 0.1068, 8.3: 0.0003, 9: 0.2065, 9.3: 0.2845,
                         10: 0.0086, 10.3: 0.0006, 11: 0.0003, "pmin": 0.0003},
                "TPOX": {5: 0.0003, 6: 0.0003, 7: 0.0013, 8: 0.5551, 9: 0.088, 10: 0.0553, 11: 0.2666, 12: 0.0323,
                         13: 0.0011, "pmin": 0.0003},
                "vWA": {12: 0.0003, 13: 0.0022, 14: 0.0956, 15: 0.1058, 16: 0.1882, 17: 0.2905, 18: 0.2211, 19: 0.082,
                        20: 0.0131, 21: 0.0013, "pmin": 0.0003},
                "D1S1656": {8: 0.0004, 9: 0.0003, 10: 0.0016, 11: 0.1033, 12: 0.1192, 13: 0.0698, 13.2: 0.0003,
                            14: 0.0763, 14.2: 0.0003, 14.3: 0.0011, 15: 0.1403, 15.2: 0.0003, 15.3: 0.0362, 16: 0.1274,
                            16.1: 0.0003, 16.3: 0.0418, 17: 0.0604, 17.1: 0.0003, 17.3: 0.1399, 18: 0.0072,
                            18.3: 0.0594, 19: 0.0003, 19.1: 0.0003, 19.3: 0.0134, 20.3: 0.0013, "pmin": 0.0003},
                "D2S441": {8: 0.0008, 9: 0.0015, 9.1: 0.0006, 10: 0.2266, 10.1: 0.0003, 10.3: 0.0003, 10.4: 0.0003,
                           11: 0.3561, 11.3: 0.0552, 12: 0.0423, 12.3: 0.0012, 13: 0.0227, 14: 0.256, 15: 0.0329,
                           16: 0.0039, "pmin": 0.0003},
                "D2S1338": {12: 0.0003, 13: 0.0003, 14: 0.0003, 15: 0.0013, 16: 0.0318, 17: 0.1815, 18: 0.0897,
                            19: 0.1369, 20: 0.1393, 21: 0.034, 22: 0.0277, 23: 0.1003, 24: 0.1119, 25: 0.12, 26: 0.0213,
                            27: 0.0027, 28: 0.0009, 29: 0.0003, "pmin": 0.0003},
                "D10S1248": {8: 0.0003, 9: 0.0003, 10: 0.0003, 11: 0.003, 12: 0.0228, 13: 0.2347, 14: 0.3307, 15: 0.231,
                             16: 0.1394, 17: 0.0354, 18: 0.0022, 19: 0.0003, "pmin": 0.0003},
                "D12S391": {14: 0.0003, 15: 0.03, 16: 0.0214, 17: 0.1003, 17.3: 0.0117, 18: 0.1996, 18.3: 0.0172,
                            19: 0.1368, 19.1: 0.0003, 19.2: 0.0003, 19.3: 0.0076, 20: 0.1266, 20.3: 0.0009, 21: 0.1094,
                            22: 0.1086, 22.3: 0.0003, 23: 0.0813, 24: 0.0333, 25: 0.0115, 26: 0.0028, 27: 0.0007,
                            "pmin": 0.0003},
                "D19S433": {6.2: 0.0003, 9: 0.0003, 10: 0.0004, 10.2: 0.0003, 11: 0.0049, 11.2: 0.0003, 12: 0.0847,
                            12.2: 0.0023, 13: 0.2307, 13.2: 0.0223, 14: 0.3331, 14.2: 0.0295, 15: 0.1562, 15.2: 0.0577,
                            16: 0.0398, 16.2: 0.0232, 17: 0.0032, 17.2: 0.0067, 18: 0.0003, 18.2: 0.0039, 19.2: 0.0004,
                            "pmin": 0.0003},
                "D22S1045": {8: 0.0003, 9: 0.0003, 10: 0.0009, 11: 0.1724, 12: 0.0285, 13: 0.0023, 14: 0.0459,
                             15: 0.3291, 16: 0.3105, 17: 0.0938, 18: 0.0128, 19: 0.0035, 20: 0.0003, "pmin": 0.0003},
                "D13S317": {6: 0.0003, 7:  0.0003, 8: 0.1489, 9: 0.0928, 10: 0.0777, 11: 0.3376, 11.1: 0.0003,
                            12: 0.2224, 13: 0.0839, 14: 0.0352, 15: 0.0012, "pmin": 0.0003},
                "D5S818": {7: 0.0088, 8: 0.0014, 9: 0.0489, 10: 0.0850, 11: 0.3294, 12: 0.3528, 13: 0.1477, 14: 0.0113,
                           15: 0.0013, "pmin": 0.0003}}
    codis_old = ["CSF1PO", "D3S1358", "D5S818", "D7S820", "D8S1179", "D13S317", "D16S539", "D18S51", "D21S11", "FGA",
                 "TH01", "TPOX", "vWA"]
    codis_new = ["CSF1PO", "D3S1358", "D5S818", "D7S820", "D8S1179", "D13S317", "D16S539", "D18S51", "D21S11", "FGA",
                 "TH01", "TPOX", "vWA", "D1S1656", "D2S441", "D2S1338", "D10S1248", "D12S391", "D19S433", "D22S1045"]

    # Upload data of trio and potential fathers
    trio_df = pd.read_excel("yak_v12_trio.xlsx")    # Enter the file name containing trio's data
    pf_df = pd.read_excel("yak_v12_pf.xlsx")    # Enter the file name with potential fathers' data
    columns_list = trio_df.columns.values.tolist()
    output_df = pd.DataFrame(columns=columns_list, index=pf_df.index)
    n = 0
    while n != len(output_df):
        for i in range(0, len(trio_df.index), 5):
            real_mother_df = trio_df.iloc[[i + 1]]    # real mother's data
            for ch in range(3):
                child_df = trio_df.iloc[[i + 2 + ch]]    # df with 1st, 2nd and 3rd child's data
                temp_pf_df = pf_df[(pf_df['Child_ID'] == child_df.iloc[0]["№"])]    # temporal df with all PF data
                for ii in range(0, len(temp_pf_df.index)):
                    one_pf_df = temp_pf_df.iloc[[ii]]    # df with one PF data
                    name = one_pf_df.iloc[0]['№'].split("_")
                    q_old_codis_trio_ref = []
                    q_new_codis_trio_ref = []
                    q_old_codis_duo_ref = []
                    q_new_codis_duo_ref = []
                    q_old_codis_trio_rus = []
                    q_new_codis_trio_rus = []
                    q_old_codis_duo_rus = []
                    q_new_codis_duo_rus = []
                    # OLD CODIS and TRIO
                    if ("old" in name) and ("trio" in name):
                        for loc in range(1, len(columns_list) - 27, 2):
                            el = columns_list[loc].split('_')[0]
                            if el in codis_old:
                                pf = [one_pf_df.iloc[0][columns_list[loc]], one_pf_df.iloc[0][columns_list[loc + 1]]]
                                rm = [real_mother_df.iloc[0][columns_list[loc]],
                                      real_mother_df.iloc[0][columns_list[loc + 1]]]
                                rc = [child_df.iloc[0][columns_list[loc]], child_df.iloc[0][columns_list[loc + 1]]]
                                pf.sort()
                                rm.sort()
                                rc.sort()
                                # check if there is mutation:
                                lack_of_mut = father_trio(pf[0], pf[1], rm[0], rm[1], rc[0], rc[1])
                                if lack_of_mut:
                                    f_allele_trio, other_allele_trio, knowledge_trio = \
                                        trio_freq_father_allele(rm[0], rm[1], rc[0], rc[1])
                                    q_old_trio_ref = calculate_q(codis_old, el, ref_dict, f_allele_trio,
                                                                 other_allele_trio, knowledge_trio)
                                    q_old_trio_rus = calculate_q_rus(codis_old, el, rus_dict, f_allele_trio,
                                                                     other_allele_trio, knowledge_trio)
                                else:
                                    q_old_trio_ref = None
                                    q_old_trio_rus = None
                                if q_old_trio_ref:
                                    q_old_codis_trio_ref.append(q_old_trio_ref)
                                if q_old_trio_rus:
                                    q_old_codis_trio_rus.append(q_old_trio_rus)
                        multiplication_old_trio_ref = 1
                        multiplication_old_trio_rus = 1
                        for m in q_old_codis_trio_ref:
                            multiplication_old_trio_ref *= m
                        for m in q_old_codis_trio_rus:
                            multiplication_old_trio_rus *= m
                        output_df.iloc[n]["p_father_PI_ref"] = 1 / multiplication_old_trio_ref
                        output_df.iloc[n]["p_father_PP_ref"] = 1 / (1 + multiplication_old_trio_ref)
                        output_df.iloc[n]["p_father_PI_rus"] = 1 / multiplication_old_trio_rus
                        output_df.iloc[n]["p_father_PP_rus"] = 1 / (1 + multiplication_old_trio_rus)
                        output_df.iloc[n]["Child_ID"] = child_df.iloc[0]["№"]
                        n += 1
                    # OLD CODIS and DUO
                    if "old" in name and "duo" in name:
                        for loc in range(1, len(columns_list) - 27, 2):
                            el = columns_list[loc].split('_')[0]
                            if el in codis_old:
                                pf = [one_pf_df.iloc[0][columns_list[loc]], one_pf_df.iloc[0][columns_list[loc + 1]]]
                                rc = [child_df.iloc[0][columns_list[loc]], child_df.iloc[0][columns_list[loc + 1]]]
                                pf.sort()
                                rc.sort()
                                # check if there is mutation:
                                lack_of_mut = father_duo(pf[0], pf[1], rc[0], rc[1])
                                if lack_of_mut:
                                    f_allele_duo, other_allele_duo, knowledge_duo = duo_freq_father_allele(rc[0], rc[1])
                                    q_old_duo_ref = calculate_q(codis_old, el, ref_dict, f_allele_duo, other_allele_duo,
                                                                knowledge_duo)
                                    q_old_duo_rus = calculate_q_rus(codis_old, el, rus_dict, f_allele_duo,
                                                                    other_allele_duo, knowledge_duo)
                                else:
                                    q_old_duo_ref = None
                                    q_old_duo_rus = None
                                if q_old_duo_ref:
                                    q_old_codis_duo_ref.append(q_old_duo_ref)
                                if q_old_duo_rus:
                                    q_old_codis_duo_rus.append(q_old_duo_rus)
                        multiplication_old_duo_ref = 1
                        multiplication_old_duo_rus = 1
                        for m in q_old_codis_duo_ref:
                            multiplication_old_duo_ref *= m
                        for m in q_old_codis_duo_rus:
                            multiplication_old_duo_rus *= m
                        output_df.iloc[n]["p_father_PI_ref"] = 1 / multiplication_old_duo_ref
                        output_df.iloc[n]["p_father_PP_ref"] = 1 / (1 + multiplication_old_duo_ref)
                        output_df.iloc[n]["p_father_PI_rus"] = 1 / multiplication_old_duo_rus
                        output_df.iloc[n]["p_father_PP_rus"] = 1 / (1 + multiplication_old_duo_rus)
                        output_df.iloc[n]["Child_ID"] = child_df.iloc[0]["№"]
                        n += 1
                    # NEW CODIS and TRIO
                    if "new" in name and "trio" in name:
                        for loc in range(1, len(columns_list) - 27, 2):
                            el = columns_list[loc].split('_')[0]
                            if el in codis_new:
                                pf = [one_pf_df.iloc[0][columns_list[loc]], one_pf_df.iloc[0][columns_list[loc + 1]]]
                                rm = [real_mother_df.iloc[0][columns_list[loc]],
                                      real_mother_df.iloc[0][columns_list[loc + 1]]]
                                rc = [child_df.iloc[0][columns_list[loc]], child_df.iloc[0][columns_list[loc + 1]]]
                                pf.sort()
                                rm.sort()
                                rc.sort()
                                # check if there is mutation:
                                lack_of_mut = father_trio(pf[0], pf[1], rm[0], rm[1], rc[0], rc[1])
                                if lack_of_mut:
                                    f_allele_trio, other_allele_trio, knowledge_trio = \
                                        trio_freq_father_allele(rm[0], rm[1], rc[0], rc[1])
                                    q_new_trio_ref = calculate_q(codis_new, el, ref_dict, f_allele_trio,
                                                                 other_allele_trio, knowledge_trio)
                                    q_new_trio_rus = calculate_q_rus(codis_new, el, rus_dict, f_allele_trio,
                                                                     other_allele_trio, knowledge_trio)
                                else:
                                    q_new_trio_ref = None
                                    q_new_trio_rus = None
                                if q_new_trio_ref:
                                    q_new_codis_trio_ref.append(q_new_trio_ref)
                                if q_new_trio_rus:
                                    q_new_codis_trio_rus.append(q_new_trio_rus)
                        multiplication_new_trio_ref = 1
                        multiplication_new_trio_rus = 1
                        for m in q_new_codis_trio_ref:
                            multiplication_new_trio_ref *= m
                        for m in q_new_codis_trio_rus:
                            multiplication_new_trio_rus *= m
                        output_df.iloc[n]["p_father_PI_ref"] = 1 / multiplication_new_trio_ref
                        output_df.iloc[n]["p_father_PP_ref"] = 1 / (1 + multiplication_new_trio_ref)
                        output_df.iloc[n]["p_father_PI_rus"] = 1 / multiplication_new_trio_rus
                        output_df.iloc[n]["p_father_PP_rus"] = 1 / (1 + multiplication_new_trio_rus)
                        output_df.iloc[n]["Child_ID"] = child_df.iloc[0]["№"]
                        n += 1
                    # NEW CODIS and DUO
                    if "new" in name and "duo" in name:
                        for loc in range(1, len(columns_list) - 27, 2):
                            el = columns_list[loc].split('_')[0]
                            if el in codis_new:
                                pf = [one_pf_df.iloc[0][columns_list[loc]], one_pf_df.iloc[0][columns_list[loc + 1]]]
                                rc = [child_df.iloc[0][columns_list[loc]], child_df.iloc[0][columns_list[loc + 1]]]
                                pf.sort()
                                rc.sort()
                                # check if there is mutation:
                                lack_of_mut = father_duo(pf[0], pf[1], rc[0], rc[1])
                                if lack_of_mut:
                                    f_allele_duo, other_allele_duo, knowledge_duo = duo_freq_father_allele(rc[0], rc[1])
                                    q_new_duo_ref = calculate_q(codis_new, el, ref_dict, f_allele_duo, other_allele_duo,
                                                                knowledge_duo)
                                    q_new_duo_rus = calculate_q_rus(codis_new, el, rus_dict, f_allele_duo,
                                                                    other_allele_duo, knowledge_duo)
                                else:
                                    q_new_duo_ref = None
                                    q_new_duo_rus = None
                                if q_new_duo_ref:
                                    q_new_codis_duo_ref.append(q_new_duo_ref)
                                if q_new_duo_rus:
                                    q_new_codis_duo_rus.append(q_new_duo_rus)
                        multiplication_new_duo_ref = 1
                        multiplication_new_duo_rus = 1
                        for m in q_new_codis_duo_ref:
                            multiplication_new_duo_ref *= m
                        for m in q_new_codis_duo_rus:
                            multiplication_new_duo_rus *= m
                        output_df.iloc[n]["p_father_PI_ref"] = 1 / multiplication_new_duo_ref
                        output_df.iloc[n]["p_father_PP_ref"] = 1 / (1 + multiplication_new_duo_ref)
                        output_df.iloc[n]["p_father_PI_rus"] = 1 / multiplication_new_duo_rus
                        output_df.iloc[n]["p_father_PP_rus"] = 1 / (1 + multiplication_new_duo_rus)
                        output_df.iloc[n]["Child_ID"] = child_df.iloc[0]["№"]
                        n += 1
    output_df.to_excel("YAK_output_PI_PP_pf.xlsx", index=False)
    print(round(time.time() - start, 2), 's')


if __name__ == "__main__":
    main()
