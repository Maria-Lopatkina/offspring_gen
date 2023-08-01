import pandas as pd
import time


def empty_freq_table():    # Create empty dictionaries for allele frequencies
    freq_df = pd.read_excel("ASTR22_main.xlsx")
    f_columns = freq_df.columns.values.tolist()
    key_dict = []
    for el in range(1, len(f_columns) - 2, 2):
        key_dict.append(list(f_columns[el].split("_"))[0])
    freq_dict = dict.fromkeys(key_dict)
    return freq_dict


def calculate_frequencies(p, d):    # Calculate alleles frequencies for three populations
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


def hyp_trio_wo_mut(m1, m2, k1, k2, our_dict, rus_dict, str_elem):
    #  1 and 2 - for our data, 3 and 4 - for rus
    p1 = p3 = 1
    p2 = calculate_p_trio(m1, m2, k1, k2, our_dict, str_elem)
    p4 = calculate_p_trio(m1, m2, k1, k2, rus_dict, str_elem)
    return p1, p2, p3, p4


def hyp_duo_wo_mut(f1, f2, k1, k2, our_dict, rus_dict, str_elem):
    #  1 and 2 - for our data, 3 and 4 - for rus
    p1, p2 = calculate_p_wo_mut_duo(f1, f2, k1, k2, our_dict, str_elem)
    p3, p4 = calculate_p_wo_mut_duo(f1, f2, k1, k2, rus_dict, str_elem)
    return p1, p2, p3, p4


def calculate_p_trio(m1, m2, k1, k2, dic, allele):
    p = None
    s_al = None
    if k1 == k2:
        if k1 in dic[allele]:
            p = dic[allele][k1] * (2 - dic[allele][k1])
        else:
            p = dic[allele]["pmin"] * (2 - dic[allele]["pmin"])
    if k1 != k2:
        if k1 == m1 and k2 == m2:
            if k1 in dic[allele] and k2 in dic[allele]:
                p = (dic[allele][k1] + dic[allele][k2]) * (2 - (dic[allele][k1] + dic[allele][k2]))
            elif k1 not in dic[allele] and k2 in dic[allele]:
                p = (dic[allele]["pmin"] + dic[allele][k2]) * (2 - (dic[allele]["pmin"] + dic[allele][k2]))
            elif k1 in dic[allele] and k2 not in dic[allele]:
                p = (dic[allele][k1] + dic[allele]["pmin"]) * (2 - (dic[allele][k1] + dic[allele]["pmin"]))
            elif k1 not in dic[allele] and k2 not in dic[allele]:
                p = (dic[allele]["pmin"] + dic[allele]["pmin"]) * (2 - (dic[allele]["pmin"] + dic[allele]["pmin"]))
        else:
            if k1 == m1 or k1 == m2:
                s_al = k2
            if k2 == m1 or k2 == m2:
                s_al = k1
            if s_al in dic[allele]:
                p = dic[allele][s_al] * (2 - dic[allele][s_al])
            else:
                p = dic[allele]["pmin"] * (2 - dic[allele]["pmin"])
    return p


def calculate_p_wo_mut_duo(f1, f2, k1, k2, dic, allele):
    p1 = p2 = None
    s_al = None
    if k1 == k2:    # Child is homozygous
        if k1 in dic[allele]:
            p1 = dic[allele][k1] * (2 - dic[allele][k1])
            p2 = (dic[allele][k1] * (2 - dic[allele][k1])) ** 2
        else:
            p1 = dic[allele]["pmin"] * (2 - dic[allele]["pmin"])
            p2 = (dic[allele]["pmin"] * (2 - dic[allele]["pmin"])) ** 2
    else:
        if k1 == f1 and k2 == f2:
            if k1 in dic[allele] and k2 in dic[allele]:
                p1 = (dic[allele][k1] + dic[allele][k2]) * (2 - (dic[allele][k1] + dic[allele][k2]))
                p2 = 2 * dic[allele][k1] * (2 - dic[allele][k1]) * dic[allele][k2] * (2 - dic[allele][k2]) - \
                     (2 * dic[allele][k1] * dic[allele][k2]) ** 2
            elif k1 not in dic[allele] and k2 in dic[allele]:
                p1 = (dic[allele]["pmin"] + dic[allele][k2]) * (2 - (dic[allele]["pmin"] + dic[allele][k2]))
                p2 = 2 * dic[allele]["pmin"] * (2 - dic[allele]["pmin"]) * dic[allele][k2] * (2 - dic[allele][k2]) - \
                     (2 * dic[allele]["pmin"] * dic[allele][k2]) ** 2
            elif k1 in dic[allele] and k2 not in dic[allele]:
                p1 = (dic[allele][k1] + dic[allele]["pmin"]) * (2 - (dic[allele][k1] + dic[allele]["pmin"]))
                p2 = 2 * dic[allele][k1] * (2 - dic[allele][k1]) * dic[allele]["pmin"] * (2 - dic[allele]["pmin"]) - \
                     (2 * dic[allele][k1] * dic[allele]["pmin"]) ** 2
            elif k1 not in dic[allele] and k2 not in dic[allele]:
                p1 = (dic[allele]["pmin"] + dic[allele]["pmin"]) * (2 - (dic[allele]["pmin"] + dic[allele]["pmin"]))
                p2 = 2 * dic[allele]["pmin"] * (2 - dic[allele]["pmin"]) * dic[allele]["pmin"] * \
                     (2 - dic[allele]["pmin"]) - (2 * dic[allele]["pmin"] * dic[allele]["pmin"]) ** 2
        else:
            for i in k1, k2:
                for j in f1, f2:
                    if i == j:
                        if k1 == i:
                            s_al = k2
                        elif k2 == i:
                            s_al = k1
            if s_al in dic[allele]:
                p1 = dic[allele][s_al] * (2 - dic[allele][s_al])
            else:
                p1 = dic[allele]["pmin"] * (2 - dic[allele]["pmin"])
            if k1 in dic[allele] and k2 in dic[allele]:
                p2 = 2 * dic[allele][k1] * (2 - dic[allele][k1]) * dic[allele][k2] * (2 - dic[allele][k2]) - \
                     (2 * dic[allele][k1] * dic[allele][k2]) ** 2
            elif k1 not in dic[allele] and k2 in dic[allele]:
                p2 = 2 * dic[allele]["pmin"] * (2 - dic[allele]["pmin"]) * dic[allele][k2] * (2 - dic[allele][k2]) - \
                     (2 * dic[allele]["pmin"] * dic[allele][k2]) ** 2
            elif k1 in dic[allele] and k2 not in dic[allele]:
                p2 = 2 * dic[allele][k1] * (2 - dic[allele][k1]) * dic[allele]["pmin"] * (2 - dic[allele]["pmin"]) - \
                     (2 * dic[allele][k1] * dic[allele]["pmin"]) ** 2
            elif k1 not in dic[allele] and k2 not in dic[allele]:
                p2 = 2 * dic[allele]["pmin"] * (2 - dic[allele]["pmin"]) * dic[allele]["pmin"] * \
                     (2 - dic[allele]["pmin"]) - (2 * dic[allele]["pmin"] * dic[allele]["pmin"]) ** 2
    return p1, p2


def main():
    start = time.time()
    # Create dictionary with alleles frequencies
    population = input("Enter the name of the population: ")
    our_dict = empty_freq_table()
    our_dict = calculate_frequencies(population, our_dict)
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
    columns_list = trio_df.columns.values.tolist()
    columns_list.append("LR_ref_trio_old_codis_f")
    columns_list.append("LR_rus_trio_old_codis_f")
    columns_list.append("LR_ref_trio_new_codis_f")
    columns_list.append("LR_rus_trio_new_codis_f")
    columns_list.append("LR_ref_duo_old_codis_f")
    columns_list.append("LR_rus_duo_old_codis_f")
    columns_list.append("LR_ref_duo_new_codis_f")
    columns_list.append("LR_rus_duo_new_codis_f")
    columns_list.append("LR_ref_trio_old_codis_m")
    columns_list.append("LR_rus_trio_old_codis_m")
    columns_list.append("LR_ref_trio_new_codis_m")
    columns_list.append("LR_rus_trio_new_codis_m")
    columns_list.append("LR_ref_duo_old_codis_m")
    columns_list.append("LR_rus_duo_old_codis_m")
    columns_list.append("LR_ref_duo_new_codis_m")
    columns_list.append("LR_rus_duo_new_codis_m")
    output_df = pd.DataFrame(columns=columns_list, index=trio_df.index)
    new_index = 0
    for i in range(0, len(trio_df.index), 5):
        real_father_df = trio_df.iloc[[i]]    # real father data
        real_mother_df = trio_df.iloc[[i + 1]]    # real mother data
        new_index += 2
        for ch in range(3):
            child_df = trio_df.iloc[[i + 2 + ch]]    # df with 1st, 2nd and 3rd child's data
            hyp1_ref_trio_old_f = []
            hyp2_ref_trio_old_f = []
            hyp1_rus_trio_old_f = []
            hyp2_rus_trio_old_f = []
            hyp1_ref_trio_new_f = []
            hyp2_ref_trio_new_f = []
            hyp1_rus_trio_new_f = []
            hyp2_rus_trio_new_f = []
            hyp1_ref_duo_old_f = []
            hyp2_ref_duo_old_f = []
            hyp1_rus_duo_old_f = []
            hyp2_rus_duo_old_f = []
            hyp1_ref_duo_new_f = []
            hyp2_ref_duo_new_f = []
            hyp1_rus_duo_new_f = []
            hyp2_rus_duo_new_f = []
            hyp1_ref_trio_old_m = []
            hyp2_ref_trio_old_m = []
            hyp1_rus_trio_old_m = []
            hyp2_rus_trio_old_m = []
            hyp1_ref_trio_new_m = []
            hyp2_ref_trio_new_m = []
            hyp1_rus_trio_new_m = []
            hyp2_rus_trio_new_m = []
            hyp1_ref_duo_old_m = []
            hyp2_ref_duo_old_m = []
            hyp1_rus_duo_old_m = []
            hyp2_rus_duo_old_m = []
            hyp1_ref_duo_new_m = []
            hyp2_ref_duo_new_m = []
            hyp1_rus_duo_new_m = []
            hyp2_rus_duo_new_m = []
            for loc in range(1, len(columns_list) - 43, 2):
                el = columns_list[loc].split('_')[0]
                rf = [real_father_df.iloc[0][columns_list[loc]], real_father_df.iloc[0][columns_list[loc + 1]]]
                rm = [real_mother_df.iloc[0][columns_list[loc]], real_mother_df.iloc[0][columns_list[loc + 1]]]
                rc = [child_df.iloc[0][columns_list[loc]], child_df.iloc[0][columns_list[loc + 1]]]
                rf.sort()
                rm.sort()
                rc.sort()
                if el in codis_old:
                    p1, p2, p3, p4 = hyp_trio_wo_mut(rm[0], rm[1], rc[0], rc[1], our_dict, rus_dict, el)
                    hyp1_ref_trio_old_f.append(p1)
                    hyp2_ref_trio_old_f.append(p2)
                    hyp1_rus_trio_old_f.append(p3)
                    hyp2_rus_trio_old_f.append(p4)
                    p1, p2, p3, p4 = hyp_trio_wo_mut(rf[0], rf[1], rc[0], rc[1], our_dict, rus_dict, el)
                    hyp1_ref_trio_old_m.append(p1)
                    hyp2_ref_trio_old_m.append(p2)
                    hyp1_rus_trio_old_m.append(p3)
                    hyp2_rus_trio_old_m.append(p4)
                    p1, p2, p3, p4 = hyp_duo_wo_mut(rf[0], rf[1], rc[0], rc[1], our_dict, rus_dict, el)
                    hyp1_ref_duo_old_f.append(p1)
                    hyp2_ref_duo_old_f.append(p2)
                    hyp1_rus_duo_old_f.append(p3)
                    hyp2_rus_duo_old_f.append(p4)
                    p1, p2, p3, p4 = hyp_duo_wo_mut(rm[0], rm[1], rc[0], rc[1], our_dict, rus_dict, el)
                    hyp1_ref_duo_old_m.append(p1)
                    hyp2_ref_duo_old_m.append(p2)
                    hyp1_rus_duo_old_m.append(p3)
                    hyp2_rus_duo_old_m.append(p4)
                if el in codis_new:
                    p1, p2, p3, p4 = hyp_trio_wo_mut(rm[0], rm[1], rc[0], rc[1], our_dict, rus_dict, el)
                    hyp1_ref_trio_new_f.append(p1)
                    hyp2_ref_trio_new_f.append(p2)
                    hyp1_rus_trio_new_f.append(p3)
                    hyp2_rus_trio_new_f.append(p4)
                    p1, p2, p3, p4 = hyp_trio_wo_mut(rf[0], rf[1], rc[0], rc[1], our_dict, rus_dict, el)
                    hyp1_ref_trio_new_m.append(p1)
                    hyp2_ref_trio_new_m.append(p2)
                    hyp1_rus_trio_new_m.append(p3)
                    hyp2_rus_trio_new_m.append(p4)
                    p1, p2, p3, p4 = hyp_duo_wo_mut(rf[0], rf[1], rc[0], rc[1], our_dict, rus_dict, el)
                    hyp1_ref_duo_new_f.append(p1)
                    hyp2_ref_duo_new_f.append(p2)
                    hyp1_rus_duo_new_f.append(p3)
                    hyp2_rus_duo_new_f.append(p4)
                    p1, p2, p3, p4 = hyp_duo_wo_mut(rm[0], rm[1], rc[0], rc[1], our_dict, rus_dict, el)
                    hyp1_ref_duo_new_m.append(p1)
                    hyp2_ref_duo_new_m.append(p2)
                    hyp1_rus_duo_new_m.append(p3)
                    hyp2_rus_duo_new_m.append(p4)
            multiplication_1 = 1
            multiplication_2 = 1
            multiplication_3 = 1
            multiplication_4 = 1
            multiplication_5 = 1
            multiplication_6 = 1
            multiplication_7 = 1
            multiplication_8 = 1
            multiplication_9 = 1
            multiplication_10 = 1
            multiplication_11 = 1
            multiplication_12 = 1
            multiplication_13 = 1
            multiplication_14 = 1
            multiplication_15 = 1
            multiplication_16 = 1
            multiplication_17 = 1
            multiplication_18 = 1
            multiplication_19 = 1
            multiplication_20 = 1
            multiplication_21 = 1
            multiplication_22 = 1
            multiplication_23 = 1
            multiplication_24 = 1
            multiplication_25 = 1
            multiplication_26 = 1
            multiplication_27 = 1
            multiplication_28 = 1
            multiplication_29 = 1
            multiplication_30 = 1
            multiplication_31 = 1
            multiplication_32 = 1
            for m in hyp1_ref_trio_old_f:
                multiplication_1 *= m
            for m in hyp2_ref_trio_old_f:
                multiplication_2 *= m
            for m in hyp1_rus_trio_old_f:
                multiplication_3 *= m
            for m in hyp2_rus_trio_old_f:
                multiplication_4 *= m
            for m in hyp1_ref_trio_new_f:
                multiplication_5 *= m
            for m in hyp2_ref_trio_new_f:
                multiplication_6 *= m
            for m in hyp1_rus_trio_new_f:
                multiplication_7 *= m
            for m in hyp2_rus_trio_new_f:
                multiplication_8 *= m
            for m in hyp1_ref_duo_old_f:
                multiplication_9 *= m
            for m in hyp2_ref_duo_old_f:
                multiplication_10 *= m
            for m in hyp1_rus_duo_old_f:
                multiplication_11 *= m
            for m in hyp2_rus_duo_old_f:
                multiplication_12 *= m
            for m in hyp1_ref_duo_new_f:
                multiplication_13 *= m
            for m in hyp2_ref_duo_new_f:
                multiplication_14 *= m
            for m in hyp1_rus_duo_new_f:
                multiplication_15 *= m
            for m in hyp2_rus_duo_new_f:
                multiplication_16 *= m
            for m in hyp1_ref_trio_old_m:
                multiplication_17 *= m
            for m in hyp2_ref_trio_old_m:
                multiplication_18 *= m
            for m in hyp1_rus_trio_old_m:
                multiplication_19 *= m
            for m in hyp2_rus_trio_old_m:
                multiplication_20 *= m
            for m in hyp1_ref_trio_new_m:
                multiplication_21 *= m
            for m in hyp2_ref_trio_new_m:
                multiplication_22 *= m
            for m in hyp1_rus_trio_new_m:
                multiplication_23 *= m
            for m in hyp2_rus_trio_new_m:
                multiplication_24 *= m
            for m in hyp1_ref_duo_old_m:
                multiplication_25 *= m
            for m in hyp2_ref_duo_old_m:
                multiplication_26 *= m
            for m in hyp1_rus_duo_old_m:
                multiplication_27 *= m
            for m in hyp2_rus_duo_old_m:
                multiplication_28 *= m
            for m in hyp1_ref_duo_new_m:
                multiplication_29 *= m
            for m in hyp2_ref_duo_new_m:
                multiplication_30 *= m
            for m in hyp1_rus_duo_new_m:
                multiplication_31 *= m
            for m in hyp2_rus_duo_new_m:
                multiplication_32 *= m
            output_df.iloc[new_index]["LR_ref_trio_old_codis_f"] = multiplication_1 / multiplication_2
            output_df.iloc[new_index]["LR_rus_trio_old_codis_f"] = multiplication_3 / multiplication_4
            output_df.iloc[new_index]["LR_ref_trio_new_codis_f"] = multiplication_5 / multiplication_6
            output_df.iloc[new_index]["LR_rus_trio_new_codis_f"] = multiplication_7 / multiplication_8
            output_df.iloc[new_index]["LR_ref_duo_old_codis_f"] = multiplication_9 / multiplication_10
            output_df.iloc[new_index]["LR_rus_duo_old_codis_f"] = multiplication_11 / multiplication_12
            output_df.iloc[new_index]["LR_ref_duo_new_codis_f"] = multiplication_13 / multiplication_14
            output_df.iloc[new_index]["LR_rus_duo_new_codis_f"] = multiplication_15 / multiplication_16
            output_df.iloc[new_index]["LR_ref_trio_old_codis_m"] = multiplication_17 / multiplication_18
            output_df.iloc[new_index]["LR_rus_trio_old_codis_m"] = multiplication_19 / multiplication_20
            output_df.iloc[new_index]["LR_ref_trio_new_codis_m"] = multiplication_21 / multiplication_22
            output_df.iloc[new_index]["LR_rus_trio_new_codis_m"] = multiplication_23 / multiplication_24
            output_df.iloc[new_index]["LR_ref_duo_old_codis_m"] = multiplication_25 / multiplication_26
            output_df.iloc[new_index]["LR_rus_duo_old_codis_m"] = multiplication_27 / multiplication_28
            output_df.iloc[new_index]["LR_ref_duo_new_codis_m"] = multiplication_29 / multiplication_30
            output_df.iloc[new_index]["LR_rus_duo_new_codis_m"] = multiplication_31 / multiplication_32
            new_index += 1
    output_df.to_excel("output_LR_yak_trio.xlsx", index=False)    # Enter the name of the output file
    print(round(time.time() - start, 2), 's')


if __name__ == "__main__":
    main()
