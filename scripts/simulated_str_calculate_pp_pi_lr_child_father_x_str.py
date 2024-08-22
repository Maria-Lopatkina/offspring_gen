import copy
import pandas as pd
import random
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

def inx_mutation(n_p, k, m_r, one, two, three):
    one_stp = []
    two_stp = []
    three_stp = []
    list_of_all_mut = []
    mutation = int(n_p * m_r * k * 22)
    for m in range(round(mutation * one)):
        a = random.choice(range(n_p * k))
        if a not in list_of_all_mut:
            if a not in one_stp:
                one_stp.append(a)
    for m in range(round(mutation * two)):
        a = random.choice(range(n_p * k))
        if a not in list_of_all_mut:
            if a not in two_stp:
                two_stp.append(a)
    for m in range(round(mutation * three)):
        a = random.choice(range(n_p * k))
        if a not in list_of_all_mut:
            if a not in three_stp:
                three_stp.append(a)
    return one_stp, two_stp, three_stp

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

def find_father(i_list, c_list, t_df, k_df, family, n_of_mis, codis, father_ind, mother_ind):
    n_fathers = 0
    ids_father = {}
    answer = None
    for i in i_list:
        if i != father_ind and i != mother_ind:
            count_of_mismatches = 0
            mismatch = []
            pot_father_df = t_df.loc[[i]]
            real_mother_df = t_df.loc[[mother_ind]]
            for loc in range(1, len(c_list) - 26, 2):
                el = c_list[loc].split('_')[0]
                if el in codis:
                    pf1, pf2 = pot_father_df.iloc[0][c_list[loc]], \
                        pot_father_df.iloc[0][c_list[loc + 1]]
                    rm1, rm2 = real_mother_df.iloc[0][c_list[loc]], \
                        real_mother_df.iloc[0][c_list[loc + 1]]
                    rk1, rk2 = k_df.iloc[0][c_list[loc]], \
                        k_df.iloc[0][c_list[loc + 1]]
                    if family == 3:
                        answer = father_trio(pf1, pf2, rm1, rm2, rk1, rk2)
                    elif family == 2:
                        answer = father_duo(pf1, pf2, rk1, rk2)
                    if not answer:
                        count_of_mismatches += 1
                        mismatch.append(el)
            if count_of_mismatches <= n_of_mis:
                n_fathers += 1
                ids_father[i] = mismatch
    return n_fathers, ids_father

def trio_freq_father_allele(m1, m2, k1, k2):
    other_all = None
    if k1 == k2:    # Child is homozygous
        f_all = k1
        kn = True
        return f_all, other_all, kn
    if k1 != k2:    # Child is heterozygous
        if k1 == m1 and k2 == m2:
            f_all = k1
            other_all = k2
            kn = False
            return f_all, other_all, kn
        elif (k1 == m1 and k2 != m2) or (k1 == m2 and k2 != m1) or (k2 == m1 and k1 != m2) or (k2 == m2 and k1 != m1):
            for i in k1, k2:
                for j in m1, m2:
                    if i != j and i not in [m1, m2]:
                        f_all = i
                        kn = True
                        return f_all, other_all, kn

def duo_freq_father_allele(k1, k2):
    other_all = None
    if k1 == k2:    # Child is homozygote
        f_all = k1
        kn = True
        return f_all, other_all, kn
    elif k1 != k2:    # Child is heterozygote
        f_all = k1
        other_all = k2
        kn = False
        return f_all, other_all, kn

def calculate_q(str_list, allele_name, dic, f_a, o_a, know):
    q = None
    if allele_name in str_list:
        if know:
            if f_a in dic[allele_name]:
                q = dic[allele_name][f_a] * (2 - dic[allele_name][f_a])
            else:
                q = dic[allele_name]["pmin"] * (2 - dic[allele_name]["pmin"])
        else:
            if f_a in dic[allele_name] and o_a in dic[allele_name]:
                q = (dic[allele_name][f_a] + dic[allele_name][o_a]) * (2 - (dic[allele_name][f_a] +
                                                                            dic[allele_name][o_a]))
            elif f_a not in dic[allele_name] and o_a in dic[allele_name]:
                q = (dic[allele_name]["pmin"] + dic[allele_name][o_a]) * (2 - (dic[allele_name]["pmin"]
                                                                               + dic[allele_name][o_a]))
            elif f_a in dic[allele_name] and o_a not in dic[allele_name]:
                q = (dic[allele_name][f_a] + dic[allele_name]["pmin"]) * (2 - (dic[allele_name][f_a] +
                                                                               dic[allele_name]["pmin"]))
            elif f_a not in dic[allele_name] and o_a not in dic[allele_name]:
                q = (dic[allele_name]["pmin"] + dic[allele_name]["pmin"]) * \
                    (2 - (dic[allele_name]["pmin"] + dic[allele_name]["pmin"]))
        return q
    else:
        return None

# functions for LR calculation
def hyp_trio_wo_mut(m1, m2, k1, k2, our_dict, rus_dict, str_elem):
    #  1 and 2 - for our data, 3 and 4 - for rus
    p1 = p3 = 1
    p2 = calculate_p_trio(m1, m2, k1, k2, our_dict, str_elem)
    p4 = calculate_p_trio(m1, m2, k1, k2, rus_dict, str_elem)
    return p1, p2, p3, p4

def hyp_trio_mut(m1, m2, k1, k2, our_dict, rus_dict, str_elem, mut):
    #  1 and 2 - for our data, 3 and 4 - for rus
    p1 = p3 = mut
    p2 = calculate_p_trio(m1, m2, k1, k2, our_dict, str_elem)
    p4 = calculate_p_trio(m1, m2, k1, k2, rus_dict, str_elem)
    return p1, p2, p3, p4

def hyp_duo_wo_mut(f1, f2, k1, k2, our_dict, rus_dict, str_elem):
    #  1 and 2 - for our data, 3 and 4 - for rus
    p1, p2 = calculate_p_wo_mut_duo(f1, f2, k1, k2, our_dict, str_elem)
    p3, p4 = calculate_p_wo_mut_duo(f1, f2, k1, k2, rus_dict, str_elem)
    return p1, p2, p3, p4

def hyp_duo_mut(k1, k2, our_dict, rus_dict, str_elem, mut):
    #  1 and 2 - for our data, 3 and 4 - for rus
    p1, p2 = calculate_p_mut_duo(k1, k2, our_dict, str_elem, mut)
    p3, p4 = calculate_p_mut_duo(k1, k2, rus_dict, str_elem, mut)
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
                p2 = (2 * dic[allele][k1] * (2 - dic[allele][k1]) * dic[allele][k2] * (2 - dic[allele][k2]) -
                      (2 * dic[allele][k1] * dic[allele][k2]) ** 2)
            elif k1 not in dic[allele] and k2 in dic[allele]:
                p1 = (dic[allele]["pmin"] + dic[allele][k2]) * (2 - (dic[allele]["pmin"] + dic[allele][k2]))
                p2 = (2 * dic[allele]["pmin"] * (2 - dic[allele]["pmin"]) * dic[allele][k2] * (2 - dic[allele][k2]) -
                      (2 * dic[allele]["pmin"] * dic[allele][k2]) ** 2)
            elif k1 in dic[allele] and k2 not in dic[allele]:
                p1 = (dic[allele][k1] + dic[allele]["pmin"]) * (2 - (dic[allele][k1] + dic[allele]["pmin"]))
                p2 = (2 * dic[allele][k1] * (2 - dic[allele][k1]) * dic[allele]["pmin"] * (2 - dic[allele]["pmin"]) -
                      (2 * dic[allele][k1] * dic[allele]["pmin"]) ** 2)
            elif k1 not in dic[allele] and k2 not in dic[allele]:
                p1 = (dic[allele]["pmin"] + dic[allele]["pmin"]) * (2 - (dic[allele]["pmin"] + dic[allele]["pmin"]))
                p2 = (2 * dic[allele]["pmin"] * (2 - dic[allele]["pmin"]) * dic[allele]["pmin"] *
                      (2 - dic[allele]["pmin"]) - (2 * dic[allele]["pmin"] * dic[allele]["pmin"]) ** 2)
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
                p2 = (2 * dic[allele]["pmin"] * (2 - dic[allele]["pmin"]) * dic[allele]["pmin"] *
                      (2 - dic[allele]["pmin"]) - (2 * dic[allele]["pmin"] * dic[allele]["pmin"]) ** 2)
    return p1, p2

def calculate_p_mut_duo(k1, k2, dic, allele, mut):
    p1 = p2 = None
    if k1 == k2:
        if k1 in dic[allele]:
            p1 = dic[allele][k1] * (2 - dic[allele][k1]) * mut
            p2 = (dic[allele][k1] * (2 - dic[allele][k1])) ** 2
        else:
            p1 = dic[allele]["pmin"] * (2 - dic[allele]["pmin"]) * mut
            p2 = (dic[allele]["pmin"] * (2 - dic[allele]["pmin"])) ** 2
    else:
        if k1 in dic[allele] and k2 in dic[allele]:
            p1 = (dic[allele][k1] + dic[allele][k2]) * (2 - (dic[allele][k1] + dic[allele][k2])) * mut
            p2 = (2 * dic[allele][k1] * (2 - dic[allele][k1]) * dic[allele][k2] * (2 - dic[allele][k2]) -
                  (2 * dic[allele][k1] * dic[allele][k2]) ** 2)
        elif k1 not in dic[allele] and k2 in dic[allele]:
            p1 = (dic[allele]["pmin"] + dic[allele][k2]) * (2 - (dic[allele]["pmin"] + dic[allele][k2])) * mut
            p2 = (2 * dic[allele]["pmin"] * (2 - dic[allele]["pmin"]) * dic[allele][k2] * (2 - dic[allele][k2]) -
                  (2 * dic[allele]["pmin"] * dic[allele][k2]) ** 2)
        elif k1 in dic[allele] and k2 not in dic[allele]:
            p1 = (dic[allele][k1] + dic[allele]["pmin"]) * (2 - (dic[allele][k1] + dic[allele]["pmin"])) * mut
            p2 = (2 * dic[allele][k1] * (2 - dic[allele][k1]) * dic[allele]["pmin"] * (2 - dic[allele]["pmin"]) -
                  (2 * dic[allele][k1] * dic[allele]["pmin"]) ** 2)
        elif k1 not in dic[allele] and k2 not in dic[allele]:
            p1 = (dic[allele]["pmin"] + dic[allele]["pmin"]) * (2 - (dic[allele]["pmin"] + dic[allele]["pmin"])) * mut
            p2 = (2 * dic[allele]["pmin"] * (2 - dic[allele]["pmin"]) * dic[allele]["pmin"] * (2 - dic[allele]["pmin"])
                  - (2 * dic[allele]["pmin"] * dic[allele]["pmin"]) ** 2)
    return p1, p2

def count_step_trio(f1, f2, m1, m2, k1, k2):
    if k1 == m1 and k2 == m2:
        d = []
        for i in f1, f2:
            for j in k1, k2:
                d.append(abs(i - j))
        d.sort()
        if d[0] == 1:
            return 1
        else:
            return 0
    else:
        for i in m1, m2:
            for j in k1, k2:
                if i == j:
                    if i == k2:
                        f_all = k1
                    else:
                        f_all = k2
                    for f in f1, f2:
                        if abs(f - f_all) == 1:
                            return 1
                        else:
                            return 0

def count_step_duo(f1, f2, k1, k2):
    for i in f1, f2:
        for j in k1, k2:
            if abs(i - j) == 1:
                return 1
    return 0

def main():
    with open("../config_file.yaml", "r") as yaml_file:
        config = yaml.load(yaml_file, Loader=yaml.FullLoader)
    # Create dataframe with offsprings
    parents_df = pd.read_excel(config["main_table_path"])
    columns_list = parents_df.columns.values.tolist()
    columns_list = columns_list[:-1]
    columns_list.append("Child_ID")
    columns_list.append("Status")
    columns_list.append("PI_old_CODIS_trio_ref")
    columns_list.append("PP_old_CODIS_trio_ref")
    columns_list.append("PI_old_CODIS_duo_ref")
    columns_list.append("PP_old_CODIS_duo_ref")
    columns_list.append("PI_new_CODIS_trio_ref")
    columns_list.append("PP_new_CODIS_trio_ref")
    columns_list.append("PI_new_CODIS_duo_ref")
    columns_list.append("PP_new_CODIS_duo_ref")
    # Columns for russian reference dictionary
    columns_list.append("PI_old_CODIS_trio_rus")
    columns_list.append("PP_old_CODIS_trio_rus")
    columns_list.append("PI_old_CODIS_duo_rus")
    columns_list.append("PP_old_CODIS_duo_rus")
    columns_list.append("PI_new_CODIS_trio_rus")
    columns_list.append("PP_new_CODIS_trio_rus")
    columns_list.append("PI_new_CODIS_duo_rus")
    columns_list.append("PP_new_CODIS_duo_rus")
    # Columns for other locus from article about Chinese population
    columns_list.append("PI_CODIS_15_trio")
    # Add columns for the number of potential fathers
    columns_list.append("number_of_p_fathers_old_CODIS_trio")
    columns_list.append("number_of_p_fathers_old_CODIS_duo")
    columns_list.append("number_of_p_fathers_new_CODIS_trio")
    columns_list.append("number_of_p_fathers_new_CODIS_duo")
    columns_list.append("number_of_mutations")
    columns_list.append("step_mutation")
    # Add columns for LR calculations (for real parents)
    columns_list.append("LR_ref_trio_old_codis")
    columns_list.append("LR_rus_trio_old_codis")
    columns_list.append("LR_ref_trio_new_codis")
    columns_list.append("LR_rus_trio_new_codis")
    columns_list.append("LR_ref_duo_old_codis")
    columns_list.append("LR_rus_duo_old_codis")
    columns_list.append("LR_ref_duo_new_codis")
    columns_list.append("LR_rus_duo_new_codis")
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
        child_counter = 0
        # Get the indexes of kids with mutation:
        # inx_mut_1, inx_mut_2, inx_mut_3 = inx_mutation(config["number_of_pairs"], config["kids"], config["mut_rate"],
        #                                                config["one_step"], config["two_step"], config["three_step"])
        step_f = 0
        n = copy.deepcopy(config["number_of_pairs"])
        while n != 0:
            father = random.choice(index_list_f)
            mother = random.choice(index_list_m)
            new_pair = [father, mother]
            pair = set(new_pair)
            if pair not in list_for_pairs and father != mother:
                list_for_pairs.append(pair)
                pair_df = temp_parents_df.loc[[father, mother]]    # create df for one random pair
                kids_df = pd.DataFrame(columns=columns_list, index=[i for i in range(0, config["kids"])])
                for k in range(config["kids"]):
                    kids_df.iloc[k]["Status"] = "kid"
                    kids_df.iloc[k]["№"] = (str(pair_df.iloc[0]['№']) + "-" + str(pair_df.iloc[1]['№']) + "-" +
                                            str(k + 1))
                    # offsprings_df = pd.concat([offsprings_df, pair_df])  # parents data in output
                    # random_str = None
                    q_old_codis_trio_ref = []
                    q_new_codis_trio_ref = []
                    q_old_codis_duo_ref = []
                    q_new_codis_duo_ref = []
                    q_old_codis_trio_rus = []
                    q_new_codis_trio_rus = []
                    q_old_codis_duo_rus = []
                    q_new_codis_duo_rus = []
                    q_codis_15_trio = []
                    count_of_mismatch_rf = 0
                    # if (child_counter in inx_mut_1) or (child_counter in inx_mut_2) or (child_counter in inx_mut_3):
                    #     random_str = random.choice(config["all_str"])
                    #     count_of_mismatch_rf += 1
                    for elem in range(1, len(columns_list) - 34, 2):
                        allele = columns_list[elem].split('_')[0]
                        f_alleles = [pair_df.iloc[0][columns_list[elem]], pair_df.iloc[0][columns_list[elem + 1]]]
                        f_allele_random = random.choice(f_alleles)
                        # if (allele == random_str) and (child_counter in inx_mut_1):
                        #     a = random.choice([0, 1])
                        #     if a == 0:
                        #         f_allele_random -= 1
                        #     else:
                        #         f_allele_random += 1
                        #     step_f = 1
                        # if (allele == random_str) and (child_counter in inx_mut_2):
                        #     a = random.choice([0, 1])
                        #     if a == 0:
                        #         f_allele_random -= 2
                        #     else:
                        #         f_allele_random += 2
                        #     step_f = 2
                        # if (allele == random_str) and (child_counter in inx_mut_3):
                        #     a = random.choice([0, 1])
                        #     if a == 0:
                        #         f_allele_random -= 3
                        #     else:
                        #         f_allele_random += 3
                        #     step_f = 3
                        m_alleles = [pair_df.iloc[1][columns_list[elem]], pair_df.iloc[1][columns_list[elem + 1]]]
                        m_allele_random = random.choice(m_alleles)
                        possible_alleles = [f_allele_random, m_allele_random]
                        kid_alleles = random.sample(possible_alleles, 2)
                        kid_alleles.sort()
                        kids_df.iloc[k][columns_list[elem]] = kid_alleles[0]
                        kids_df.iloc[k][columns_list[elem+1]] = kid_alleles[1]
                        f_alleles.sort()
                        m_alleles.sort()
                        # Check the father allele in a case of trio
                        f_allele_trio, other_allele_trio, knowledge_trio = trio_freq_father_allele(m_alleles[0],
                                                                                                   m_alleles[1],
                                                                                                   kid_alleles[0],
                                                                                                   kid_alleles[1])
                        # Choose the dictionary to calculate Q and Q_all:
                        q_old_trio_ref = calculate_q(config["codis_old"], allele, ref_dict, f_allele_trio,
                                                     other_allele_trio, knowledge_trio)
                        if q_old_trio_ref:
                            q_old_codis_trio_ref.append(q_old_trio_ref)
                        q_new_trio_ref = calculate_q(config["codis_new"], allele, ref_dict, f_allele_trio,
                                                     other_allele_trio, knowledge_trio)
                        if q_new_trio_ref:
                            q_new_codis_trio_ref.append(q_new_trio_ref)
                        q_trio_15 = calculate_q(config["codis_15"], allele, ref_dict, f_allele_trio, other_allele_trio,
                                                knowledge_trio)
                        if q_trio_15:
                            q_codis_15_trio.append(q_trio_15)
                        # Check the father allele in a case of duo
                        f_allele_duo, other_allele_duo, knowledge_duo = duo_freq_father_allele(kid_alleles[0],
                                                                                               kid_alleles[1])
                        # Choose the dictionary to calculate Q and Q_all:
                        q_old_duo_ref = calculate_q(config["codis_old"], allele, ref_dict, f_allele_duo,
                                                    other_allele_duo, knowledge_duo)
                        if q_old_duo_ref:
                            q_old_codis_duo_ref.append(q_old_duo_ref)
                        q_new_duo_ref = calculate_q(config["codis_new"], allele, ref_dict, f_allele_duo,
                                                    other_allele_duo, knowledge_duo)
                        if q_new_duo_ref:
                            q_new_codis_duo_ref.append(q_new_duo_ref)
                        # Choose the dictionary to calculate Q and Q_all for new ref in trio:
                        q_old_trio_rus = calculate_q(config["codis_old"], allele, config["rus_dict"], f_allele_trio,
                                                     other_allele_trio, knowledge_trio)
                        if q_old_trio_rus:
                            q_old_codis_trio_rus.append(q_old_trio_rus)
                        q_new_trio_rus = calculate_q(config["codis_new"], allele, config["rus_dict"], f_allele_trio,
                                                     other_allele_trio, knowledge_trio)
                        if q_new_trio_rus:
                            q_new_codis_trio_rus.append(q_new_trio_rus)
                        # Choose the dictionary to calculate Q and Q_all for new ref in duo:
                        q_old_duo_rus = calculate_q(config["codis_old"], allele, config["rus_dict"], f_allele_duo,
                                                    other_allele_duo, knowledge_duo)
                        if q_old_duo_rus:
                            q_old_codis_duo_rus.append(q_old_duo_rus)
                        q_new_duo_rus = calculate_q(config["codis_new"], allele, config["rus_dict"], f_allele_duo,
                                                    other_allele_duo, knowledge_duo)
                        if q_new_duo_rus:
                            q_new_codis_duo_rus.append(q_new_duo_rus)
                    multiplication_old_trio_ref = 1
                    multiplication_old_duo_ref = 1
                    multiplication_new_trio_ref = 1
                    multiplication_new_duo_ref = 1
                    multiplication_old_trio_rus = 1
                    multiplication_old_duo_rus = 1
                    multiplication_new_trio_rus = 1
                    multiplication_new_duo_rus = 1
                    multiplication_15_trio = 1
                    for m in q_old_codis_trio_ref:
                        multiplication_old_trio_ref *= m
                    for m in q_old_codis_duo_ref:
                        multiplication_old_duo_ref *= m
                    for m in q_new_codis_trio_ref:
                        multiplication_new_trio_ref *= m
                    for m in q_new_codis_duo_ref:
                        multiplication_new_duo_ref *= m
                    for m in q_old_codis_trio_rus:
                        multiplication_old_trio_rus *= m
                    for m in q_old_codis_duo_rus:
                        multiplication_old_duo_rus *= m
                    for m in q_new_codis_trio_rus:
                        multiplication_new_trio_rus *= m
                    for m in q_new_codis_duo_rus:
                        multiplication_new_duo_rus *= m
                    for m in q_codis_15_trio:
                        multiplication_15_trio *= m
                    new_pair_df.iloc[0]["PI_old_CODIS_trio_ref"] = 1 / multiplication_old_trio_ref
                    new_pair_df.iloc[0]["PP_old_CODIS_trio_ref"] = 1 / (1 + multiplication_old_trio_ref)
                    new_pair_df.iloc[0]["PI_old_CODIS_duo_ref"] = 1 / multiplication_old_duo_ref
                    new_pair_df.iloc[0]["PP_old_CODIS_duo_ref"] = 1 / (1 + multiplication_old_duo_ref)
                    new_pair_df.iloc[0]["PI_new_CODIS_trio_ref"] = 1 / multiplication_new_trio_ref
                    new_pair_df.iloc[0]["PP_new_CODIS_trio_ref"] = 1 / (1 + multiplication_new_trio_ref)
                    new_pair_df.iloc[0]["PI_new_CODIS_duo_ref"] = 1 / multiplication_new_duo_ref
                    new_pair_df.iloc[0]["PP_new_CODIS_duo_ref"] = 1 / (1 + multiplication_new_duo_ref)
                    new_pair_df.iloc[0]["PI_old_CODIS_trio_rus"] = 1 / multiplication_old_trio_rus
                    new_pair_df.iloc[0]["PP_old_CODIS_trio_rus"] = 1 / (1 + multiplication_old_trio_rus)
                    new_pair_df.iloc[0]["PI_old_CODIS_duo_rus"] = 1 / multiplication_old_duo_rus
                    new_pair_df.iloc[0]["PP_old_CODIS_duo_rus"] = 1 / (1 + multiplication_old_duo_rus)
                    new_pair_df.iloc[0]["PI_new_CODIS_trio_rus"] = 1 / multiplication_new_trio_rus
                    new_pair_df.iloc[0]["PP_new_CODIS_trio_rus"] = 1 / (1 + multiplication_new_trio_rus)
                    new_pair_df.iloc[0]["PI_new_CODIS_duo_rus"] = 1 / multiplication_new_duo_rus
                    new_pair_df.iloc[0]["PP_new_CODIS_duo_rus"] = 1 / (1 + multiplication_new_duo_rus)
                    new_pair_df.iloc[0]["PI_CODIS_15_trio"] = 1 / multiplication_15_trio
                    new_pair_df.iloc[0]["number_of_mutations"] = count_of_mismatch_rf
                    # if step_f != 0:
                    #     new_pair_df.iloc[0]["step_mutation"] = step_f
                    # else:
                    #     new_pair_df.iloc[0]["step_mutation"] = "-"
                    kids_df.iloc[k]["population"] = pair_df.iloc[0]["population"]
                    kids_df.iloc[k]["Child_ID"] = kids_df.iloc[k]["№"]
                    kids_df.iloc[k]["number_of_mutations"] = "-"
                    # Find potential fathers for every child:
                    kids_df.iloc[k]["number_of_p_fathers_old_CODIS_trio"], ids1 = find_father(index_list_f,
                                                                                              columns_list,
                                                                                              temp_parents_df,
                                                                                              kids_df.loc[[k]], 3, 2,
                                                                                              config["codis_old"],
                                                                                              father, mother)
                    kids_df.iloc[k]["number_of_p_fathers_old_CODIS_duo"], ids2 = find_father(index_list_f,
                                                                                             columns_list,
                                                                                             temp_parents_df,
                                                                                             kids_df.loc[[k]], 2, 2,
                                                                                             config["codis_old"],
                                                                                             father, mother)
                    kids_df.iloc[k]["number_of_p_fathers_new_CODIS_trio"], ids3 = find_father(index_list_f,
                                                                                              columns_list,
                                                                                              temp_parents_df,
                                                                                              kids_df.loc[[k]], 3, 2,
                                                                                              config["codis_new"],
                                                                                              father, mother)
                    kids_df.iloc[k]["number_of_p_fathers_new_CODIS_duo"], ids4 = find_father(index_list_f,
                                                                                             columns_list,
                                                                                             temp_parents_df,
                                                                                             kids_df.loc[[k]], 2, 2,
                                                                                             config["codis_new"],
                                                                                             father, mother)
                    for elem in range(len(columns_list) - 33):
                        new_pair_df.iloc[0][columns_list[elem]] = pair_df.iloc[0][columns_list[elem]]
                        new_pair_df.iloc[1][columns_list[elem]] = pair_df.iloc[1][columns_list[elem]]
                    new_pair_df.iloc[0]["Status"] = "real_father"
                    new_pair_df.iloc[1]["Status"] = "real_mother"
                    new_pair_df.iloc[0]["Child_ID"] = (str(pair_df.iloc[0]['№']) + "-" + str(pair_df.iloc[1]['№']) +
                                                       "-" + str(k + 1))
                    new_pair_df.iloc[1]["Child_ID"] = (str(pair_df.iloc[0]['№']) + "-" + str(pair_df.iloc[1]['№']) +
                                                       "-" + str(k + 1))
                    new_pair_df.iloc[1]["number_of_mutations"] = 0
                    offsprings_df = pd.concat([offsprings_df, new_pair_df], ignore_index=True)
                    offsprings_df = pd.concat([offsprings_df, kids_df.loc[[k]]], ignore_index=True)

                    #################################################
                    # CALCULATE PI AND PP FOR FALSE POSITIVE FATHERS
                    # 1 old CODIS trio
                    for ii in ids1:
                        # create df for a pair of real mother and p father
                        p_pair_df = temp_parents_df.loc[[ii, mother]]
                        pot_father_df = pd.DataFrame(columns=columns_list, index=[0])
                        pot_father_df.iloc[0]["№"] = temp_parents_df.loc[ii]["№"] + "_old_codis_trio"
                        q_old_codis_trio_ref = []
                        q_old_codis_trio_rus = []
                        for elem in range(1, len(columns_list) - 34, 2):
                            allele = columns_list[elem].split('_')[0]
                            pot_father_df.iloc[0][columns_list[elem]] = p_pair_df.iloc[0][columns_list[elem]]
                            pot_father_df.iloc[0][columns_list[elem + 1]] = p_pair_df.iloc[0][columns_list[elem + 1]]
                            if (allele not in ids1[ii]) and (allele in config["codis_old"]):
                                m_alleles = [p_pair_df.iloc[1][columns_list[elem]],
                                             p_pair_df.iloc[1][columns_list[elem + 1]]]
                                kid_alleles = [kids_df.iloc[k][columns_list[elem]],
                                               kids_df.iloc[k][columns_list[elem + 1]]]
                                m_alleles.sort()
                                kid_alleles.sort()
                                f_allele_trio, other_allele_trio, knowledge_trio = \
                                    trio_freq_father_allele(m_alleles[0], m_alleles[1], kid_alleles[0], kid_alleles[1])
                                q_old_trio_ref = calculate_q(config["codis_old"], allele, ref_dict, f_allele_trio,
                                                             other_allele_trio, knowledge_trio)
                                q_old_trio_rus = calculate_q(config["codis_old"], allele, config["rus_dict"],
                                                             f_allele_trio, other_allele_trio, knowledge_trio)
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
                        pot_father_df.iloc[0]["PI_old_CODIS_trio_ref"] = 1 / multiplication_old_trio_ref
                        pot_father_df.iloc[0]["PP_old_CODIS_trio_ref"] = 1 / (1 + multiplication_old_trio_ref)
                        pot_father_df.iloc[0]["PI_old_CODIS_trio_rus"] = 1 / multiplication_old_trio_rus
                        pot_father_df.iloc[0]["PP_old_CODIS_trio_rus"] = 1 / (1 + multiplication_old_trio_rus)
                        pot_father_df.iloc[0]["population"] = p_pair_df.iloc[0]["population"]
                        pot_father_df.iloc[0]["Child_ID"] = kids_df.iloc[k]["№"]
                        pot_father_df.iloc[0]["Status"] = "false_positive_father"
                        pot_father_df.iloc[0]["number_of_mutations"] = len(ids1[ii])
                        offsprings_df = pd.concat([offsprings_df, pot_father_df], ignore_index=True)
                    # 2 old CODIS duo
                    for ii in ids2:
                        p_pair_df = temp_parents_df.loc[[ii, mother]]
                        pot_father_df = pd.DataFrame(columns=columns_list, index=[0])
                        pot_father_df.iloc[0]["№"] = temp_parents_df.loc[ii]["№"] + "_old_codis_duo"
                        q_old_codis_duo_ref = []
                        q_old_codis_duo_rus = []
                        for elem in range(1, len(columns_list) - 34, 2):
                            allele = columns_list[elem].split('_')[0]
                            pot_father_df.iloc[0][columns_list[elem]] = p_pair_df.iloc[0][columns_list[elem]]
                            pot_father_df.iloc[0][columns_list[elem + 1]] = p_pair_df.iloc[0][columns_list[elem + 1]]
                            if (allele not in ids2[ii]) and (allele in config["codis_old"]):
                                kid_alleles = [kids_df.iloc[k][columns_list[elem]],
                                               kids_df.iloc[k][columns_list[elem + 1]]]
                                kid_alleles.sort()
                                f_allele_duo, other_allele_duo, knowledge_duo = duo_freq_father_allele(kid_alleles[0],
                                                                                                       kid_alleles[1])
                                q_old_duo_ref = calculate_q(config["codis_old"], allele, ref_dict, f_allele_duo,
                                                            other_allele_duo, knowledge_duo)
                                q_old_duo_rus = calculate_q(config["codis_old"], allele, config["rus_dict"],
                                                            f_allele_duo, other_allele_duo, knowledge_duo)
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
                        pot_father_df.iloc[0]["PI_old_CODIS_duo_ref"] = 1 / multiplication_old_duo_ref
                        pot_father_df.iloc[0]["PP_old_CODIS_duo_ref"] = 1 / (1 + multiplication_old_duo_ref)
                        pot_father_df.iloc[0]["PI_old_CODIS_duo_rus"] = 1 / multiplication_old_duo_rus
                        pot_father_df.iloc[0]["PP_old_CODIS_duo_rus"] = 1 / (1 + multiplication_old_duo_rus)
                        pot_father_df.iloc[0]["population"] = p_pair_df.iloc[0]["population"]
                        pot_father_df.iloc[0]["Child_ID"] = kids_df.iloc[k]["№"]
                        pot_father_df.iloc[0]["Status"] = "false_positive_father"
                        pot_father_df.iloc[0]["number_of_mutations"] = len(ids2[ii])
                        offsprings_df = pd.concat([offsprings_df, pot_father_df], ignore_index=True)
                    # 3 new CODIS trio
                    for ii in ids3:
                        p_pair_df = temp_parents_df.loc[[ii, mother]]
                        pot_father_df = pd.DataFrame(columns=columns_list, index=[0])
                        pot_father_df.iloc[0]["№"] = temp_parents_df.loc[ii]["№"] + "_new_codis_trio"
                        q_new_codis_trio_ref = []
                        q_new_codis_trio_rus = []
                        for elem in range(1, len(columns_list) - 34, 2):
                            allele = columns_list[elem].split('_')[0]
                            pot_father_df.iloc[0][columns_list[elem]] = p_pair_df.iloc[0][columns_list[elem]]
                            pot_father_df.iloc[0][columns_list[elem + 1]] = p_pair_df.iloc[0][columns_list[elem + 1]]
                            if (allele not in ids3[ii]) and (allele in config["codis_new"]):
                                m_alleles = [p_pair_df.iloc[1][columns_list[elem]],
                                             p_pair_df.iloc[1][columns_list[elem + 1]]]
                                kid_alleles = [kids_df.iloc[k][columns_list[elem]],
                                               kids_df.iloc[k][columns_list[elem + 1]]]
                                m_alleles.sort()
                                kid_alleles.sort()
                                f_allele_trio, other_allele_trio, knowledge_trio = \
                                    trio_freq_father_allele(m_alleles[0], m_alleles[1], kid_alleles[0], kid_alleles[1])
                                q_new_trio_ref = calculate_q(config["codis_new"], allele, ref_dict, f_allele_trio,
                                                             other_allele_trio, knowledge_trio)
                                q_new_trio_rus = calculate_q(config["codis_new"], allele, config["rus_dict"],
                                                             f_allele_trio, other_allele_trio, knowledge_trio)
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
                        pot_father_df.iloc[0]["PI_new_CODIS_trio_ref"] = 1 / multiplication_new_trio_ref
                        pot_father_df.iloc[0]["PP_new_CODIS_trio_ref"] = 1 / (1 + multiplication_new_trio_ref)
                        pot_father_df.iloc[0]["PI_new_CODIS_trio_rus"] = 1 / multiplication_new_trio_rus
                        pot_father_df.iloc[0]["PP_new_CODIS_trio_rus"] = 1 / (1 + multiplication_new_trio_rus)
                        pot_father_df.iloc[0]["population"] = p_pair_df.iloc[0]["population"]
                        pot_father_df.iloc[0]["Child_ID"] = kids_df.iloc[k]["№"]
                        pot_father_df.iloc[0]["Status"] = "false_positive_father"
                        pot_father_df.iloc[0]["number_of_mutations"] = len(ids3[ii])
                        offsprings_df = pd.concat([offsprings_df, pot_father_df], ignore_index=True)
                    # 4 new CODIS duo
                    for ii in ids4:
                        p_pair_df = temp_parents_df.loc[[ii, mother]]
                        pot_father_df = pd.DataFrame(columns=columns_list, index=[0])
                        pot_father_df.iloc[0]["№"] = temp_parents_df.loc[ii]["№"] + "_new_codis_duo"
                        q_new_codis_duo_ref = []
                        q_new_codis_duo_rus = []
                        for elem in range(1, len(columns_list) - 34, 2):
                            allele = columns_list[elem].split('_')[0]
                            pot_father_df.iloc[0][columns_list[elem]] = p_pair_df.iloc[0][columns_list[elem]]
                            pot_father_df.iloc[0][columns_list[elem + 1]] = p_pair_df.iloc[0][columns_list[elem + 1]]
                            if (allele not in ids4[ii]) and (allele in config["codis_new"]):
                                kid_alleles = [kids_df.iloc[k][columns_list[elem]],
                                               kids_df.iloc[k][columns_list[elem + 1]]]
                                kid_alleles.sort()
                                f_allele_duo, other_allele_duo, knowledge_duo = duo_freq_father_allele(kid_alleles[0],
                                                                                                       kid_alleles[1])
                                q_new_duo_ref = calculate_q(config["codis_new"], allele, ref_dict, f_allele_duo,
                                                            other_allele_duo, knowledge_duo)
                                q_new_duo_rus = calculate_q(config["codis_new"], allele, config["rus_dict"],
                                                            f_allele_duo, other_allele_duo, knowledge_duo)
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
                        pot_father_df.iloc[0]["PI_new_CODIS_duo_ref"] = 1 / multiplication_new_duo_ref
                        pot_father_df.iloc[0]["PP_new_CODIS_duo_ref"] = 1 / (1 + multiplication_new_duo_ref)
                        pot_father_df.iloc[0]["PI_new_CODIS_duo_rus"] = 1 / multiplication_new_duo_rus
                        pot_father_df.iloc[0]["PP_new_CODIS_duo_rus"] = 1 / (1 + multiplication_new_duo_rus)
                        pot_father_df.iloc[0]["population"] = p_pair_df.iloc[0]["population"]
                        pot_father_df.iloc[0]["Child_ID"] = kids_df.iloc[k]["№"]
                        pot_father_df.iloc[0]["Status"] = "false_positive_father"
                        pot_father_df.iloc[0]["number_of_mutations"] = len(ids4[ii])
                        offsprings_df = pd.concat([offsprings_df, pot_father_df], ignore_index=True)
                    child_counter += 1

                    #################################################
                    # CALCULATE LR FOR REAL PARENTS
                    df_for_lr = offsprings_df.loc[(offsprings_df['Child_ID'] == kids_df.iloc[k]["№"])]
                    inx_lr = df_for_lr.index.values.tolist()
                    real_father_df = df_for_lr.iloc[[0]]  # real father data
                    real_mother_df = df_for_lr.iloc[[1]]  # real mother data
                    child_df = df_for_lr.iloc[[2]]  # df with 1st, 2nd or 3rd child's data
                    hyp1_ref_trio_old = []
                    hyp2_ref_trio_old = []
                    hyp1_rus_trio_old = []
                    hyp2_rus_trio_old = []
                    hyp1_ref_trio_new = []
                    hyp2_ref_trio_new = []
                    hyp1_rus_trio_new = []
                    hyp2_rus_trio_new = []
                    hyp1_ref_duo_old = []
                    hyp2_ref_duo_old = []
                    hyp1_rus_duo_old = []
                    hyp2_rus_duo_old = []
                    hyp1_ref_duo_new = []
                    hyp2_ref_duo_new = []
                    hyp1_rus_duo_new = []
                    hyp2_rus_duo_new = []
                    for loc in range(1, len(columns_list) - 34, 2):
                        el = columns_list[loc].split('_')[0]
                        rf = [real_father_df.iloc[0][columns_list[loc]], real_father_df.iloc[0][columns_list[loc + 1]]]
                        rm = [real_mother_df.iloc[0][columns_list[loc]], real_mother_df.iloc[0][columns_list[loc + 1]]]
                        rc = [child_df.iloc[0][columns_list[loc]], child_df.iloc[0][columns_list[loc + 1]]]
                        rf.sort()
                        rm.sort()
                        rc.sort()
                        # check if there is a mutation:
                        lack_of_mut_trio = father_trio(rf[0], rf[1], rm[0], rm[1], rc[0], rc[1])
                        if lack_of_mut_trio:
                            p1, p2, p3, p4 = hyp_trio_wo_mut(rm[0], rm[1], rc[0], rc[1], ref_dict,
                                                             config["rus_dict"], el)
                        else:
                            p1, p2, p3, p4 = hyp_trio_mut(rm[0], rm[1], rc[0], rc[1], ref_dict, config["rus_dict"],
                                                          el, config["mut_rate"])
                        lack_of_mut_duo = father_duo(rf[0], rf[1], rc[0], rc[1])  # check if there is mutation
                        if lack_of_mut_duo:
                            p5, p6, p7, p8 = hyp_duo_wo_mut(rf[0], rf[1], rc[0], rc[1], ref_dict,
                                                            config["rus_dict"], el)
                        else:
                            p5, p6, p7, p8 = hyp_duo_mut(rc[0], rc[1], ref_dict, config["rus_dict"], el,
                                                         config["mut_rate"])
                        if el in config["codis_old"]:
                            hyp1_ref_trio_old.append(p1)
                            hyp2_ref_trio_old.append(p2)
                            hyp1_rus_trio_old.append(p3)
                            hyp2_rus_trio_old.append(p4)
                            hyp1_ref_duo_old.append(p5)
                            hyp2_ref_duo_old.append(p6)
                            hyp1_rus_duo_old.append(p7)
                            hyp2_rus_duo_old.append(p8)
                        if el in config["codis_new"]:
                            hyp1_ref_trio_new.append(p1)
                            hyp2_ref_trio_new.append(p2)
                            hyp1_rus_trio_new.append(p3)
                            hyp2_rus_trio_new.append(p4)
                            hyp1_ref_duo_new.append(p5)
                            hyp2_ref_duo_new.append(p6)
                            hyp1_rus_duo_new.append(p7)
                            hyp2_rus_duo_new.append(p8)
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
                    for m in hyp1_ref_trio_old:
                        multiplication_1 *= m
                    for m in hyp2_ref_trio_old:
                        multiplication_2 *= m
                    for m in hyp1_rus_trio_old:
                        multiplication_3 *= m
                    for m in hyp2_rus_trio_old:
                        multiplication_4 *= m
                    for m in hyp1_ref_trio_new:
                        multiplication_5 *= m
                    for m in hyp2_ref_trio_new:
                        multiplication_6 *= m
                    for m in hyp1_rus_trio_new:
                        multiplication_7 *= m
                    for m in hyp2_rus_trio_new:
                        multiplication_8 *= m
                    for m in hyp1_ref_duo_old:
                        multiplication_9 *= m
                    for m in hyp2_ref_duo_old:
                        multiplication_10 *= m
                    for m in hyp1_rus_duo_old:
                        multiplication_11 *= m
                    for m in hyp2_rus_duo_old:
                        multiplication_12 *= m
                    for m in hyp1_ref_duo_new:
                        multiplication_13 *= m
                    for m in hyp2_ref_duo_new:
                        multiplication_14 *= m
                    for m in hyp1_rus_duo_new:
                        multiplication_15 *= m
                    for m in hyp2_rus_duo_new:
                        multiplication_16 *= m
                    offsprings_df.iloc[inx_lr[0]]["LR_ref_trio_old_codis"] = multiplication_1 / multiplication_2
                    offsprings_df.iloc[inx_lr[0]]["LR_rus_trio_old_codis"] = multiplication_3 / multiplication_4
                    offsprings_df.iloc[inx_lr[0]]["LR_ref_trio_new_codis"] = multiplication_5 / multiplication_6
                    offsprings_df.iloc[inx_lr[0]]["LR_rus_trio_new_codis"] = multiplication_7 / multiplication_8
                    offsprings_df.iloc[inx_lr[0]]["LR_ref_duo_old_codis"] = multiplication_9 / multiplication_10
                    offsprings_df.iloc[inx_lr[0]]["LR_rus_duo_old_codis"] = multiplication_11 / multiplication_12
                    offsprings_df.iloc[inx_lr[0]]["LR_ref_duo_new_codis"] = multiplication_13 / multiplication_14
                    offsprings_df.iloc[inx_lr[0]]["LR_rus_duo_new_codis"] = multiplication_15 / multiplication_16

                    #################################################
                    # CALCULATE LR FOR FALSE POSITIVE FATHERS
                    for ii in range(3, len(df_for_lr.index)):
                        p_hyp1_list_ref = []
                        p_hyp2_list_ref = []
                        p_hyp1_list_rus = []
                        p_hyp2_list_rus = []
                        steps = []
                        one_pf_df = offsprings_df.iloc[[df_for_lr.index[ii]]]  # df with one PF data
                        name = one_pf_df.iloc[0]['№'].split("_")
                        if ("old" in name) and ("trio" in name):
                            for loc in range(1, len(columns_list) - 34, 2):
                                el = columns_list[loc].split('_')[0]
                                if el in config["codis_old"]:
                                    pf = [one_pf_df.iloc[0][columns_list[loc]],
                                          one_pf_df.iloc[0][columns_list[loc + 1]]]
                                    rm = [real_mother_df.iloc[0][columns_list[loc]],
                                          real_mother_df.iloc[0][columns_list[loc + 1]]]
                                    rc = [child_df.iloc[0][columns_list[loc]], child_df.iloc[0][columns_list[loc + 1]]]
                                    pf.sort()
                                    rm.sort()
                                    rc.sort()
                                    # check if there is a mutation:
                                    lack_of_mut = father_trio(pf[0], pf[1], rm[0], rm[1], rc[0], rc[1])
                                    if lack_of_mut:
                                        p_1, p_2, p_3, p_4 = hyp_trio_wo_mut(rm[0], rm[1], rc[0], rc[1], ref_dict,
                                                                             config["rus_dict"], el)
                                    else:
                                        steps.append(count_step_trio(pf[0], pf[1], rm[0], rm[1], rc[0], rc[1]))
                                        p_1, p_2, p_3, p_4 = hyp_trio_mut(rm[0], rm[1], rc[0], rc[1], ref_dict,
                                                                          config["rus_dict"], el, config["mut_rate"])
                                    p_hyp1_list_ref.append(p_1)
                                    p_hyp2_list_ref.append(p_2)
                                    p_hyp1_list_rus.append(p_3)
                                    p_hyp2_list_rus.append(p_4)
                            multiplication_1 = 1
                            multiplication_2 = 1
                            multiplication_3 = 1
                            multiplication_4 = 1
                            for m in p_hyp1_list_ref:
                                multiplication_1 *= m
                            for m in p_hyp2_list_ref:
                                multiplication_2 *= m
                            for m in p_hyp1_list_rus:
                                multiplication_3 *= m
                            for m in p_hyp2_list_rus:
                                multiplication_4 *= m
                            offsprings_df.iloc[df_for_lr.index[ii]]["LR_ref_trio_old_codis"] = (multiplication_1 /
                                                                                                multiplication_2)
                            offsprings_df.iloc[df_for_lr.index[ii]]["LR_rus_trio_old_codis"] = (multiplication_3 /
                                                                                                multiplication_4)
                            offsprings_df.iloc[df_for_lr.index[ii]]["step_mutation"] = steps
                        if "old" in name and "duo" in name:
                            for loc in range(1, len(columns_list) - 34, 2):
                                el = columns_list[loc].split('_')[0]
                                if el in config["codis_old"]:
                                    pf = [one_pf_df.iloc[0][columns_list[loc]],
                                          one_pf_df.iloc[0][columns_list[loc + 1]]]
                                    rc = [child_df.iloc[0][columns_list[loc]], child_df.iloc[0][columns_list[loc + 1]]]
                                    pf.sort()
                                    rc.sort()
                                    lack_of_mut = father_duo(pf[0], pf[1], rc[0], rc[1])  # check if there is mutation
                                    if lack_of_mut:
                                        p_1, p_2, p_3, p_4 = hyp_duo_wo_mut(pf[0], pf[1], rc[0], rc[1], ref_dict,
                                                                            config["rus_dict"], el)
                                    else:
                                        steps.append(count_step_duo(pf[0], pf[1], rc[0], rc[1]))
                                        p_1, p_2, p_3, p_4 = hyp_duo_mut(rc[0], rc[1], ref_dict, config["rus_dict"], el,
                                                                         config["mut_rate"])
                                    p_hyp1_list_ref.append(p_1)
                                    p_hyp2_list_ref.append(p_2)
                                    p_hyp1_list_rus.append(p_3)
                                    p_hyp2_list_rus.append(p_4)
                            multiplication_1 = 1
                            multiplication_2 = 1
                            multiplication_3 = 1
                            multiplication_4 = 1
                            for m in p_hyp1_list_ref:
                                multiplication_1 *= m
                            for m in p_hyp2_list_ref:
                                multiplication_2 *= m
                            for m in p_hyp1_list_rus:
                                multiplication_3 *= m
                            for m in p_hyp2_list_rus:
                                multiplication_4 *= m
                            offsprings_df.iloc[df_for_lr.index[ii]]["LR_ref_duo_old_codis"] = (multiplication_1 /
                                                                                               multiplication_2)
                            offsprings_df.iloc[df_for_lr.index[ii]]["LR_rus_duo_old_codis"] = (multiplication_3 /
                                                                                               multiplication_4)
                            offsprings_df.iloc[df_for_lr.index[ii]]["step_mutation"] = steps
                        if "new" in name and "trio" in name:
                            for loc in range(1, len(columns_list) - 34, 2):
                                el = columns_list[loc].split('_')[0]
                                if el in config["codis_new"]:
                                    pf = [one_pf_df.iloc[0][columns_list[loc]],
                                          one_pf_df.iloc[0][columns_list[loc + 1]]]
                                    rm = [real_mother_df.iloc[0][columns_list[loc]],
                                          real_mother_df.iloc[0][columns_list[loc + 1]]]
                                    rc = [child_df.iloc[0][columns_list[loc]], child_df.iloc[0][columns_list[loc + 1]]]
                                    pf.sort()
                                    rm.sort()
                                    rc.sort()
                                    # check if there is mutation:
                                    lack_of_mut = father_trio(pf[0], pf[1], rm[0], rm[1], rc[0], rc[1])
                                    if lack_of_mut:
                                        p_1, p_2, p_3, p_4 = hyp_trio_wo_mut(rm[0], rm[1], rc[0], rc[1], ref_dict,
                                                                             config["rus_dict"], el)
                                    else:
                                        steps.append(count_step_trio(pf[0], pf[1], rm[0], rm[1], rc[0], rc[1]))
                                        p_1, p_2, p_3, p_4 = hyp_trio_mut(rm[0], rm[1], rc[0], rc[1], ref_dict,
                                                                          config["rus_dict"], el, config["mut_rate"])
                                    p_hyp1_list_ref.append(p_1)
                                    p_hyp2_list_ref.append(p_2)
                                    p_hyp1_list_rus.append(p_3)
                                    p_hyp2_list_rus.append(p_4)
                            multiplication_1 = 1
                            multiplication_2 = 1
                            multiplication_3 = 1
                            multiplication_4 = 1
                            for m in p_hyp1_list_ref:
                                multiplication_1 *= m
                            for m in p_hyp2_list_ref:
                                multiplication_2 *= m
                            for m in p_hyp1_list_rus:
                                multiplication_3 *= m
                            for m in p_hyp2_list_rus:
                                multiplication_4 *= m
                            offsprings_df.iloc[df_for_lr.index[ii]]["LR_ref_trio_new_codis"] = (multiplication_1 /
                                                                                                multiplication_2)
                            offsprings_df.iloc[df_for_lr.index[ii]]["LR_rus_trio_new_codis"] = (multiplication_3 /
                                                                                                multiplication_4)
                            offsprings_df.iloc[df_for_lr.index[ii]]["step_mutation"] = steps
                        if "new" in name and "duo" in name:
                            for loc in range(1, len(columns_list) - 34, 2):
                                el = columns_list[loc].split('_')[0]
                                if el in config["codis_new"]:
                                    pf = [one_pf_df.iloc[0][columns_list[loc]],
                                          one_pf_df.iloc[0][columns_list[loc + 1]]]
                                    rc = [child_df.iloc[0][columns_list[loc]], child_df.iloc[0][columns_list[loc + 1]]]
                                    pf.sort()
                                    rc.sort()
                                    lack_of_mut = father_duo(pf[0], pf[1], rc[0], rc[1])  # check if there is mutation
                                    if lack_of_mut:
                                        p_1, p_2, p_3, p_4 = hyp_duo_wo_mut(pf[0], pf[1], rc[0], rc[1], ref_dict,
                                                                            config["rus_dict"], el)
                                    else:
                                        steps.append(count_step_duo(pf[0], pf[1], rc[0], rc[1]))
                                        p_1, p_2, p_3, p_4 = hyp_duo_mut(rc[0], rc[1], ref_dict, config["rus_dict"], el,
                                                                         config["mut_rate"])
                                    p_hyp1_list_ref.append(p_1)
                                    p_hyp2_list_ref.append(p_2)
                                    p_hyp1_list_rus.append(p_3)
                                    p_hyp2_list_rus.append(p_4)
                            multiplication_1 = 1
                            multiplication_2 = 1
                            multiplication_3 = 1
                            multiplication_4 = 1
                            for m in p_hyp1_list_ref:
                                multiplication_1 *= m
                            for m in p_hyp2_list_ref:
                                multiplication_2 *= m
                            for m in p_hyp1_list_rus:
                                multiplication_3 *= m
                            for m in p_hyp2_list_rus:
                                multiplication_4 *= m
                            offsprings_df.iloc[df_for_lr.index[ii]]["LR_ref_duo_new_codis"] = (multiplication_1 /
                                                                                               multiplication_2)
                            offsprings_df.iloc[df_for_lr.index[ii]]["LR_rus_duo_new_codis"] = (multiplication_3 /
                                                                                               multiplication_4)
                            offsprings_df.iloc[df_for_lr.index[ii]]["step_mutation"] = steps
                n -= 1
    # offsprings_df.drop(columns=["groups"], axis=1, inplace=True)    # parents data in output
    offsprings_df.to_excel("output_gen_gen_TESTTEST.xlsx", index=True)

if __name__ == "__main__":
    main()
