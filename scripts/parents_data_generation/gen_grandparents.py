import copy
import pandas as pd
import random
import time
import yaml


# Generate and count the LR for

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
def hyp_grandparents(g1, g2, k1, k2, our_dict, rus_dict, str_elem):
    #  1 - for our frequencies, 2 - for rus
    p1, p2 = calculate_p(g1, g2, k1, k2, our_dict, str_elem)
    p3, p4 = calculate_p(g1, g2, k1, k2, rus_dict, str_elem)
    return p1, p2, p3, p4


def calculate_p(g1, g2, k1, k2, dic, allele):
    print(g1, g2, k1, k2, allele)
    p1 = p2 = None
    if k1 == k2:
        if g1 == g2 and g1 == k1:
            if k1 in dic[allele]:
                p1 = dic[allele][k1] * (2 - dic[allele][k1])
                p2 = (dic[allele][k1] * (2 - dic[allele][k1])) ** 2
            else:
                p1 = dic[allele]["pmin"] * (2 - dic[allele]["pmin"])
                p2 = (dic[allele]["pmin"] * (2 - dic[allele]["pmin"])) ** 2
        elif k1 == g1 or k1 == g2:
            if k1 == g1:
                a = g1
            elif k1 == g2:
                a = g2
            if a in dic[allele]:
                p1 = (0.5 + 0.5 * dic[allele][a]) * (dic[allele][a] * (2 - dic[allele][a]))
                p2 = (dic[allele][a] * (2 - dic[allele][a])) ** 2
            elif a not in dic[allele]:
                p1 = (0.5 + 0.5 * dic[allele]["pmin"]) * (dic[allele]["pmin"] * (2 - dic[allele]["pmin"]))
                p2 = (dic[allele]["pmin"] * (2 - dic[allele]["pmin"])) ** 2
        elif k1 != g2 and k1 != g2:
            a = k1
            if a in dic[allele]:
                p1 = dic[allele][a] * (dic[allele][a] * (2 - dic[allele][a]))
                p2 = (dic[allele][a] * (2 - dic[allele][a])) ** 2
            elif a not in dic[allele]:
                p1 = dic[allele]["pmin"] * (dic[allele]["pmin"] * (2 - dic[allele]["pmin"]))
                p2 = (dic[allele]["pmin"] * (2 - dic[allele]["pmin"])) ** 2
    else:
        if g1 == g2 and (g1 == k1 or g1 == k2):
            a = g1
            b = None
            if g1 == k1:
                b = k2
            elif g1 == k2:
                b = k1
            if a in dic[allele] and b in dic[allele]:
                p1 = dic[allele][b] * (2 - dic[allele][b]) + dic[allele][b] \
                     * (dic[allele][a] * (2 - dic[allele][a]) - 2 * dic[allele][a] * dic[allele][b])
                p2 = 2 * dic[allele][a] * (2 - dic[allele][a]) * dic[allele][b] * (2 - dic[allele][b]) - \
                     2 * dic[allele][a] * dic[allele][b] * 2 * dic[allele][a] * dic[allele][b]
            elif a not in dic[allele] and b in dic[allele]:
                p1 = dic[allele][b] * (2 - dic[allele][b]) + dic[allele][b] \
                     * (dic[allele]["pmin"] * (2 - dic[allele]["pmin"]) - 2 * dic[allele]["pmin"] * dic[allele][b])
                p2 = 2 * dic[allele]["pmin"] * (2 - dic[allele]["pmin"]) * dic[allele][b] * (2 - dic[allele][b]) - \
                     2 * dic[allele]["pmin"] * dic[allele][b] * 2 * dic[allele]["pmin"] * dic[allele][b]
            elif a in dic[allele] and b not in dic[allele]:
                p1 = dic[allele]["pmin"] * (2 - dic[allele]["pmin"]) + dic[allele]["pmin"] \
                     * (dic[allele][a] * (2 - dic[allele][a]) - 2 * dic[allele][a] * dic[allele]["pmin"])
                p2 = 2 * dic[allele][a] * (2 - dic[allele][a]) * dic[allele]["pmin"] * (2 - dic[allele]["pmin"]) - \
                     2 * dic[allele][a] * dic[allele]["pmin"] * 2 * dic[allele][a] * dic[allele]["pmin"]
            elif a not in dic[allele] and b not in dic[allele]:
                p1 = dic[allele]["pmin"] * (2 - dic[allele]["pmin"]) + dic[allele]["pmin"] \
                     * (dic[allele]["pmin"] * (2 - dic[allele]["pmin"]) - 2 * dic[allele]["pmin"] * dic[allele]["pmin"])
                p2 = 2 * dic[allele]["pmin"] * (2 - dic[allele]["pmin"]) * dic[allele]["pmin"] * (2 - dic[allele]["pmin"]) - \
                     2 * dic[allele]["pmin"] * dic[allele]["pmin"] * 2 * dic[allele]["pmin"] * dic[allele]["pmin"]
        if g1 != g2 and k1 == g1 and k2 == g2:
            a = k1
            b = k2
            if a in dic[allele] and b in dic[allele]:
                p1 = (0.5 + 0.5 * dic[allele][a]) * dic[allele][b] * (2 - dic[allele][b]) + \
                     (0.5 + 0.5 * dic[allele][b]) * dic[allele][a] * (2 - dic[allele][a]) - \
                     0.5 * (dic[allele][a] + dic[allele][b]) * 2 * dic[allele][a] * dic[allele][b]
                p2 = 2 * dic[allele][a] * (2 - dic[allele][a]) * dic[allele][b] * (2 - dic[allele][b]) - \
                     2 * dic[allele][a] * dic[allele][b] * 2 * dic[allele][a] * dic[allele][b]
            elif a not in dic[allele] and b in dic[allele]:
                p1 = (0.5 + 0.5 * dic[allele]["pmin"]) * dic[allele][b] * (2 - dic[allele][b]) + \
                     (0.5 + 0.5 * dic[allele][b]) * dic[allele]["pmin"] * (2 - dic[allele]["pmin"]) - \
                     0.5 * (dic[allele]["pmin"] + dic[allele][b]) * 2 * dic[allele]["pmin"] * dic[allele][b]
                p2 = 2 * dic[allele]["pmin"] * (2 - dic[allele]["pmin"]) * dic[allele][b] * (2 - dic[allele][b]) - \
                     2 * dic[allele]["pmin"] * dic[allele][b] * 2 * dic[allele]["pmin"] * dic[allele][b]
            elif a in dic[allele] and b not in dic[allele]:
                p1 = (0.5 + 0.5 * dic[allele][a]) * dic[allele]["pmin"] * (2 - dic[allele]["pmin"]) + \
                     (0.5 + 0.5 * dic[allele]["pmin"]) * dic[allele][a] * (2 - dic[allele][a]) - \
                     0.5 * (dic[allele][a] + dic[allele]["pmin"]) * 2 * dic[allele][a] * dic[allele]["pmin"]
                p2 = 2 * dic[allele][a] * (2 - dic[allele][a]) * dic[allele]["pmin"] * (2 - dic[allele]["pmin"]) - \
                     2 * dic[allele][a] * dic[allele]["pmin"] * 2 * dic[allele][a] * dic[allele]["pmin"]
            elif a not in dic[allele] and b not in dic[allele]:
                p1 = (0.5 + 0.5 * dic[allele]["pmin"]) * dic[allele]["pmin"] * (2 - dic[allele]["pmin"]) + \
                     (0.5 + 0.5 * dic[allele]["pmin"]) * dic[allele]["pmin"] * (2 - dic[allele]["pmin"]) - \
                     0.5 * (dic[allele]["pmin"] + dic[allele]["pmin"]) * 2 * dic[allele]["pmin"] * dic[allele]["pmin"]
                p2 = 2 * dic[allele]["pmin"] * (2 - dic[allele]["pmin"]) * dic[allele]["pmin"] * \
                     (2 - dic[allele]["pmin"]) - 2 * dic[allele]["pmin"] * dic[allele]["pmin"] * \
                     2 * dic[allele]["pmin"] * dic[allele]["pmin"]
        if g1 != g2 and (k1 == g1 and k2 != g2 or k1 == g2 and k2 != g1 or k1 != g1 and k2 == g2 or k1 != g2 and k2 == g1):
            a = None
            b = None
            if k1 == g1 and k2 != g2:
                a = k1
                b = k2
            if k1 == g2 and k2 != g1:
                a = k1
                b = k2
            if k1 != g1 and k2 == g2:
                a = k2
                b = k1
            if k1 != g2 and k2 == g1:
                a = k2
                b = k1
            if a in dic[allele] and b in dic[allele]:
                p1 = (0.5 + 0.5 * dic[allele][a]) * dic[allele][b] * (2 - dic[allele][b]) + \
                     dic[allele][b] * dic[allele][a] * (2 - dic[allele][a]) - 0.5 * dic[allele][b] * \
                     2 * dic[allele][a] * dic[allele][b]
                p2 = 2 * dic[allele][a] * (2 - dic[allele][a]) * dic[allele][b] * (2 - dic[allele][b]) - \
                     2 * dic[allele][a] * dic[allele][b] * 2 * dic[allele][a] * dic[allele][b]
            elif a not in dic[allele] and b in dic[allele]:
                p1 = (0.5 + 0.5 * dic[allele]["pmin"]) * dic[allele][b] * (2 - dic[allele][b]) + \
                     dic[allele][b] * dic[allele]["pmin"] * (2 - dic[allele]["pmin"]) - 0.5 * dic[allele][b] * \
                     2 * dic[allele]["pmin"] * dic[allele][b]
                p2 = 2 * dic[allele]["pmin"] * (2 - dic[allele]["pmin"]) * dic[allele][b] * (2 - dic[allele][b]) - \
                     2 * dic[allele]["pmin"] * dic[allele][b] * 2 * dic[allele]["pmin"] * dic[allele][b]
            elif a in dic[allele] and b not in dic[allele]:
                p1 = (0.5 + 0.5 * dic[allele][a]) * dic[allele]["pmin"] * (2 - dic[allele]["pmin"]) + \
                     dic[allele]["pmin"] * dic[allele][a] * (2 - dic[allele][a]) - 0.5 * dic[allele]["pmin"] * \
                     2 * dic[allele][a] * dic[allele]["pmin"]
                p2 = 2 * dic[allele][a] * (2 - dic[allele][a]) * dic[allele]["pmin"] * (2 - dic[allele]["pmin"]) - \
                     2 * dic[allele][a] * dic[allele]["pmin"] * 2 * dic[allele][a] * dic[allele]["pmin"]
            elif a not in dic[allele] and b not in dic[allele]:
                p1 = (0.5 + 0.5 * dic[allele]["pmin"]) * dic[allele]["pmin"] * (2 - dic[allele]["pmin"]) + \
                     dic[allele]["pmin"] * dic[allele]["pmin"] * (2 - dic[allele]["pmin"]) - 0.5 * dic[allele]["pmin"] * \
                     2 * dic[allele]["pmin"] * dic[allele]["pmin"]
                p2 = 2 * dic[allele]["pmin"] * (2 - dic[allele]["pmin"]) * dic[allele]["pmin"] * (2 - dic[allele]["pmin"]) - \
                     2 * dic[allele]["pmin"] * dic[allele]["pmin"] * 2 * dic[allele]["pmin"] * dic[allele]["pmin"]
        if k1 != g1 and k2 != g2 and k1 != g2 and k2 != g1:
            a = k1
            b = k2
            if a in dic[allele] and b in dic[allele]:
                p1 = dic[allele][a] * dic[allele][b] * (2 - dic[allele][b]) + dic[allele][b] * \
                     dic[allele][a] * (2 - dic[allele][a])
                p2 = 2 * dic[allele][a] * (2 - dic[allele][a]) * dic[allele][b] * (2 - dic[allele][b]) -\
                     2 * dic[allele][a] * dic[allele][b] * 2 * dic[allele][a] * dic[allele][b]
            elif a not in dic[allele] and b in dic[allele]:
                p1 = dic[allele]["pmin"] * dic[allele][b] * (2 - dic[allele][b]) + dic[allele][b] * \
                     dic[allele]["pmin"] * (2 - dic[allele]["pmin"])
                p2 = 2 * dic[allele]["pmin"] * (2 - dic[allele]["pmin"]) * dic[allele][b] * (2 - dic[allele][b]) -\
                     2 * dic[allele]["pmin"] * dic[allele][b] * 2 * dic[allele]["pmin"] * dic[allele][b]
            elif a in dic[allele] and b not in dic[allele]:
                p1 = dic[allele][a] * dic[allele]["pmin"] * (2 - dic[allele]["pmin"]) + dic[allele]["pmin"] * \
                     dic[allele][a] * (2 - dic[allele][a])
                p2 = 2 * dic[allele][a] * (2 - dic[allele][a]) * dic[allele]["pmin"] * (2 - dic[allele]["pmin"]) -\
                     2 * dic[allele][a] * dic[allele]["pmin"] * 2 * dic[allele][a] * dic[allele]["pmin"]
            elif a not in dic[allele] and b not in dic[allele]:
                p1 = dic[allele]["pmin"] * dic[allele]["pmin"] * (2 - dic[allele]["pmin"]) + dic[allele]["pmin"] * \
                     dic[allele]["pmin"] * (2 - dic[allele]["pmin"])
                p2 = 2 * dic[allele]["pmin"] * (2 - dic[allele]["pmin"]) * dic[allele]["pmin"] * (2 - dic[allele]["pmin"]) -\
                     2 * dic[allele]["pmin"] * dic[allele]["pmin"] * 2 * dic[allele]["pmin"] * dic[allele]["pmin"]
    return p1, p2


def main():
    start = time.time()
    with open("../config_file_grandparents.yaml", "r") as yaml_file:
        config = yaml.load(yaml_file, Loader=yaml.FullLoader)
    # Create dataframe with offsprings
    grandparents_df = pd.read_excel(config["main_table_path"])
    columns_list = grandparents_df.columns.values.tolist()
    columns_list = columns_list[:-1]
    columns_list.append("Family_ID")
    columns_list.append("Status")
    offsprings_df = pd.DataFrame(columns=columns_list)    # Create main df
    for pop in config["populations"]:
        ref_dict = empty_freq_table(config["main_table_path"])
        ref_dict = calculate_frequencies(pop, ref_dict, config["main_table_path"])
        new_pair_df = pd.DataFrame(columns=columns_list, index=[0, 1])
        temp_parents_df = grandparents_df.loc[(grandparents_df['population'] == pop)]
        index_list = temp_parents_df.index.values.tolist()
        index_list_f = index_list[:4000]
        index_list_m = index_list[4000:]
        list_for_pairs = []
        n = copy.deepcopy(config["number_of_pairs_grands"])
        while n != 0:
            father = random.choice(index_list_f)
            mother = random.choice(index_list_m)
            if father not in list_for_pairs and mother not in list_for_pairs and father != mother:
                list_for_pairs.append(father)
                list_for_pairs.append(mother)
                pair_df = temp_parents_df.loc[[father, mother]]    # create df for one random pair
                kids_df = pd.DataFrame(columns=columns_list, index=[i for i in range(0, config["kids"])])
                for k in range(config["kids"]):
                    kids_df.iloc[k]["Status"] = random.choice(["father", "mother"])
                    kids_df.iloc[k]["№"] = (str(pair_df.iloc[0]['№']) + "-" + str(pair_df.iloc[1]['№']) + "-" +
                                            str(k + 1))
                    for elem in range(1, len(columns_list) - 3, 2):
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
                        for elem in range(len(columns_list) - 2):
                            new_pair_df.iloc[0][columns_list[elem]] = pair_df.iloc[0][columns_list[elem]]
                            new_pair_df.iloc[1][columns_list[elem]] = pair_df.iloc[1][columns_list[elem]]
                        new_pair_df.iloc[0]["Status"] = "grandfather"
                        new_pair_df.iloc[1]["Status"] = "grandmother"
                        offsprings_df = pd.concat([offsprings_df, new_pair_df], ignore_index=True)
                    kids_df.iloc[k]["population"] = pop
                    offsprings_df = pd.concat([offsprings_df, kids_df.loc[[k]]], ignore_index=True)
                n -= 1
        # Generate grandkids
        fathers_df = offsprings_df.loc[(offsprings_df['Status'] == "father")]
        mothers_df = offsprings_df.loc[(offsprings_df['Status'] == "mother")]
        columns_list = offsprings_df.columns.values.tolist()
        temp_parents_df = offsprings_df.loc[(offsprings_df['population'] == pop)]
        index_list_f = fathers_df.index.values.tolist()
        index_list_m = mothers_df.index.values.tolist()
        list_for_pairs = []
        n = copy.deepcopy(config["number_of_pairs_parents"])
        w = 0
        while w != n:
            father = random.choice(index_list_f)
            mother = random.choice(index_list_m)
            if father not in list_for_pairs and mother not in list_for_pairs and father != mother:
                list_for_pairs.append(father)
                list_for_pairs.append(mother)
                pair_df = temp_parents_df.loc[[father, mother]]    # create df for one random pair
                kids_df = pd.DataFrame(columns=columns_list, index=[i for i in range(0, config["kids"])])
                for k in range(config["kids"]):
                    kids_df.iloc[k]["Status"] = "kid"
                    kids_df.iloc[k]["№"] = (str(pair_df.iloc[0]['№']) + "-" + str(pair_df.iloc[1]['№']) + "-" +
                                            str(k + 1))
                    gf1 = list(kids_df.iloc[k]["№"].split("-"))[0]
                    gm1 = list(kids_df.iloc[k]["№"].split("-"))[1]
                    gf2 = list(kids_df.iloc[k]["№"].split("-"))[3]
                    gm2 = list(kids_df.iloc[k]["№"].split("-"))[4]
                    f = list(kids_df.iloc[k]["№"].split("-1-"))[0] + "-1"
                    m = list(kids_df.iloc[k]["№"].split("-1-"))[1] + "-1"
                    offsprings_df.loc[offsprings_df['№'] == gf1, 'Family_ID'] = w + 1
                    offsprings_df.loc[offsprings_df['№'] == gm1, 'Family_ID'] = w + 1
                    offsprings_df.loc[offsprings_df['№'] == gf2, 'Family_ID'] = w + 1
                    offsprings_df.loc[offsprings_df['№'] == gm2, 'Family_ID'] = w + 1
                    offsprings_df.loc[offsprings_df['№'] == f, 'Family_ID'] = w + 1
                    offsprings_df.loc[offsprings_df['№'] == m, 'Family_ID'] = w + 1
                    kids_df.iloc[k]["Family_ID"] = w + 1
                    kids_df.iloc[k]["population"] = pop
                    kids_df.iloc[k]["№"] = (str(pair_df.iloc[0]['№']) + "-" + str(pair_df.iloc[1]['№']) + "-" +
                                            str(k + 1))
                    for elem in range(1, len(columns_list) - 4, 2):
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
                    offsprings_df = pd.concat([offsprings_df, kids_df.loc[[k]]], ignore_index=True)
                w += 1
        offsprings_df.to_excel("generated_families_GRAND_240329.xlsx", index=False)

        ######################
        # Count LR
        lr_all_df = offsprings_df.loc[(offsprings_df['Status'] != "father")]
        lr_all_kids_grands_df = lr_all_df.loc[(offsprings_df['Status'] != "mother")]
        kidd = lr_all_kids_grands_df.loc[(offsprings_df['Status'] == "kid")]
        kids_index = kidd.index.values.tolist()
        grandparents = lr_all_kids_grands_df.loc[(offsprings_df['Status'] != "kid")]
        grandparents_index = grandparents.index.values.tolist()
        kids_list = lr_all_kids_grands_df.loc[(offsprings_df['Status'] == "kid")]['№'].tolist()
        grandparents_list = lr_all_kids_grands_df.loc[(offsprings_df['Status'] != "kid")]['№'].tolist()
        lr_old_codis_ref_df = pd.DataFrame(columns=grandparents_list, index=kids_list)
        lr_new_codis_ref_df = pd.DataFrame(columns=grandparents_list, index=kids_list)
        lr_old_codis_rus_df = pd.DataFrame(columns=grandparents_list, index=kids_list)
        lr_new_codis_rus_df = pd.DataFrame(columns=grandparents_list, index=kids_list)
        count_el = 0
        for el in kids_index:
            prob_df = kidd.loc[[el]]
            count_k = 0
            for kk in grandparents_index:
                hyp1_ref_old = []
                hyp2_ref_old = []
                hyp1_rus_old = []
                hyp2_rus_old = []
                hyp1_ref_new = []
                hyp2_ref_new = []
                hyp1_rus_new = []
                hyp2_rus_new = []
                gr_df = grandparents.loc[[kk]]
                for loc in range(1, len(columns_list) - 3, 2):
                    el = columns_list[loc].split('_')[0]
                    kid = [prob_df.iloc[0][columns_list[loc]], prob_df.iloc[0][columns_list[loc + 1]]]
                    gr = [gr_df.iloc[0][columns_list[loc]], gr_df.iloc[0][columns_list[loc + 1]]]
                    kid.sort()
                    gr.sort()
                    p1, p2, p3, p4 = hyp_grandparents(gr[0], gr[1], kid[0], kid[1], ref_dict, config["rus_dict"], el)
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
                lr_old_codis_ref_df.iat[count_el, count_k] = multiplication_1 / multiplication_2
                lr_new_codis_ref_df.iat[count_el, count_k] = multiplication_3 / multiplication_4
                lr_old_codis_rus_df.iat[count_el, count_k] = multiplication_5 / multiplication_6
                lr_new_codis_rus_df.iat[count_el, count_k] = multiplication_7 / multiplication_8
                count_k += 1
            count_el += 1
        lr_old_codis_ref_df.to_excel("generated_lr_old_codis_ref_df_240329.xlsx", index=True)
        lr_new_codis_ref_df.to_excel("generated_lr_new_codis_ref_df_240329.xlsx", index=True)
        lr_old_codis_rus_df.to_excel("generated_lr_old_codis_rus_df_240329.xlsx", index=True)
        lr_new_codis_rus_df.to_excel("generated_lr_new_codis_rus_df_240329.xlsx", index=True)
    print(round(time.time() - start, 2), 's')


if __name__ == "__main__":
    main()
