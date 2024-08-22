import pandas as pd
import random
import time
import yaml


def empty_freq_table(path):  # Create empty dictionaries for autosomal str allele frequencies
    freq_df = pd.read_excel(path)
    f_columns = freq_df.columns.values.tolist()
    key_dict = []
    for el in range(1, len(f_columns) - 2, 2):
        key_dict.append(list(f_columns[el].split("_"))[0])
    freq_dict = dict.fromkeys(key_dict)
    return freq_dict


def calculate_frequencies(d, path, gr):  # Calculate autosomal str allele frequencies
    freq_df = pd.read_excel(path)
    freq_df = freq_df.loc[(freq_df['groups'] == gr)]
    f_columns = freq_df.columns.values.tolist()
    for al in range(1, len(f_columns) - 2, 2):
        allele_name = f_columns[al].split('_')[0]
        allele_df_1 = freq_df.iloc[:, [al]]
        allele_df_1 = allele_df_1.set_axis([allele_name], axis=1)
        allele_df_2 = freq_df.iloc[:, [al + 1]]
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


def main():
    start = time.time()
    with open("../config_file_x_str.yaml", "r") as yaml_file:
        config = yaml.load(yaml_file, Loader=yaml.FullLoader)
    ref_dict = empty_freq_table(config["main_table_path"])
    ref_dict = calculate_frequencies(ref_dict, config["main_table_path"], config["group"])
    n_parents = config["number_of_parents"]
    temp_df = pd.read_excel(config["main_table_path"])
    columns_list = temp_df.columns.values.tolist()
    columns_list = columns_list[:-2]
    columns_list.extend(config["add_columns"])
    parents_df = pd.DataFrame(columns=columns_list, index=range(0, n_parents))
    for f in range(round(n_parents/2)):
        parents_df.iloc[f]["№"] = "f" + str(f + 1)
    for m in range(round(n_parents/2)):
        parents_df.iloc[m + round(n_parents/2)]["№"] = "m" + str(m + 1)
    for el in ref_dict.keys():
        s = 0
        for i in ref_dict[el].keys():
            ref_dict[el][i] = round(ref_dict[el][i] * n_parents * 2)
            s += ref_dict[el][i]
        if s != n_parents * 2:
            while s != n_parents * 2:
                ii = random.randint(0, len(ref_dict[el]) - 1)
                if n_parents * 2 - s > 0:
                    ref_dict[el][list(ref_dict[el].keys())[ii]] += 1
                    s += 1
                else:
                    ref_dict[el][list(ref_dict[el].keys())[ii]] -= 1
                    s -= 1
    for el in ref_dict.keys():
        loc = []
        elem_1 = el + "_1"
        elem_2 = el + "_2"
        for i in ref_dict[el].keys():
            for k in range(int(round(ref_dict[el][i]))):
                loc.append(i)
        random.shuffle(loc)
        for k in range(n_parents):
            parents_df.iloc[k][elem_1] = loc[-1]
            loc.pop()
            parents_df.iloc[k][elem_2] = loc[-1]
            loc.pop()

    # Add X-STR data for women
    bur_dict_m = config["bur_kurumkan_f"]
    for el in bur_dict_m.keys():
        s = 0
        for i in bur_dict_m[el].keys():
            bur_dict_m[el][i] = round(bur_dict_m[el][i] * n_parents)
            s += bur_dict_m[el][i]
        if s != n_parents:
            while s != n_parents:
                ii = random.randint(0, len(bur_dict_m[el]) - 1)
                if n_parents - s > 0:
                    bur_dict_m[el][list(bur_dict_m[el].keys())[ii]] += 1
                    s += 1
                else:
                    bur_dict_m[el][list(bur_dict_m[el].keys())[ii]] -= 1
                    s -= 1
    for el in bur_dict_m.keys():
        loc = []
        elem_1 = el + "_1"
        elem_2 = el + "_2"
        for i in bur_dict_m[el].keys():
            for k in range(int(round(bur_dict_m[el][i]))):
                loc.append(i)
        random.shuffle(loc)
        for k in range(int(n_parents/2)):
            parents_df.iloc[k][elem_1] = loc[-1]
            loc.pop()
            parents_df.iloc[k][elem_2] = loc[-1]

    # Add X-STR data for men
    bur_dict_m = config["bur_kurumkan_m"]
    for el in bur_dict_m.keys():
        s = 0
        for i in bur_dict_m[el].keys():
            bur_dict_m[el][i] = round(bur_dict_m[el][i] * n_parents / 2)
            s += bur_dict_m[el][i]
        if s != n_parents / 2:
            while s != n_parents / 2:
                ii = random.randint(0, len(bur_dict_m[el]) - 1)
                if n_parents / 2 - s > 0:
                    bur_dict_m[el][list(bur_dict_m[el].keys())[ii]] += 1
                    s += 1
                else:
                    bur_dict_m[el][list(bur_dict_m[el].keys())[ii]] -= 1
                    s -= 1
    for el in bur_dict_m.keys():
        loc = []
        elem_1 = el + "_1"
        elem_2 = el + "_2"
        for i in bur_dict_m[el].keys():
            for k in range(int(round(bur_dict_m[el][i]))):
                loc.append(i)
        random.shuffle(loc)
        for k in range(int(n_parents/2), n_parents):
            parents_df.iloc[k][elem_1] = loc[-1]
            loc.pop()
            parents_df.iloc[k][elem_2] = "-"

    # Add data for DXS6809-DXS6789 for women
    gapl_dict = config["DXS6809-DXS6789"]
    s = 0
    for el in gapl_dict.keys():
        gapl_dict[el] = round(gapl_dict[el] * n_parents)
        s += gapl_dict[el]
    if s != n_parents:
        while s != n_parents:
            ii = random.randint(0, len(gapl_dict) - 1)
            if n_parents - s > 0:
                gapl_dict[list(gapl_dict.keys())[ii]] += 1
                s += 1
            else:
                gapl_dict[list(gapl_dict.keys())[ii]] -= 1
                s -= 1
    el11 = "DXS6809_1"
    el12 = "DXS6809_2"
    el21 = "DXS6789_1"
    el22 = "DXS6789_2"
    for k in range(int(n_parents/2)):
        a = False
        while a != True:
            rand1 = random.choice(list(gapl_dict.keys()))
            rand2 = random.choice(list(gapl_dict.keys()))
            if gapl_dict[rand1] != 0 and gapl_dict[rand2] != 0 and not (gapl_dict[rand1] == gapl_dict[rand2] == 1):
                el11_v = int(rand1.split("-")[0])
                el21_v = int(rand1.split("-")[1])
                el12_v = int(rand2.split("-")[0])
                el22_v = int(rand2.split("-")[1])
                parents_df.iloc[k][el11] = el11_v
                parents_df.iloc[k][el21] = el21_v
                parents_df.iloc[k][el12] = el12_v
                parents_df.iloc[k][el22] = el22_v
                gapl_dict[rand1] -= 1
                gapl_dict[rand2] -= 1
                a = True

    # Add data for DXS6809-DXS6789 for men
    gapl_dict = config["DXS6809-DXS6789"]
    s = 0
    for el in gapl_dict.keys():
        gapl_dict[el] = round(gapl_dict[el] * n_parents / 2)
        s += gapl_dict[el]
    if s != n_parents / 2:
        while s != n_parents / 2:
            ii = random.randint(0, len(gapl_dict) - 1)
            if n_parents / 2 - s > 0:
                gapl_dict[list(gapl_dict.keys())[ii]] += 1
                s += 1
            else:
                gapl_dict[list(gapl_dict.keys())[ii]] -= 1
                s -= 1
    el11 = "DXS6809_1"
    el12 = "DXS6809_2"
    el21 = "DXS6789_1"
    el22 = "DXS6789_2"
    for k in range(int(n_parents/2), n_parents):
        a = False
        while a != True:
            rand = random.choice(list(gapl_dict.keys()))
            if gapl_dict[rand] != 0:
                el11_v = int(rand.split("-")[0])
                el21_v = int(rand.split("-")[1])
                parents_df.iloc[k][el11] = el11_v
                parents_df.iloc[k][el21] = el21_v
                parents_df.iloc[k][el12] = "-"
                parents_df.iloc[k][el22] = "-"
                gapl_dict[rand] -= 1
                a = True
    for j in range(n_parents):
        parents_df.iloc[j]["population"] = config["population"]
        parents_df.iloc[j]["groups"] = config["group"]
    parents_df.to_excel("generated_families.xlsx", index=True)
    print(round(time.time() - start, 2), 's')


if __name__ == "__main__":
    main()
