import pandas as pd
import random
import time
import yaml

def empty_freq_table(path):
    """
    Loads data from a specified file and creates an empty dictionary for allele frequencies.

    :param path: The file path from which the data will be read to initialize the dictionary.
    :type path: str

    :return: An empty dictionary intended for loading allele frequencies.
    :rtype: dict
    """
    freq_df = pd.read_excel(path)
    f_columns = freq_df.columns.values.tolist()
    key_dict = []
    for el in range(1, len(f_columns) - 2, 2):
        key_dict.append(list(f_columns[el].split("_"))[0])
    freq_dict = dict.fromkeys(key_dict)
    return freq_dict

def calculate_frequencies(d, path):
    """
    Loads data from a file, calculates STR frequencies, and fills the dictionary with these frequencies.

    :param d: An empty dictionary intended for storing STR allele frequencies.
    :type d: dict

    :param path: The file path from which the data will be read to calculate and populate the dictionary.
    :type path: str

    :return: A dictionary populated with STR allele frequencies.
    :rtype: dict
    """
    freq_df = pd.read_excel(path)
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
    """
    The main function of the program.

    Generates a table with the genotypes of individuals based on STR allele frequencies and saves it as an .xlsx file.
    """

    start = time.time()
    with open("./config_file.yaml", "r") as yaml_file:
        config = yaml.load(yaml_file, Loader=yaml.FullLoader)
    ref_dict = empty_freq_table(config["main_table_path"])
    ref_dict = calculate_frequencies(ref_dict, config["main_table_path"])
    temp_df = pd.read_excel(config["main_table_path"])
    columns_list = temp_df.columns.values.tolist()
    n_parents = config["number_of_gen_individuals"]
    final_df = pd.DataFrame(columns=columns_list)
    for pp in range(len(config["populations"])):
        parents_df = pd.DataFrame(columns=columns_list, index=range(0, n_parents))
        pop_name = config["populations"][pp]
        pop_name_short = config["names"][pop_name]
        for m in range(round(n_parents/2)):    # m - male
            parents_df.iloc[m]["№"] = pop_name_short + "m" + str(m + 1)
            parents_df.iloc[m]["population"] = pop_name
        for f in range(round(n_parents/2)):    # f - female
            parents_df.iloc[f + round(n_parents/2)]["№"] = pop_name_short + "f" + str(f + 1)
            parents_df.iloc[f + round(n_parents/2)]["population"] = pop_name
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
        final_df = pd.concat([final_df, parents_df], ignore_index=True)
    final_df.to_excel(config["output_file_gen_path"], index=False)
    print(round(time.time() - start, 2), 's')

if __name__ == "__main__":
    main()
