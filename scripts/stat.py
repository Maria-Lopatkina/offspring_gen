from scipy.stats import chisquare
from scipy.stats import chi2_contingency
import pandas as pd
import time
import yaml


def empty_freq_table(path):  # Create empty dictionaries for allele frequencies
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

def calculate_frequencies_for_gen_data (d, table):  # Calculate alleles frequencies for generated population
    f_columns = table.columns.values.tolist()
    for al in range(2, len(f_columns) - 33, 2):
        allele_name = f_columns[al].split('_')[0]
        allele_df_1 = table.iloc[:, [al]]
        allele_df_1 = allele_df_1.set_axis([allele_name], axis=1)
        allele_df_2 = table.iloc[:, [al + 1]]
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


def main():
    start = time.time()
    with open("config_file.yaml", "r") as yaml_file:
        config = yaml.load(yaml_file, Loader=yaml.FullLoader)
    # Create dataframe with offsprings
    obs_freq_df = pd.read_excel("freq_output.xlsx")
    for pop in config["populations"]:
        ref_dict = empty_freq_table(config["main_table_path"])
        ref_dict = calculate_frequencies(pop, ref_dict, config["main_table_path"])
        temp_obs_df = obs_freq_df.loc[(obs_freq_df['population'] == pop) & (obs_freq_df['Status'] == "kid")]
        obs_dict = empty_freq_table(config["main_table_path"])
        obs_dict = calculate_frequencies_for_gen_data(obs_dict, temp_obs_df)
        all_p_val = []
        for locus in obs_dict.keys():
            exp = []
            obs = []
            for repeat in ref_dict[locus].keys():
                exp.append(ref_dict[locus][repeat] * config["number_of_pairs"] * config["kids"])
                obs.append(obs_dict[locus][repeat] * config["number_of_pairs"] * config["kids"])
            # add check for cases of different len(dict[locus])
            p_value = chisquare(exp, obs)[1]
            significance_level = 0.05
            if p_value >= significance_level:
                all_p_val.append(0)
            else:
                all_p_val.append(1)
        # if 1 in all_p_val:


if __name__ == "__main__":
    main()
