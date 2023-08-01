import pandas as pd
import time


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


def main():
    start = time.time()
    codis_old = ["CSF1PO", "D3S1358", "D5S818", "D7S820", "D8S1179", "D13S317", "D16S539", "D18S51", "D21S11", "FGA",
                 "TH01", "TPOX", "vWA"]
    codis_new = ["CSF1PO", "D3S1358", "D5S818", "D7S820", "D8S1179", "D13S317", "D16S539", "D18S51", "D21S11", "FGA",
                 "TH01", "TPOX", "vWA", "D1S1656", "D2S441", "D2S1338", "D10S1248", "D12S391", "D19S433", "D22S1045"]
    # Upload data of trio and potential fathers
    trio_df = pd.read_excel("output_v12_trio.xlsx")
    pf_df = pd.read_excel("output_v12_p_fathers.xlsx")
    columns_list = trio_df.columns.values.tolist()
    for i in range(0, len(trio_df.index), 5):
        real_mother_df = trio_df.iloc[[i + 1]]
        for ch in range(3):
            child_df = trio_df.iloc[[i + 2 + ch]]
            temp_pf_df = pf_df[(pf_df['Child_ID'] == child_df.iloc[0]["№"])]
            for ii in range(0, len(temp_pf_df.index)):
                one_pf_df = temp_pf_df.iloc[[ii]]
                name = one_pf_df.iloc[0]['№'].split("_")
                count_of_mismatches = 0
                if ("old" in name) and ("trio" in name):
                    for loc in range(1, len(columns_list) - 27, 2):
                        el = columns_list[loc].split('_')[0]
                        if el in codis_old:
                            pf1, pf2 = one_pf_df.iloc[0][columns_list[loc]], one_pf_df.iloc[0][columns_list[loc + 1]]
                            rm1, rm2 = real_mother_df.iloc[0][columns_list[loc]], \
                                real_mother_df.iloc[0][columns_list[loc + 1]]
                            rc1, rc2 = child_df.iloc[0][columns_list[loc]], child_df.iloc[0][columns_list[loc + 1]]
                            answer = father_trio(pf1, pf2, rm1, rm2, rc1, rc2)
                            if not answer:
                                count_of_mismatches += 1
                    print(count_of_mismatches)
                if "old" in name and "duo" in name:
                    for loc in range(1, len(columns_list) - 27, 2):
                        el = columns_list[loc].split('_')[0]
                        if el in codis_old:
                            pf1, pf2 = one_pf_df.iloc[0][columns_list[loc]], one_pf_df.iloc[0][columns_list[loc + 1]]
                            rc1, rc2 = child_df.iloc[0][columns_list[loc]], child_df.iloc[0][columns_list[loc + 1]]
                            answer = father_duo(pf1, pf2, rc1, rc2)
                            if not answer:
                                count_of_mismatches += 1
                    print(count_of_mismatches)
                    # one_pf_df.iloc[0]["p_father_PP_rus"] = count_of_mismatches
                    # print(one_pf_df)
                    # p_f_output_df = pd.concat([p_f_output_df, one_pf_df])
                if "new" in name and "trio" in name:
                    for loc in range(1, len(columns_list) - 27, 2):
                        el = columns_list[loc].split('_')[0]
                        if el in codis_new:
                            pf1, pf2 = one_pf_df.iloc[0][columns_list[loc]], one_pf_df.iloc[0][columns_list[loc + 1]]
                            rm1, rm2 = real_mother_df.iloc[0][columns_list[loc]], \
                                real_mother_df.iloc[0][columns_list[loc + 1]]
                            rc1, rc2 = child_df.iloc[0][columns_list[loc]], child_df.iloc[0][columns_list[loc + 1]]
                            answer = father_trio(pf1, pf2, rm1, rm2, rc1, rc2)
                            if not answer:
                                count_of_mismatches += 1
                    print(count_of_mismatches)
                if "new" in name and "duo" in name:
                    for loc in range(1, len(columns_list) - 27, 2):
                        el = columns_list[loc].split('_')[0]
                        if el in codis_new:
                            pf1, pf2 = one_pf_df.iloc[0][columns_list[loc]], one_pf_df.iloc[0][columns_list[loc + 1]]
                            rc1, rc2 = child_df.iloc[0][columns_list[loc]], child_df.iloc[0][columns_list[loc + 1]]
                            answer = father_duo(pf1, pf2, rc1, rc2)
                            if not answer:
                                count_of_mismatches += 1
                    print(count_of_mismatches)
    print(round(time.time() - start, 2), 's')


if __name__ == "__main__":
    main()
