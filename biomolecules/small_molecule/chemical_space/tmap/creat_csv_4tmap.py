import pandas as pd
import argparse


def read_txt(file):
    with open(file, "r") as f:
        smiles_list = [line.strip() for line in f]
    return smiles_list

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--smiles_txt', type=str, default='used_mols.txt')
    parser.add_argument('--save_file', type=str, default='saved.csv')

    args = parser.parse_args()

    smiles_list = read_txt(args.smiles_txt)

    # Creating a DataFrame
    df = pd.DataFrame({
        'DrugBank ID': ['DB00006'] * len(smiles_list),
        'Name': ['Not Provided'] * len(smiles_list),
        'CAS Number': ['000-00-0'] * len(smiles_list),
        'Drug Groups': ['approved; investigational'] * len(smiles_list),
        'InChIKey': ['OIRCOABEOLEUMC-GEJPAHFPSA-N'] * len(smiles_list),
        'InChI': ['InChI=1S/C98H138N24O33/c1-5-52(4)82(96(153)122'] * len(smiles_list),
        'SMILES': smiles_list
    })

    # Saving the DataFrame to CSV
    df.to_csv(args.save_file, index=False)

