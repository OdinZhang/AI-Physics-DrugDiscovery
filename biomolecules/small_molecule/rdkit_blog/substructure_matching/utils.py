from rdkit import Chem
from rdkit import DataStructs
from PyPDF2 import PdfMerger 
import io
from PIL import Image
import numpy as np
import os
import os.path as osp

def compute_sim(ref_mol, gen_mol):
    fp_refmol = Chem.RDKFingerprint(ref_mol)
    fp_genmol = Chem.RDKFingerprint(gen_mol)
    sim = DataStructs.TanimotoSimilarity(fp_refmol, fp_genmol)
    return sim

def compute_sims(ref_mols, gen_mols):
    num_refmols = len(ref_mols)
    num_genmols = len(gen_mols)
    sims = np.ones([num_refmols, num_genmols])
    for i in range(num_refmols):
        for j in range(num_genmols):
            sims[i][j] = compute_sim(ref_mols[i], gen_mols[j])
    return sims

def generalize(core):
    query_params = Chem.AdjustQueryParameters()
    query_params.makeBondsGeneric = True
    query_params.aromatizeIfPossible = False
    query_params.adjustDegree = False
    query_params.adjustHeavyDegree = False
    generic_core = Chem.AdjustQueryProperties(core,query_params)
    return generic_core

def find_match(target,query_mols):
    for query_mol in query_mols:
        query = generalize(query_mol)
        match = target.GetSubstructMatch(query)
        if len(match) > 0:
            return match
    return ()

def sort_lists_by_first_list(lists, ascending=True):
    if not lists or not all(lists):
        return []  # Handle empty lists or lists with empty sublists

    # Determine the sorting order
    sort_order = sorted(range(len(lists[0])), key=lists[0].__getitem__, reverse=not ascending)

    # Sort each list in lists according to the sort_order
    sorted_lists = []
    for lst in lists:
        if len(lst) != len(lists[0]):
            raise ValueError("All lists must be of the same length")
        sorted_lists.append([lst[i] for i in sort_order])

    return sorted_lists

def imgs2singlePDF(imgs, out_pdf):
    # img is Draw.MolsToGridImage
    pdf_files = []
    for idx, img in enumerate(imgs):
        img_byte_arr = io.BytesIO(img.data)
        img_pil = Image.open(img_byte_arr)
        pdf_file = f"./group_{idx}.pdf"
        img_pil.save(pdf_file)
        pdf_files.append(pdf_file)
    merger = PdfMerger()
    for pdf in pdf_files:
        merger.append(pdf)
    merger.write(out_pdf)
    merger.close()
    for pdf_file in pdf_files:
        os.remove(pdf_file)

def read_sdf(sdf_file, sanitize=False):
    supp = Chem.SDMolSupplier(sdf_file, sanitize=sanitize)
    mols_list = [i for i in supp]
    return mols_list

def write_sdf(mol_list,file, voice=False):
    writer = Chem.SDWriter(file)
    mol_cnt = 0
    for i in mol_list:
        try:
            writer.write(i)
            mol_cnt+=1
        except:
            pass
    writer.close()
    if voice: 
        print('Write {} molecules to {}'.format(mol_cnt,file))

def save_img(ipy_Image, out_file):
    img_byte_arr = io.BytesIO(ipy_Image.data)
    img_pil = Image.open(img_byte_arr)
    png_file = searched_file.replace('.sdf','.png')
    img_pil.save(out_file)
    print('saved at {}'.format(out_file))