This repo contains four strategies collected by Odin, and used for training Delete model,  **D**eep **Le**ad Op**t**imization **E**nveloped in Protein Pocket through Unified Deleting Strategies and a Structure-aware Network. 

I provide the source code here, hope this collection could help you!

A detailed explanation could be found at `fragmentation.ipynb`

In `lead_decomp.py`, I provide a more comprehensive way to decompose molecules for linker design. fragment elaboration, scaffold hopping, and side-chain decoration. 

```python
from lead_decomp import linker_decomp, fragment_decomp, scaffold_decompo, side_chains_decomp

ligand_nm = './1djy_A_rec_1djz_ip2_lig_tt_min_0.sdf'
mol = read_sdf(ligand_nm)[0]

linker_decomp_infos = linker_decomp(mol)
fragment_decomp_infos = fragment_decomp(mol)
scaffold_decomp_infos = scaffold_decompo(mol)
side_chains_decomp_infos = side_chains_decomp(mol)
```

