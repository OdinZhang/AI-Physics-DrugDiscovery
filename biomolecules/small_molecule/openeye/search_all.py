import subprocess
from glob import glob 
import os.path as osp
import time

sdf_files = glob('./BRD4_delete_gen_select/*')
save_base = './BRD4_delete_gen_in_chemdiv'

for sdf_file in sdf_files:
    sdf_name = osp.basename(sdf_file)
    searched_name = sdf_name.replace('.sdf','_chemdiv.sdf')
    searched_file = osp.join(save_base, searched_name)
    command = f'python searchfastfp.py -fpdb chemdiv.fpbin -molfname chemdiv.sdf -query {sdf_file} -out {searched_file} -memorytype in-memory'
    print(command)
    result = subprocess.run(command, shell=True, capture_output=True, text=True)

    if result.returncode == 0:
        print('executed successfully.')
        print(result.stdout)
    else:
        print('execution failed.')
        print(result.stderr)
