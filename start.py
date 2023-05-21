#!/usr/bin/env python3
import numpy as np
import pandas as pd
import os
import subprocess
import shutil
from Bio.Data.IUPACData import protein_letters_1to3 as acids
import argparse
import logging
import sys
import json
from copy import copy


def my_custom_logger(logger_name, level=logging.DEBUG):
    """
    Method to return a custom logger with the given name and level
    """
    logger = logging.getLogger(logger_name)
    logger.setLevel(level)
    format_string = (
        "%(asctime)s — %(levelname)s: %(lineno)d — %(message)s"
    )
    log_format = logging.Formatter(format_string)
    # Creating and adding the console handler``
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setFormatter(log_format)
    logger.addHandler(console_handler)
    # Creating and adding the file handler
    file_handler = logging.FileHandler(logger_name, mode='w')
    file_handler.setFormatter(log_format)
    logger.addHandler(file_handler)
    return logger

def clean_pdb(cur_dir, file_name):
    # cleans water molecules and leaves only one chain (in case of a number subunits in file )
    with open(os.path.join(cur_dir, f'{file_name}.pdb')) as a:
        for line in a:
            if line.startswith('DBREF'):
                chain = line.split()[2]
                break
    args = subprocess.Popen(['pdb_selchain', f'-{chain}', f'{file_name}.pdb'], stdout=subprocess.PIPE)
    args_ = subprocess.Popen(['pdb_delhetatm'], stdin=args.stdout, stdout = subprocess.PIPE)
    args__ = subprocess.Popen(['pdb_tidy'], stdin=args_.stdout, stdout=open(f'{file_name}_clean.pdb', 'w'))
    args__.communicate()
    args = subprocess.Popen(['pdb_wc', f'{file_name}_clean.pdb'], stdout=subprocess.PIPE)
    args_ = subprocess.Popen(['grep', 'chain'], stdin=args.stdout, stdout = subprocess.PIPE)
    out, _ = args_.communicate()
    return int(out.decode().split('\t')[1]) ==1

def check_mutagenesis(pos, pdb_path, dir_name, res):
    try:
        args = subprocess.Popen(['pdb_selres', f'-{pos}', os.path.join(pdb_path, f'{dir_name}')], stdout=subprocess.PIPE)
        args_ = subprocess.Popen(['grep', 'ATOM'], stdin=args.stdout, stdout=subprocess.PIPE)
        args__ = subprocess.Popen(['tail', '-1'], stdout=subprocess.PIPE, stdin=args_.stdout)
        out, _ = args__.communicate()
        return acids[res] == out.decode().split()[3]
    except:
        return False
    
def get_mutation_dikt(df, force_reload, cur_dir):
    lz = df[['pdb_id', 'position', 'wild_type', 'mutation']].groupby('pdb_id').apply(lambda x: x[['position', 'wild_type', 'mutation']].to_dict(orient = 'records')).to_dict()
    lz = {i: [k['wild_type']+str(k['position'])+k['mutation'] for k in j] for i, j in lz.items()}
    md_dikt = copy(lz)
    mutation_dikt = copy(md_dikt)
    em_dikt = copy(md_dikt)
    eq_dikt = copy(md_dikt)
    skipped = []
    for prot in lz:
        mutation_reload = force_reload['mutation'].get(prot, [])
        em_reload = force_reload['em'].get(prot, [])
        eq_reload = force_reload['eq'].get(prot, [])
        md_reload = force_reload['md'].get(prot, [])
        for mut in lz[prot]:
            if os.path.isdir(os.path.join(cur_dir, prot, mut, f'{mut}.pdb')) and mutation_reload!='all' and mut not in mutation_reload:
                mutation_dikt[prot].remove(mut)
            else:
                continue
            if os.path.isdir(os.path.join(cur_dir, prot, mut, f'em.gro')) and em_reload!='all' and mut not in em_reload:
                em_dikt[prot].remove(mut)
            else:
                continue
            if os.path.isdir(os.path.join(cur_dir, prot, mut, f'eq.gro')) and eq_reload!='all' and mut not in eq_reload:
                eq_dikt[prot].remove(mut)
            else:
                continue
            if os.path.isdir(os.path.join(cur_dir, prot, mut, f'md.gro')) and md_reload!='all' and mut not in md_reload:
                md_dikt[prot].remove(mut)
                skipped.append(' '.join([prot, mut]))
            else:
                continue
    return mutation_dikt, em_dikt, eq_dikt, md_dikt, skipped
            

def make_rosetta_xml_for_point_mutation(cur_dir, mut_dir, dir_name, pdb_id, position, mutation, args_p):
    #create mutate.xml
    template_string = f'''<MutateResidue name="one_point_mutation" target="{str(position)}A" new_res="{acids[mutation].upper()}" preserve_atom_coords="false"
    mutate_self="false" update_polymer_bond_dependent="false"/>'''
    temp_string = '<Add mover="one_point_mutation"/>'
    new_text = []
    with open(os.path.join(cur_dir, 'resorces', 'mutate.xml'), 'r') as f:
        for line in f:
            if '/MOVERS' in line:
                new_text.append(template_string)
            elif '</PROTOCOLS>' in line:
                new_text.append(temp_string)
            new_text.append(line)
    with open(os.path.join(mut_dir, f'mutate_{dir_name}.xml'), 'w') as f:
        for i in new_text:
            f.write(i)
    args = [args_p.rosetta_path, '-s', os.path.join(cur_dir, pdb_id, f'{pdb_id}_clean.pdb'), '-parser:protocol', 
            os.path.join(mut_dir, f'mutate_{dir_name}.xml')] 
    subprocess.call(args, stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)
    args = ['mv', os.path.join(mut_dir, f'{pdb_id}_clean_0001.pdb'), os.path.join(mut_dir, f'{dir_name}.pdb')]
    subprocess.call(args, stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)

def fit_to_box(dir_name, mut_dir, args_p):
    args = subprocess.Popen([ 'printf', '6\n1\n'], stdout=subprocess.PIPE)
    args_ = subprocess.Popen([args_p.gmx_path, 'pdb2gmx', '-f', os.path.join(mut_dir,f'{dir_name}.pdb'), '-o', os.path.join(mut_dir, 'pep_.gro'), '-p', os.path.join(mut_dir, 'pep.top'), '-ignh', '-q'], 
                            stdin = args.stdout, stdout= subprocess.DEVNULL)#, stderr = subprocess.DEVNULL)
    oerr = args_.communicate()

    args = [args_p.gmx_path, 'editconf', '-f', os.path.join(mut_dir, 'pep_.gro'), '-o', os.path.join(mut_dir,'pep.gro'), '-d', '0.3', '-bt', 'cubic',  '-c']
    subprocess.call(args, stdout = subprocess.DEVNULL)#, stderr = subprocess.DEVNULL)

def prepare_topology(mut_dir, cur_dir):
    shutil.copyfile(os.path.join(cur_dir, 'sys.top'), os.path.join(mut_dir, 'sys.top')) 
    text= []
    mark = False
    with open(os.path.join(mut_dir, 'pep.top')) as file:
        for line in file:
            if 'amber' not in line and not mark:
                text.append(line)
            if 'posre.itp' in line:
                mark = True
    text.append('#endif')
    with open(os.path.join(mut_dir, 'pep.itp'), 'w') as f:
        for i in text:
            f.write(i)

def minimization_run(mut_dir, cur_dir, args_p):
    args = [args_p.gmx_path, 'grompp', '-f', os.path.join(cur_dir, 'resorces', 'mdp', 'em.mdp'), '-p', os.path.join(mut_dir, 'sys.top'), '-o', os.path.join(mut_dir, 'emp.tpr'), '-c', os.path.join(mut_dir, 'pep.gro')]
    subprocess.call(args, stdout = subprocess.DEVNUL)#L, stderr=subprocess.DEVNULL)
    args = [args_p.gmx_path, 'mdrun', '-deffnm', os.path.join(mut_dir, 'emp'), '-v']
    subprocess.call(args, stdout = subprocess.DEVNULL)#, stderr=subprocess.DEVNULL)

def system_solvate(mut_dir, cur_dir, args_p):
    args = [args_p.gmx_path, 'solvate', '-cp', os.path.join(mut_dir, 'emp.gro'), '-cs', '-o', os.path.join(mut_dir,'s0.gro'), '-p', os.path.join(mut_dir, 'sys.top')]
    subprocess.call(args, stdout = subprocess.DEVNULL)#, stderr=subprocess.DEVNULL)
    args = [args_p.gmx_path, 'grompp', '-f', os.path.join(cur_dir, 'resorces', 'mdp',  'em.mdp'), '-c', os.path.join(mut_dir, 's0.gro'), '-p', os.path.join(mut_dir, 'sys.top'), '-o', os.path.join(mut_dir, 'ion.tpr')]
    subprocess.call(args, stdout = subprocess.DEVNULL)#, stderr=subprocess.DEVNULL)
    args = subprocess.Popen([ 'echo', '13'], stdout=subprocess.PIPE)#, stderr=subprocess.DEVNULL)
    args_ = subprocess.Popen([args_p.gmx_path, 'genion', '-s', os.path.join(mut_dir, 'ion.tpr'), '-neutral', '-conc', '0.1', '-p', os.path.join(mut_dir, 'sys.top'), '-o', os.path.join(mut_dir, 'start.gro')], 
                            stdin = args.stdout, stdout= subprocess.DEVNULL)#, stderr=subprocess.DEVNULL)
    args_.communicate()
    args = [args_p.gmx_path, 'grompp', '-f', os.path.join(cur_dir, 'resorces', 'mdp', 'em.mdp'), '-c', os.path.join(mut_dir,'start.gro'), '-p', os.path.join(mut_dir, 'sys.top'), '-o', os.path.join(mut_dir, 'em.tpr')]
    subprocess.call(args, stdout= subprocess.DEVNULL)#, stderr=subprocess.DEVNULL)
    args = [args_p.gmx_path, 'mdrun', '-deffnm', os.path.join(mut_dir, 'em'), '-v']
    subprocess.call(args, stdout= subprocess.DEVNULL)#, stderr=subprocess.DEVNULL)

def equilibration(mut_dir, cur_dir, args_p):
    args = [args_p.gmx_path, 'grompp', '-f', os.path.join(cur_dir, 'resorces', 'mdp', 'eq.mdp'), '-c', os.path.join(mut_dir, 'em.gro'), '-p', os.path.join(mut_dir, 'sys.top'), '-o', os.path.join(mut_dir, 'eq.tpr')]
    subprocess.call(args, stdout = subprocess.DEVNULL)#, stderr= subprocess.DEVNULL)
    args = [args_p.gmx_path, 'mdrun', '-deffnm', os.path.join(mut_dir, 'eq'), '-v']
    subprocess.call(args, stdout = subprocess.DEVNULL)#, stderr=subprocess.DEVNULL)

def run_simulation(cur_dir, mut_dir, args_p):
    args = [args_p.gmx_path, 'grompp', '-f', os.path.join(cur_dir, 'resorces', 'mdp', 'md.mdp'), '-c', os.path.join(mut_dir, 'eq.gro'), '-p', os.path.join(mut_dir, 'sys.top'), '-o', os.path.join(mut_dir, 'md.tpr')]
    subprocess.call(args, stdout = subprocess.DEVNULL)#, stderr= subprocess.DEVNULL)
    args= [args_p.gmx_path, 'mdrun', '-deffnm', os.path.join(mut_dir, 'md')]
    subprocess.call(args, stdout = subprocess.DEVNULL)#, stderr=subprocess.DEVNULL)


def process_point_mutation(dir_name, cur_dir, pdb_id, logger_all, mutation_dikt, em_dikt, eq_dikt, md_dikt, args_p):
    position = dir_name[1:-1]
    wild_type = dir_name[0]
    mutation = dir_name[-1]
    mut_dir = os.path.join(cur_dir, pdb_id, dir_name)
    if not os.path.isdir(os.path.join(cur_dir, pdb_id, dir_name)):        
        os.mkdir(mut_dir)
    logger = my_custom_logger(os.path.join(mut_dir, 'logger.log'))          
    logger_all.info(os.path.join(pdb_id, dir_name))
    if dir_name in mutation_dikt[pdb_id]:
        make_rosetta_xml_for_point_mutation(cur_dir, mut_dir, dir_name, pdb_id, position, mutation, args_p)
        if check_mutagenesis(position, os.path.join(mut_dir), mutation):
            logger.info(f'{dir_name}.pdb is sucessfully created.')
        else:
            logger.warning(f'Mutagenesis for {dir_name}.pdb was unsucessfull, continuing to the next mutation.')
            args = ['rm', os.path.join(mut_dir, f'{dir_name}.pdb')]
            subprocess.call(args)
            logger_all.warning(f'{pdb_id}, {dir_name}.  Error at Mutagenesis stage. Continuing to next mutation.')
            return
    else:
        logger.info(f'{dir_name}.pdb already exists. Mutagenesis step is skipped.')

    if dir_name in em_dikt[pdb_id]:
        fit_to_box(dir_name, mut_dir, args_p)
        if not os.path.isfile(os.path.join(mut_dir, 'pep.gro')):
            logger.warning(f'{pdb_id}, {dir_name}.There was an error attempting fitting a cubic box. Continuing to the next mutation.')
            logger_all.warning(f'{pdb_id}, {dir_name}.  Error at creating system stage. Continuing to next mutation.')
            return
        prepare_topology(mut_dir, cur_dir)
        #PP: btw check if gromacs already has python interface
        minimization_run(mut_dir, cur_dir, args_p)
        if not os.path.isfile(os.path.join(mut_dir, 'emp.tpr')):
            logger.warning(f'{pdb_id}, {dir_name} Error at first minimization stage. Continuing to next mutation.')
            logger_all.warning(f'{pdb_id}, {dir_name}.  Error at frist minimization stage. Continuing to next mutation.')
            return
        else:
            logger.info(f'{pdb_id}, {dir_name}. Sucessfully ran first minimization step.')
        system_solvate(mut_dir, cur_dir, args_p)
        if not os.path.isfile(os.path.join(mut_dir, 'em.gro')):
            logger.warning(f'{pdb_id}, {dir_name} Error at second minimization stage. Continuing to next mutation.')
            logger_all.warning(f'{pdb_id}, {dir_name}.  Error at second minimization stage. Continuing to next mutation.')
            return 
        else: 
            logger.info(f'{pdb_id}, {dir_name}. Sucessfully ran second minimization')
    else:
        logger.info(f'{pdb_id}, {dir_name} em.gro already exists. Skipping minimization step.')

    if dir_name in eq_dikt[pdb_id]:
        equilibration(mut_dir, cur_dir, args_p)
        if not os.path.isfile(os.path.join(mut_dir, 'eq.gro')):
            logger.warning(f'{pdb_id}, {dir_name} Error at equilibration stage. Continuing to next mutation.')
            logger_all.warning(f'{pdb_id}, {dir_name}.  Error at equilibration stage. Continuing to next mutation.')
            return
        else: 
            logger.info(f'{pdb_id}, {dir_name}. Sucessfully ran equilibration')

    else:
        logger.info(f'{pdb_id}, {dir_name} eq.gro already exists. Skipping equilibration step.')

    if dir_name in md_dikt[pdb_id]:
        run_simulation(cur_dir, mut_dir, args_p)
        if not os.path.isfile(os.path.join(mut_dir, 'md.gro')):
            logger.warning(f'{pdb_id}, {dir_name} Error at simulation run. Continuing to next mutation.')
            logger_all.warning(f'{pdb_id}, {dir_name} Error at simulation run. Continuing to next mutation.')
            return
        else: 
            logger.info(f'{pdb_id}, {dir_name}. Sucessfully ran simulation run')
            logger_all.info(f'{pdb_id}, {dir_name}. Sucessfully ran simulation run')


def main():
    parser = argparse.ArgumentParser(description="run")
    parser.add_argument("--dataset", type=str, default="dataset_prot.csv", help="dataset")
    parser.add_argument("--delim", type=str, default=";", help="delimiter of each line in the dataset file")
    parser.add_argument("--rosetta_path", type=str, default='/home/ekaterina_andreichuk/Downloads/rosetta_src_2021.16.61629_bundle/main/source/bin/rosetta_scripts.default.linuxgccrelease', help="path to rosetta")
    parser.add_argument("--gmx_path", type=str, default='/usr/local/bin/gmx', help="path to gmx")
    parser.add_argument("--recalc", type=str, default='./force_reload.json', help="path to configs for reload")

    args_p = parser.parse_args()
    force_reload = json.load(open('./force_reload.json'))
    cur_dir = os.getcwd()
    df = pd.read_csv(args_p.dataset, delimiter =args_p.delim)

    mutation_dikt, em_dikt, eq_dikt, md_dikt, skipped = get_mutation_dikt(df, force_reload, cur_dir)
    logger_all = my_custom_logger(os.path.join(cur_dir, 'logger.log'))
    logger_all.info(f"Following mutations will be skipped: {', '.join(skipped)}")
    for pdb_id in md_dikt:
        if not os.path.isdir(os.path.join(cur_dir, pdb_id)):
            logger_all.warning(f'There are no structure for {pdb_id}. To run simulations create corresponding directory with {pdb_id} name and load {pdb_id}.pdb. Continuing to the next structure')
            continue
        else:
            if not os.path.isfile(os.path.join(cur_dir, pdb_id, f'{pdb_id}_clean.pdb')):   
                if not clean_pdb(cur_dir, pdb_id):
                    logger_all.info(f'Error during cleaning pdb file for {pdb_id}. Continuing to the next structure.')
                    continue
                else:
                    logger_all.info(f'Cleaned {pdb_id}.pdb')
            else:
                logger_all.info(f'Skipping cleaning step for {pdb_id} as it already exists.')
        for dir_name in md_dikt[pdb_id]:
            process_point_mutation(dir_name, cur_dir, pdb_id, logger_all, mutation_dikt, em_dikt, eq_dikt, md_dikt, args_p)



if __name__ =='__main__':
    main()
