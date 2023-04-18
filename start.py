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
    # Creating and adding the console handler
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

def check_mutagenesis(pos, pdb_path, res):
    args = subprocess.Popen(['pdb_selres', f'-{pos}', f'{pdb_path}'], stdout=subprocess.PIPE)
    args_ = subprocess.Popen(['grep', 'ATOM'], stdin=args.stdout, stdout=subprocess.PIPE)
    args__ = subprocess.Popen(['tail', '-1'], stdout=subprocess.PIPE, stdin=args_.stdout)
    out, _ = args__.communicate()
    return res == out.split()[3]


def main():
    parser = argparse.ArgumentParser(description="run")
    parser.add_argument("--dataset", type=str, default="dataset_prot.csv", help="dataset")
    parser.add_argument("--delim", type=str, default=";", help="delimiter of each line in the dataset file")
    parser.add_argument("--rosetta_path", type=str, default='/home/ekaterina_andreichuk/Downloads/rosetta_src_2021.16.61629_bundle/main/source/bin/rosetta_scripts.default.linuxgccrelease', help="path to rosetta")
    parser.add_argument("--gmx_path", type=str, default='/usr/local/bin/gmx', help="path to gmx")
    #add force reload from some step for particular pdb_id
    # parser.add_argument("--lambda_u", type=float, default=0.1, help="regularization parameter lambda_u of the CF model")
    # parser.add_argument("--lambda_v", type=float, default=0.1, help="regularization parameter lambda_v of the CF model")
    # parser.add_argument("--als_iter", type=int, default=15, help="# of iterations for ALS training")
    # parser.add_argument("--debug_iter", type=int, default=20, help="# of iterations in the debugging stage")
    # parser.add_argument("--debug_lr", type=float, default=0.05, help="learning rate in the debugging stage")
    # parser.add_argument("--retrain", type=str, default="full", help="the retraining mode in the debugging stage: full/inc")
    # parser.add_argument("--process", type=int, default=4, help="# of processes in the debugging stage")
    # parser.add_argument("--mode", type=str, default="debug", help="debug/test")
    # parser.add_argument("--implicit", action='store_true', help="use implicit ALS")
    # parser.add_argument("--alpha", type=int, default=1, help="confidence scaling for implicit feedback dataset")
    # parser.add_argument("--als_threads", type=int, default=6, help="num threads during implicit ALS fit")

    args_p = parser.parse_args()

    cur_dir = os.getcwd()
    df = pd.read_csv(args_p.dataset, delimiter =args_p.delim)
    without_structure = []
    skipped = []
    logger_all = my_custom_logger(os.path.join(cur_dir, 'logger.log'))
    for _, row in df.iterrows():
        pdb_id = row.pdb_id
        position = row.position
        wild_type = row.wild_type
        mutation = row.mutation 
        if os.path.isdir(os.path.join(cur_dir, pdb_id)) and pdb_id not in without_structure:
            #PP: check that each 'if' has 'else'
            #PP: avoid jumping over directories. You can make path variable using os.path.join(workdir, pdb_id)
            if not os.path.isfile(os.path.join(cur_dir, pdb_id, f'{pdb_id}_clean.pdb')):   
                if not clean_pdb(cur_dir, pdb_id):
                    without_structure.append(pdb_id)
                    logger_all.info(f'Error during cleaning pdb file for {pdb_id}. Continuing to the next structure.')
                else:
                    logger_all.info(f'Cleaned {pdb_id}.pdb')

            dir_name = wild_type+str(position)+mutation
            mut_dir = os.path.join(cur_dir, pdb_id, dir_name)
            if not os.path.isdir(os.path.join(cur_dir, pdb_id, dir_name)):        
                os.mkdir(mut_dir) 
            if os.path.isfile(os.path.join(mut_dir, 'md.gro')):
                skipped.append(' '.join([pdb_id, dir_name]))
                continue
            logger = my_custom_logger(os.path.join(mut_dir, 'logger.log'))          
            logger.info(os.path.join(pdb_id, dir_name))
            if not os.path.isfile(os.path.join(mut_dir, f'{dir_name}.pdb')):
                #create mutate.xml
                template_string = f'''<MutateResidue name="one_point_mutation" target="{str(position)}A" new_res="{acids[mutation].upper()}" preserve_atom_coords="false"
                mutate_self="false" update_polymer_bond_dependent="false"/>'''
                temp_string = '<Add mover="one_point_mutation"/>'
                new_text = []
                with open(os.path.join(cur_dir, 'mutate.xml'), 'r') as f:
                    for line in f:
                        if '/MOVERS' in line:
                            new_text.append(template_string)
                        elif '</PROTOCOLS>' in line:
                            new_text.append(temp_string)
                        new_text.append(line)
                with open(os.path.join(mut_dir, 'mutate.xml'), 'w') as f:
                    for i in new_text:
                        f.write(i)
                args = [args_p.rosetta_path, '-s', os.path.join(cur_dir, pdb_id, f'{pdb_id}_clean.pdb'), '-parser:protocol', 
                        os.path.join(mut_dir, 'mutate.xml')] 
                subprocess.call(args, stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)
                args = ['mv', os.path.join(mut_dir, f'{pdb_id}_clean_0001.pdb'), os.path.join(mut_dir, f'{dir_name}.pdb')]
                subprocess.call(args, stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)
                if check_mutagenesis(str(position), os.path.join(mut_dir), mutation):
                    logger.info(f'{dir_name}.pdb is sucessfully created.')
                else:
                    logger.warning(f'Mutagenesis for {dir_name}.pdb was unsucessull, continuing to the next mutation.')
                    logger_all.warning(f'{dir_name}.  Error at Mutagenesis stage. Continuing to next mutation.')
                    continue
            else:
                logger.info(f'{dir_name}.pdb already exists. Mutagenesis step is skipped.')
            if not os.path.isfile(os.path.join(mut_dir,'em.gro')):
                #PP: what is emp.gro?

                args = subprocess.Popen([ 'printf', '6\n1\n'], stdout=subprocess.PIPE)
                args_ = subprocess.Popen([args_p.gmx_path, 'pdb2gmx', '-f', os.path.join(mut_dir,f'{dir_name}.pdb'), '-o', os.path.join(mut_dir, 'pep_.gro'), '-p', os.path.join(mut_dir, 'pep.top'), '-ignh', '-q'], 
                                        stdin = args.stdout, stdout= subprocess.DEVNULL)#, stderr = subprocess.DEVNULL)
                out, err = args_.communicate()
                print(err)

                #PP this subprocess will be successful, even if gmx returns error, so check that gmx worked properly, i.e. all required files were created.
                args = [args_p.gmx_path, 'editconf', '-f', os.path.join(mut_dir, 'pep_.gro'), '-o', os.path.join(mut_dir,'pep.gro'), '-d', '0.3', '-bt', 'cubic',  '-c']
                subprocess.call(args, stdout = subprocess.DEVNULL)#, stderr = subprocess.DEVNULL)
                if not os.path.isfile(os.path.join(mut_dir, 'pep.gro')):
                    logger.warning(f'{pdb_id}, {dir_name}.There was an error attempting fitting a cubic box. Continuing to the next mutation.')
                    logger_all.warning(f'{pdb_id}, {dir_name}.  Error at creating system stage. Continuing to next mutation.')
                    continue
                #PP: what is sys.top?
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

                args = [args_p.gmx_path, 'grompp', '-f', os.path.join(cur_dir, 'em.mdp'), '-p', os.path.join(mut_dir, 'sys.top'), '-o', os.path.join(mut_dir, 'emp.tpr'), '-c', os.path.join(mut_dir, 'pep.gro')]
                subprocess.call(args, stdout = subprocess.DEVNUL)#L, stderr=subprocess.DEVNULL)
                args = [args_p.gmx_path, 'mdrun', '-deffnm', os.path.join(mut_dir, 'emp'), '-v']
                subprocess.call(args, stdout = subprocess.DEVNULL)#, stderr=subprocess.DEVNULL)
                if not os.path.isfile(os.path.join(mut_dir, 'emp.tpr')) or not os.path.isfile(os.path.join(mut_dir, 'emp.trr')):
                    logger.warning(f'{pdb_id}, {dir_name} Error at first minimization stage. Continuing to next mutation.')
                    logger_all.warning(f'{pdb_id}, {dir_name}.  Error at frist minimization stage. Continuing to next mutation.')
                    continue
                else:
                    logger.info(f'{pdb_id}, {dir_name}. Sucessfully ran first minimization step.')
                args = [args_p.gmx_path, 'solvate', '-cp', os.path.join(mut_dir, 'emp.gro'), '-cs', '-o', os.path.join(mut_dir,'s0.gro'), '-p', os.path.join(mut_dir, 'sys.top')]
                subprocess.call(args, stdout = subprocess.DEVNULL)#, stderr=subprocess.DEVNULL)
                args = [args_p.gmx_path, 'grompp', '-f', os.path.join(cur_dir, 'em.mdp'), '-c', os.path.join(mut_dir, 's0.gro'), '-p', os.path.join(mut_dir, 'sys.top'), '-o', os.path.join(mut_dir, 'ion.tpr')]
                subprocess.call(args, stdout = subprocess.DEVNULL)#, stderr=subprocess.DEVNULL)
                args = subprocess.Popen([ 'echo', '13'], stdout=subprocess.PIPE)#, stderr=subprocess.DEVNULL)
                args_ = subprocess.Popen([args_p.gmx_path, 'genion', '-s', os.path.join(mut_dir, 'ion.tpr'), '-neutral', '-conc', '0.1', '-p', os.path.join(mut_dir, 'sys.top'), '-o', os.path.join(mut_dir, 'start.gro')], 
                                        stdin = args.stdout, stdout= subprocess.DEVNULL)#, stderr=subprocess.DEVNULL)
                args_.communicate()
                args = [args_p.gmx_path, 'grompp', '-f', os.path.join(cur_dir, 'em.mdp'), '-c', os.path.join(mut_dir,'start.gro'), '-p', os.path.join(mut_dir, 'sys.top'), '-o', os.path.join(mut_dir, 'em.tpr')]
                subprocess.call(args, stdout= subprocess.DEVNULL)#, stderr=subprocess.DEVNULL)
                args = [args_p.gmx_path, 'mdrun', '-deffnm', os.path.join(mut_dir, 'em'), '-v']
                subprocess.call(args, stdout= subprocess.DEVNULL)#, stderr=subprocess.DEVNULL)
                if os.path.isfile(os.path.join(mut_dir, 'em.gro')):
                    logger.info(f'{pdb_id}, {dir_name}. Sucessfully ran second minimization')
                else: 
                    logger.warning(f'{pdb_id}, {dir_name} Error at second minimization stage. Continuing to next mutation.')
                    logger_all.warning(f'{pdb_id}, {dir_name}.  Error at second minimization stage. Continuing to next mutation.')
                    continue
            else:
                logger.info(f'{pdb_id}, {dir_name} em.gro already exists. Skipping minimization step.')
            if not os.path.isfile(os.path.join(mut_dir,'eq.gro')):
                args = [args_p.gmx_path, 'grompp', '-f', os.path.join(cur_dir, 'eq.mdp'), '-c', os.path.join(mut_dir, 'em.gro'), '-p', os.path.join(mut_dir, 'sys.top'), '-o', os.path.join(mut_dir, 'eq.tpr')]
                subprocess.call(args, stdout = subprocess.DEVNULL)#, stderr= subprocess.DEVNULL)
                args = [args_p.gmx_path, 'mdrun', '-deffnm', os.path.join(mut_dir, 'eq'), '-v']
                subprocess.call(args, stdout = subprocess.DEVNULL)#, stderr=subprocess.DEVNULL)
                print(err)
                if os.path.isfile(os.path.join(mut_dir, 'eq.gro')):
                    logger.info(f'{pdb_id}, {dir_name}. Sucessfully ran equilibration')
                else: 
                    logger.warning(f'{pdb_id}, {dir_name} Error at equilibration stage. Continuing to next mutation.')
                    logger_all.warning(f'{pdb_id}, {dir_name}.  Error at equilibration stage. Continuing to next mutation.')
                    continue
            else:
                logger.info(f'{pdb_id}, {dir_name} eq.gro already exists. Skipping equilibration step.')
            if not os.path.isfile(os.path.join(mut_dir,'md.gro')):
                args = [args_p.gmx_path, 'grompp', '-f', os.path.join(cur_dir, 'md.mdp'), '-c', os.path.join(mut_dir, 'eq.gro'), '-p', os.path.join(mut_dir, 'sys.top'), '-o', os.path.join(mut_dir, 'md.tpr')]
                subprocess.call(args, stdout = subprocess.DEVNULL)#, stderr= subprocess.DEVNULL)
                args= [args_p.gmx_path, 'mdrun', '-deffnm', os.path.join(mut_dir, 'md')]
                subprocess.call(args, stdout = subprocess.DEVNULL)#, stderr=subprocess.DEVNULL)
                if os.path.isfile(os.path.join(mut_dir, 'md.gro')):
                    logger.info(f'{pdb_id}, {dir_name}. Sucessfully ran simulation run')
                    logger_all.info(f'{pdb_id}, {dir_name}. Sucessfully ran simulation run')
                else: 
                    logger.warning(f'{pdb_id}, {dir_name} Error at simulation run. Continuing to next mutation.')
                    continue
        elif pdb_id not in without_structure:
            without_structure.append(pdb_id)
            logger_all.warning(f'There are no structure for {pdb_id}. To run simulations create corresponding directory with {pdb_id} name and load {pdb_id}.pdb. Continuing to the next structure')
    logger_all.info(f"Skipped following mutations as they are already finished: {', '.join(skipped)}")

if __name__ =='__main__':
    main()