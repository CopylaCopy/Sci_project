#!/usr/bin/env python3
import numpy as np
import pandas as pd
import os
import subprocess
import shutil




#PP: it is better to use pdbtools for cleaning the file. Particularly, there are a lot of pitfalls regarding the residue numbering.
def clean_pdb(file_name):
    new_file = []
    counter_new = 0
    counter_new_res = 0
    res_counter = 0
    mark = True
    with open(file_name+'.pdb') as a:
        for line in a:
            if line.startswith('ATOM') and mark:
                new_res = line[17:20]
                counter_atom = int(line[6:11].lstrip())
                counter_res = int(line[23:26].lstrip())
                if counter_res !=res_counter:
                    counter_new_res +=1
                    temp = '    '+ str(counter_new_res)
                    temp = temp[-5:]
                counter_new+=1
                temp_atom = '    '+str(counter_new)
                temp_atom = temp_atom[-5:]

                if res_counter<= counter_res:
                    res_counter = counter_res
                    new_line = line[:6]+temp_atom+line[11:22]+temp+line[26:]
                    new_file.append(line)
                    #break
                else:
                    mark = False
            elif not line.startswith('ATOM') and not line.startswith('HETATM'):
                new_file.append(line)


    with open(f'{file_name}_clean.pdb', 'w') as f:
        for i in new_file:
            f.write(i)


#PP: you can use bulit-in dictionaries for amino acids, e.g. in biopython
acids = {'A':'ALA', 'R':'ARG', 'D': 'ASP', 'N':'ASN', 'C':'CYS', 
           'E':'GLU', 'Q':'GLN', 'G':'GLY', 'H':'HIS', 'I': 'ILE',
            'L': 'LEU', 'K':'LYS', 'M':'MET', 'F':'PHE',
           'P': 'PRO', 'S':'SER', 'T':'THR', 'W':'TRP', 'Y': 'TYR',
           'V':'VAL'}

#PP: think in advance of usage on other devices, e.g. all the pathes should be either input parameters or relative ones
gmx = '/usr/local/bin/gmx'
rosetta = '/home/ekaterina_andreichuk/Downloads/rosetta_src_2021.16.61629_bundle/main/source/bin/rosetta_scripts.default.linuxgccrelease'

df = pd.read_csv('dataset_prot.csv', delimiter = ';')
for row in df.iterrows():
    pdb_id = row[1].pdb_id
    position = row[1].position
    wild_type = row[1].wild_type
    mutation = row[1].mutation
    if os.path.isdir(pdb_id):
    #PP: check that each 'if' has 'else'

        os.chdir(pdb_id)
        #PP: avoid jumping over directories. You can make path variable using os.path.join(workdir, pdb_id)

        print(os.getcwd())
        clean_pdb(pdb_id)
        dir_name = wild_type+str(position)+mutation
        if not os.path.isdir(dir_name):
            os.mkdir(dir_name)
        os.chdir(dir_name)
        try:
            print(dir_name)
            if not os.path.isfile(f'{dir_name}.pdb'):
                cur_dir = os.getcwd()
                #create mutate.xml
                shutil.copyfile(os.path.abspath('../../mutate.xml'), cur_dir+'/mutate.xml')
                template_string = f'''<MutateResidue name="one_point_mutation" target="{str(position)}A" new_res="{acids[mutation]}" preserve_atom_coords="false"
                mutate_self="false" update_polymer_bond_dependent="false"/>
                '''
                temp_string = '<Add mover="one_point_mutation"/>'
                new_text = []
                with open('mutate.xml') as f:
                    for line in f:
                        if '/MOVERS' in line:
                            new_text.append(template_string)
                        elif '</PROTOCOLS>' in line:
                            new_text.append(temp_string)
                        new_text.append(line)
                with open('mutate.xml', 'w') as f:
                    for i in new_text:
                        f.write(i)
                args = [rosetta, '-s', f'../{pdb_id}_clean.pdb', '-parser:protocol', 'mutate.xml'] 
                subprocess.call(args)
                args = ['mv', f'{pdb_id}_clean_0001.pdb', f'{dir_name}.pdb']
                subprocess.call(args)
    #               if True: #add check
    #                   print('mutagenesis sucessfull')

            if not os.path.isfile('emp.gro'):
                #PP: what is emp.gro?

                args = subprocess.Popen([ 'printf', '6\n1\n'], stdout=subprocess.PIPE)
                args_ = subprocess.Popen(['gmx', 'pdb2gmx', '-f', dir_name, '-o', 'pep.gro', '-p', 'pep.top', '-ignh', '-q'], 
                                         stdin = args.stdout, stdout= subprocess.DEVNULL)
                output = args_.communicate()

                #PP this subprocess will be successful, even if gmx returns error, so check that gmx worked properly, i.e. all required files were created.
                args = ['gmx', 'editconf', '-f', 'pep_.gro', '-o', 'pep.gro', '-d', '0.3', '-bt', 'cubic',  '-c']
                subprocess.call(args, stdout = subprocess.DEVNULL)

                #PP: what is sys.top?
                shutil.copyfile(os.path.abspath('../../sys.top'), cur_dir+'/sys.top') 

                # todo: change the name of the protein in the name of the system
                text= []
                mark = False
                with open('pep.top') as file:
                    for line in file:
                        if 'amber' not in line and not mark:
                            text.append(line)
                        if 'posre.itp' in line:
                            mark = True
                text.append('#endif')
                with open('pep.itp', 'w') as f:
                    for i in text:
                        f.write(i)
                shutil.copyfile(os.path.abspath('../../em.mdp'), cur_dir+'/em.mdp')
                shutil.copyfile(os.path.abspath('../../eq.mdp'), cur_dir+'/eq.mdp')
                shutil.copyfile(os.path.abspath('../../md.mdp'), cur_dir+'/md.mdp')
                args = ['gmx', 'grompp', '-f', 'em', '-p', 'sys', '-o', 'emp', '-c', 'pep.gro']
                subprocess.call(args, stdout = subprocess.DEVNULL)
                args = ['gmx', 'mdrun', '-deffnm', 'emp', '-v']
                subprocess.call(args, stdout = subprocess.DEVNULL)
                if not os.path.isfile('emp.tpr') or not os.path.isfile('emp.trr'):
                    print('error at emp stage')
                    continue
            if not os.path.isfile('em.gro'):
                args = ['gmx', 'solvate', '-cp', 'emp.gro', '-cs', '-o', 's0.gro', '-p', 'sys.top']
                subprocess.call(args, stdout = subprocess.DEVNULL)
                args = ['gmx', 'grompp', '-f', 'em', '-c', 's0', '-p', 'sys', '-o', 'ion']
                subprocess.call(args, stdout = subprocess.DEVNULL)
                args = subprocess.Popen([ 'echo', '13'], stdout=subprocess.PIPE)
                args_ = subprocess.Popen(['gmx', 'genion', '-s', 'ion', '-neutral', '-conc', '0.1', '-p', 'sys.top', '-o', 'start.gro'], 
                                         stdin = args.stdout, stdout= subprocess.DEVNULL)
                args_.communicate()
                args = ['gmx', 'grompp', '-f', 'em', '-c', 'start', '-p', 'sys', '-o', 'em']
                subprocess.call(args)
                args = ['gmx', 'mdrun', '-deffnm', 'em', '-v']
                subprocess.call(args)
                if os.path.isfile('em.gro') and os.path.isfile('em.trr'):
                    print('sucessful minimization')
                else: print('error at minimization stage')
            if not os.path.isfile('eq.gro'):
                args = ['gmx', 'grompp', '-f', 'eq', '-c', 'em', '-p', 'sys', '-o', 'eq']
                subprocess.call(args)
                args = ['gmx', 'mdrun', '-deffnm', 'eq', '-v']
                subprocess.call(args)
            if not os.path.isfile('md.gro'):
                args = ['gmx', 'grompp', '-f', 'md', '-c', 'eq', '-p', 'sys', '-o', 'md']
                subprocess.call(args)
                args= ['gmx', 'mdrun', '-deffnm', 'md']
                subprocess.call(args)

            #PP: avoid using break
            break
        finally:
            os.chdir('../../')


