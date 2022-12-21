# cg_high_res_chromatin_setup

Setup scripts for high-resolution chromatin model in https://doi.org/10.1038/s41467-021-23090-3.

1. first look at https://github.com/CollepardoLab/CollepardoLab_Chromatin_Model.

2. This repo contains the python script to make the initial nucleosome structure and the lammps script to replicate it into a fiber.

3. `setup.sh` contains the steps (`bash setup.sh`):
  - create a single nucleosome:
    ```
    python creation.py
    ```
    change the value of NRL in creation.py to set the nucleosome NRL.
    this creates `nucl.txt` - the lammps data file. And `in.nucl` - the first part of the lammps input script.
    
  - setup the replication settings:
    ```
    python replicate_dna_seq.py $N > replicate_settings.txt
    cp new_DNA_sequence.txt DNA_sequence.txt
    ```
    where `$N` is the number of nucleosomes in the fiber
    
  - make the input script with the run settings in `run_settings.txt`
    ```
    cat in.nucl run_settings.txt > run.in
    ```
    
  - run the simulation:
    ```
    ./lmp -in run.in
    ```
    
 The prodecure uses the lammps `fix drag` to pull the DNA linkers, it then uses `replicate` to create more nucleosomes, finally it uses `create_bonds` to join the nucleosomes together.
