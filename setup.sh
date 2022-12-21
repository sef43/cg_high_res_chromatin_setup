N=3

python creation.py

python replicate_dna_seq.py $N > replicate_settings.txt

cp new_DNA_sequence.txt DNA_sequence.txt

cat in.nucl run_settings.txt > run.in

./lmp -in run.in



