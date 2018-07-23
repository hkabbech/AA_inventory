
# Amino acids inventory
This program generates an inventory of amino acids present on the surface of peripheral proteins of a same family.<br />
(See report.pdf and presentation.pdf for more details)

### Prerequisites
UNIX/shell environment<br />
Python3, R, PyMOL

### Running the program with the current cases
```
$./launch.sh
```
Then choose a CATH superfamily (3.20.20.190 or 2.29.30.29)

### Menu items
1. Display the structural alignment colored by secondary structure
2. Display the structural alignment in PyMOL
3. Study of a fold/loop (gives out a csv file of all amino acids) 
4. Add sequences to the structural alignment
5. Generate a plot (amino acid inventory)
6. Delete unused files and Quit

### Add a new case
1. Choose a **CATH superfamily** and create a directory _/superfamilies/NEW-FAMILY_
2. Provide a **structural alignment** file (FASTA format) _/NEW-FAMILY/fasta/cath/alignment.fasta_
3. Provide the **pdb files** corresponding _/NEW-FAMILY/pdb_
4. Provide a **file containing sequences to add** (FASTA format) : _/NEW-FAMILY/fasta/new/seq_to_align.fasta_
5. In the **header of the python script** (_/scripts/script.py_) indicate the superfamily studied, the model chosen and the different secondary structures to analyse (name, type and number in the sequence)

### Authors

* **Hélène Kabbech** - Bioinformatics Master student (University of Paris Diderot)
* Internship supervised by **Prof. Nathalie Reuter** (CBU Bergen)
