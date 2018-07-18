
# Amino acids inventory
This program generates an inventory of amino acids present on the surface of peripheral proteins of the same family.
(See report.pdf and presentation.pdf for more details)

### Prerequisites
UNIX/shell environment<br />
Python3, R, PyMOL

### Running the program
```
$./launch.sh
```
Then choose a CATH superfamily (3.20.20.190 or 2.29.30.29)

### Authors

* **Hélène Kabbech** - Bioinformatics Master student (University of Paris Diderot)

* Internship supervised by Prof. **Nathalie Reuter** (CBU Bergen)

Spring 2018

softwares used : mkdssp, hmmbuild, hmmalign & seqret
PyMol

Superfamily

1. Files to provide
Structural alignment file (/fasta/cath/alignment.fasta)
pdb files (/pdb)
Sequences to add (/fasta/new/seq_to_align.fasta)
