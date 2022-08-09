Scripts for the paper:
# "Low Complexity Induces Structure in Protein Regions Predicted as Intrinsically Disordered"

DOI: pending

In this project we evaluated the intrinsically disordered regions (IDRs) and that contain simple repeated regions composed by one amino acid (homorepeat or polyX) or two amino acids (di-repeats or polyXYs), when aligned to the sequencies of PDB structures in search for secondary structure patterns and other physical-chemical characteristics.

### Python scripts
-----

Packages versions: python 3.8.10, biopython 1.79, gemmi 0.5.5, localcider 0.1.19, lxml 4.9.0, numpy 1.22.4, pandas, 1.4.2 scipy 1.8.1, urllib3 1.26.10

**Required inputs generated by external tools:**
* MobiDB json file version 4.1 of the target organism (for canonical unioprot organisms) or all set of protein of specific organisms (check help options to verify how to download a json file through API);
* Uniprot proteome of the target organism;
* PolyX2 tab output for the target organism runned with the same number of identical residues as the local window of amino acids to generate homorepeats (http://cbdm-01.zdv.uni-mainz.de/~munoz/polyx2/);
* PolyXY tab output for the target organism (pending publication).
* blastP xml output file. The blast must be executed against a customized databased generated with our scritp. **We mask the IDR regions, prerequisite to the following analysis.**


_**IMPORTANT:** The scripts were designed to run over all protein sequences and annotations of complete proteomes. Filtering of specific proteins must be done later on the .csv output files._

**main.py:** Main function where that automate the process of calling the other scripts. Use not required but recommended. The option 0 is still not fully implmented, but options 1 to 3 are fully operational.

**idr.py:** Contain all the transformation functions between extract the IDRs from MobiDB json file and generate its physical-chemical properties;

**pdbDssp.py:** Downloads all PDB structures and annotate secondary strutures based on the PDB or CIF files. Extracts additional usefull information from the PDB files. After the first run, this function incrementally downloads new structures and appends the output to the previous file available on disk, also deleting obsolete PDB structures.

**idrPdb.py:** Merge the data provided by the blastP XML file, fitting the secondary structures annotated with DSSP to the target IDR regions. It selects the more relevant structure according with the criterias discussed in the paper and provide several csv files for further analysis.

**poly.py:** Overlap the polyX and polyXYs regions with IDRs and PDB aligned sequences to allow further analysis.

**resources.py:** Global functions required by multiple of the scripts described above.
