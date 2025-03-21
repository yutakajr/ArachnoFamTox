# ArachnoFamTox: prediction and classification of Arachnids venom and toxin families 
### Fernanda M. Abukawa, Marcela A. Ishihara, Wesley F. Mantovani, Leo K. Iwai, Alexandre K. Tashima,  Milton Y. Nishiyama-Jr

# Description

:computer:**ArachnoFamTox** is a standalone tool and can be used to efficiently search for new toxins and venom proteins from protein sequences. ArachnoFamToxDB can be downloaded and used as a reliable database for Arachnids venom proteins. We anticipate ArachnoFamTox as a pivotal tool to enhance biodiscovery in the study of arachnid venoms.

# Basic Requirements and Dependencies

* [**Python3.8**](https://www.python.org/downloads/release/python-380/)
* [**BioPython**](https://biopython.org/)
* [**Cython**](https://cython.org/)
* [**Numpy**](https://numpy.org/)
* [**Pandas**](https://pandas.pydata.org/)

* [**HMMER 3**](https://hmmer.org)  
  Used for HMM profile prediction.   
  *Eddy SR, Accelerated Profile HMM Searches. PLOS Computational Biology 2011, 10.1371/journal.pcbi.1002195*

* [**BLAST 2.11**](https://blast.ncbi.nlm.nih.gov/Blast.cgi)  
  Used for RPS-BLAST and BLASTp prediction.    
  *Altschul, Stephen F., et al. "Gapped BLAST and PSI-BLAST: a new generation of protein database search programs." Nucleic acids research 25.17 (1997): 3389-3402.*

* [**SIGNALP 5.0b**](https://www.nature.com/articles/s41587-019-0036-z)
  Used for RPS-BLAST and BLASTp prediction.    
  *Almagro Armenteros, José Juan, et al. "SignalP 5.0 improves signal peptide predictions using deep neural networks." Nature biotechnology 37.4 (2019): 420-423.*

# Installation

Before installing the software, its recommended that you create a Python Virtual Environment (Python 3.8).
For example, using python's `venv`:

- Clone the repository on GitHub, and enter the program folder:
```
git clone https://github.com/yutakajr/ArachnoFamTox
cd ArachnoFamTox
```

- (Optional) Create a virtual environment for Python 3.8 and activate it:
```
python3.8 -m venv arachnofamtox_env
source arachnofamtox_env/bin/activate
```

- Update the pip3.8 to the latest version:
```
pip3.8 install --upgrade pip
```

- Install packages in new environment
```
pip3.8 install -r requirements.txt
```

- Install ArachnoFamTox
```
python setup.py install 
```

# Usage

To run the software on a peptide file:

```
ArachnoFamTox -fasta <peptides.pep> -out out_test_dir 
```

Use ```ArachnoFamTox -h``` to see aditional documentation and options.

```
  -h, --help                show this help message and exit
  -fasta <fasta file>       Specify fasta file
  -path <string>            Specify models path. Default=db
  -model_name_PSSM <string> Specify models name. Default=ArachnidaToxinsV3
  -model_name_HMM <string>  Specify models name. Default=Arachnida
  -eHMM <evalue>            e-value for HMMSCAN. Default=1e-1
  -ePSSM <evalue>           e-value for PSSM. Default=1e-5
  -eBLASTP <evalue>         e-value for BLASTp. Default=1e-5
  -qcovsfilter <float>      Filter BLASTp output for qcovs >= <float>. Default=Off
  -pposfilter <float>       Filter BLASTp output for ppos >= <float>. Default=Off
  -pidentfilter <float>     Filter BLASTp output for pident >= <float>. Default=Off
  -evaluefilter <float>     Filter BLASTp output for evalue <= <evalue>. Default=1e-10
  -cpu <int>                Specify number of threads. Default=1
  -out <output folder>      Specify directory to output
  --force                   Force re-use output directory. Default=Off.
  --tempfiles               Maintain temporary files. Default=Off.
  --signalp                 Perform Signalp analysis. Default=Off.
 ```

## Running ArachnoFamTox
### Models 

* The default PSSM and HMM models directory is "db" (default). To use ArachnoFamTox databases, make sure to run command with flag -path /path_of_installation/ArachnoFamTox/db (if running outside directory of installation).

* To test with other PSSM, HMM and BLASTp models, use flag -path to specify the directory with new models and flags -model_name_PSSM and -model_name_HMM to specify model names for PSSM and HMM respectively. 
```
ArachnoFamTox.py -fasta <peptides.pep> -out <output_dir> -path </new_models/db> -model_name_PSSM <PSSM_MODEL> -model_name_HMM <HMM_MODEL>
```

### Output


Inside output directory, results files are created:
  - ```classification_results```: file with proteins classified by PSSM and HMM merged results;
  - ```toxprot_results```: file with proteins identified only by BLASTp search against ToxProtDB; 
  - ```out.pssm```: Output from RPS-BLAST search against ArachnoFamTox PSSM database; 
  - ```out.hmmer.domtab.parsed```: Output from HMMScan search against ArachnoFamTox HMM database; 
  - ```merged_outputs.tsv (optional)```: Temporary file with merged outputs from PSSM and HMM searches;
  - ```out.blastp.toxprot (optional)```: Temporary file with toxprot search result. 
  - ```signalp_predictions.tsv (optional)```: File with Signalp predictions.
  - ```signalp_predictions_filtered.tsv (optional)```: File with Signalp predictions filtered to show only sequences with predicted signal peptide.


# Reference

If you use or discuss **ArachnoFamTox**, its guide, or any script available at this repository, please cite:

Abukawa et al. (2024) **Improving classification of Arachnid venom proteins and toxins based on evolutionary conservation with ArachnoFamTox.** DOI:[XXXX](XXXX)

# Licence

[GPL v3](https://github.com/yutakajr/ArachnoFamTox/LICENSE)