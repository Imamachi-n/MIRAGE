# MIRAGE #
**Comprehensive miRNA target prediction software**

###-prerequisite###
Python-3.x  
ViennaRNA  
PyYAML  

For Ubuntu,  
**ViennaRNA installation**
```
sudo apt-add-repository ppa:j-4/vienna-rna
sudo apt-get update
sudo apt-get install vienna-rna
```

**PyYAML installation**
```
sudo apt-get install libyaml-dev
wget http://pyyaml.org/download/pyyaml/PyYAML-3.11.tar.gz
tar zxvf PyYAML-3.11.tar.gz
cd PyYAML-3.11
sudo python3 setup.py install
```

**Download Conservation Score file**  
(1)Wigfix files (e.g.chr1.phyloP46way.wigFix.gz) are downloaded from the following HTML.  
```
http://hgdownload.soe.ucsc.edu/goldenPath/hg19/phastCons46way/vertebrate/
http://hgdownload.soe.ucsc.edu/goldenPath/hg19/phyloP46way/vertebrate/
```

(2)The following scripts are run to prepare conservation score for RefSeq and miRBase transcriptome.  
```
python3 mirage_prepare.py phastcons_prep
python3 mirage_prepare.py phastcons_sizedown
python3 mirage_prepare.py phastcons_score_list
```

```
python3 mirage_prepare.py phylop_score_prep
python3 mirage_prepare.py phylop_sizedown
python3 mirage_prepare.py phylop_score_list
```
