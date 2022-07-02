# UKiS 

This program calculates the length distribution of the first introns and all introns based on the exon information of the transcript obtained from Refseq database. After measurement, the results are output as boxplot.

The Refseq transcript information is stored in the Source_Data directory. Alternatively, new data downloaded by yourself from refseq will also be supported.  

In addition to drawing a boxplot, add a dashed line indicating the intron length of a specific gene, in this case <i>TP53</i>.  

## Result Image
<img src="/Sample_Images/IntronSize.png" alt="Sample Image" title="Sample Image" width="200">

## Requirement

It requires an environment capable of running python3 and its two modules, matplotlib and pandas.  

### Source data
* Refseq.csv

### Python packages
#### Python 3.8.9  
* matplotlib (Version: 3.5.1)
* pandas (Version: 1.4.1)

## Installation

Install the necessary modules with the following commands.  

```bash
pip3 install matplotlib
pip3 install pandas
```

## Usage


```bash
git clone https://github.com/Taichi-Akase/UKiS/
cd UKiS
python3 IntronSizeDistribution.py
```

## Notes


## Authors
* Taichi Akase
  * E-mail: akase.t.aa@m.titech.ac.jp
* Kono Syunya
* Tomoyuki Ohno
* Yasunori Aizawa


## Reference  
T. Ohno <i>et al. bioRxiv</i> (2022)  
Biallelic and gene-wide genomic substitution for endogenous intron and retroelement mutagenesis in human cells  
doi: https://doi.org/10.1101/2022.03.09.482138  

## License  

