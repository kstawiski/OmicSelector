# DearWXpub
## wx keras version for everyone
## Web tool : https://wx.deargendev.me/
## We show the feature selection and Cancer/Normal classification result on TCGA dataset

A Keras implementation of Wx in preprint, :   
**[Wx: a neural network-based feature selection algorithm for next-generation sequencing data Sungsoo Park, Bonggun Shin, Yoonjung Choi, Kilsoo Kang, and Keunsoo Kang]**
(https://www.biorxiv.org/content/biorxiv/early/2017/11/18/221911.full.pdf)   

 
**Differences with the paper:**   
- We use the Keras as Neural Network running framework, original paper used naive tensorflow framework
- Selected Features can be different as which backend learning framework used
- Some TCGA samples updated, so we have more samples than paper.

**Contacts**
- Your contributions to the repo are always welcome. 
Open an issue or contact me with E-mail `sspark@deargen.me`


## Usage

**Step 1.**
```
Experienment Environments

$ python 3.4
$ tensorflow gpu 1.4.0
$ keras 2.1.2
```

```
Install instructions

$ pip install tensorflow-gpu
$ pip install keras
```

**Step 2. Clone this repository to local.**
```
$ git clone https://github.com/deargen/DearWXpub.git
$ cd DearWXpub
```


**Step 3. Download the TCGA RNA-seq data**

1. Download rna-seq via TCGA-Assembler tool

  - we downloaded tool (`Module_A.R`) via  http://www.compgenome.org/TCGA-Assembler/index.php
  
2. Run `tcga_download.R`

  - you may have to install releative packages
```
 - Shall Terminal
 $ sudo apt-get install libssl-dev
   * deb: libssl-dev (Debian, Ubuntu, etc)
   * rpm: openssl-devel (Fedora, CentOS, RHEL)
   * csw: libssl_dev (Solaris)
   * brew: openssl@1.1 (Mac OSX)
 
 $ sudo apt-get install libcurl4-openssl-dev
   * deb: libcurl4-openssl-dev (Debian, Ubuntu, etc)
   * rpm: libcurl-devel (Fedora, CentOS, RHEL)
   * csw: libcurl_dev (Solaris)

- R-studio
 > install.packages('curl')
 > install.packages('httr')
 > install.packages('stringr')
 > install.packages('rjson')
```
  
  
  - you can see the 'TCGA_DATAS' folder in current DearWXpub path.
 
3. Data Status

| 　   | TOTAL | Tumor | Normal | Tumor Ratio(%) |
|------|-------|-------|--------|----------------|
| TYPE | **6226**  | **5609**  | **617**    | **90.09**  |
| BLCA | 427   | 408   | 19     | 95.55          |
| BRCA | 1214  | 1101  | 113    | 90.69          |
| COAD | 327   | 286   | 41     | 87.46          |
| HNSC | 566   | 522   | 44     | 92.23          |
| KICH | 90    | 65    | 25     | 72.22          |
| KIRC | 606   | 534   | 72     | 88.12          |
| KIRP | 323   | 291   | 32     | 90.09          |
| LIHC | 424   | 374   | 50     | 88.21          |
| LUAD | 576   | 517   | 59     | 89.76          |
| LUSC | 553   | 502   | 51     | 90.78          |
| PRAD | 549   | 497   | 52     | 90.53          |
| THCA | 571   | 512   | 59     | 89.67          |

**Step 4. Do the Feature selection and Get Classification Accuracy**

```
$ python src/wx_tcga.py
```
It will generate preprocessed TCGA data set. And Select features, Get the scores.


## Results

#### Selected Gene Markers
***Keras Wx 14***
`['EEF1A1','FN1','GAPDH','SFTPC','AHNAK','KLK3','UMOD','CTSB','COL1A1','GPX3','GNAS','ATP1A1','SFTPB','ACTB']`

***Peng 14***
`['KIF4A','NUSAP1','HJURP','NEK2','FANCI','DTL','UHRF1','FEN1','IQGAP3','KIF20A','TRIM59','CENPL','C16ORF59','UBE2C']`
 
 with WX ranking
`['KIF4A'(12970),'NUSAP1'(15886),'HJURP'(11479),'NEK2'(11939),'FANCI'(17123),'DTL'(12600),'UHRF1'(11825),
'FEN1'(17497),'IQGAP3'(14173),'KIF20A'(13057),'TRIM59'(11113),'CENPL'(10344),'C16ORF59'(9463),'UBE2C'(9039)]`


***edgeR 14***
`['LCN1','UMOD','AQP2','PATE4','SLC12A1','OTOP2','ACTN3','KRT36','ATP2A1','PRH2','AGER','PYGM','PRR4','ESRRB']`

with WX ranking
`['LCN1'(3847),'UMOD'(7),'AQP2'(106),'PATE4'(5878),'SLC12A1'(165),'OTOP2'(7785),'ACTN3'(4939),'KRT36'(9441),
'ATP2A1'(14193),'PRH2'(2155),'AGER'(934),'PYGM'(13688),'PRR4'(2151),'ESRRB'(9630)]`

#### Cancer Classifiation Accuracy
TCGA data( Downloaded at Dec. 26th. 2017 ), Half of data for feature selection / Half of data for validation

**Accuracy of 14 Biomarker**

|       |         | Wx 14  |       | Peng 14 |       | EdgeR 14 |        |
|:-----:|---------|--------|-------|---------|-------|----------|--------|
| TYPE  | SAMPLES | Hit    | Acc(%)| Hit     | Acc(%)| Hit      | Acc(%) |
| **TOTAL** | **3119** | **3017** | **96.72** | **2961** | **94.93** | **2957** | **94.81** |
| BLCA  | 214     | 205    | 95.79 | 208     | 97.20 | 203      | 94.86  |
| BRCA  | 608     | 597    | 98.19 | 586     | 96.38 | 558      | 91.78  |
| COAD  | 164     | 155    | 94.51 | 143     | 87.20 | 162      | 98.78  |
| HNSC  | 283     | 275    | 97.17 | 261     | 92.23 | 267      | 94.35  |
| KICH  | 46      | 44     | 95.65 | 44      | 95.65 | 46       | 100.00 |
| KIRC  | 303     | 302    | 99.67 | 293     | 96.70 | 301      | 99.34  |
| KIRP  | 162     | 161    | 99.38 | 158     | 97.53 | 161      | 99.38  |
| LIHC  | 212     | 192    | 90.57 | 201     | 94.81 | 186      | 87.74  |
| LUAD  | 289     | 283    | 97.92 | 282     | 97.58 | 286      | 98.96  |
| LUSC  | 277     | 272    | 98.19 | 268     | 96.75 | 275      | 99.28  |
| PRAD  | 275     | 257    | 93.45 | 260     | 94.55 | 254      | 92.36  |
| THCA  | 286     | 274    | 95.80 | 257     | 89.86 | 258      | 90.21  |

**WX 14 AUC curve (LUSC, LUAD, BRCA)**

![roc_wx14_lusc](https://user-images.githubusercontent.com/4970569/36409868-60e728e8-1651-11e8-9fef-11c660fc3f68.png)
![roc_wx14_luad](https://user-images.githubusercontent.com/4970569/36409921-9d9019c6-1651-11e8-8520-453fcd4d4afa.png)
![roc_wx14_brca](https://user-images.githubusercontent.com/4970569/36410142-a0be1c3c-1652-11e8-9f2a-59fd97d263a3.png)



**Accuracy of 7 Biomarker**

***Keras wx 7***

`'EEF1A1', 'FN1', 'GAPDH', 'SFTPC', 'AHNAK', 'KLK3', 'UMOD'`

***Martinez 7***

 `BLCA : 'SMAD2', 'RUNX2', 'ABTB1', 'ST5', 'CEBPB', 'SETDB1', 'CEBPG'`
 
 `BRCA : 'JAK2', 'NFKBIA', 'TBP', 'RXRA', 'VAV1', 'HES5', 'NFKBIB'`
 
 `HNSC : 'DUSP16', 'KRT8', 'RAF1', 'MED1', 'PPARG', 'YWHAB', 'FABP1'`
 
 `KIRC : 'AR', 'HGS', 'RUNX1', 'BCL3', 'BRCA1', 'STAT2', 'ITGA8'`
 
 `LUAD : 'DOK1', 'FUT4', 'INSR', 'ITGB2', 'SHC1', 'PTPRC', 'KHDRBS1'`
 
 `LUSC : 'BRCA1', 'ETS2', 'HIF1A', 'JUN', 'LMO4', 'PIAS3', 'RBBP7'`
 

| 　    | 　   | WX 7 | 　      | Martinez 7 | 　     |
|-------|------|------|---------|------------|--------|
| TYPE  | SAMPLES | Hit    | Acc(%)| Hit     | Acc(%)|
| TOTAL | 3119 | 2986 | 95.74   | 　         | 　     |
| BLCA  | 214  | 205  | 95.79   | 206        | 96.26  |
| BRCA  | 608  | 591  | 97.20   | 556        | 91.45  |
| COAD  | 164  | 152  | 92.68   | 　         | 　     |
| HNSC  | 283  | 269  | 95.05   | 268        | 94.70  |
| KICH  | 46   | 45   | 97.83   | 　         | 　     |
| KIRC  | 303  | 299  | 98.68   | 273        | 90.10  |
| KIRP  | 162  | 162  | 100.00  | 　         | 　     |
| LIHC  | 212  | 187  | 88.21   | 　         | 　     |
| LUAD  | 289  | 283  | 97.92   | 260        | 89.97  |
| LUSC  | 277  | 271  | 97.83   | 257        | 92.78  |
| PRAD  | 275  | 249  | 90.55   | 　         | 　     |
| THCA  | 286  | 273  | 95.45   | 　         | 　     |

## Wx example from actual data file

```
$ python src/wx_examples.py
```
## Wx example for k > 2

```
$ python src/wx_gse_multi_class.py
```


## Additional Validation

***GSE72056***

```
Cancer type : melanoma

Sample size
malignant : 1257
normal : 3256
total : 4513

Dataset splits(trn,dev,tst) = (3611, 902)
We perfomred the 5 fold cross validation

Accuracy Result (test set)
WX14 : 90.71%
Peng14 : 70.22%
```

***GSE40419***

```
Cancer type : luad(lung)

Sample size
malignant : 87
normal : 77
total : 164

Dataset splits(trn,dev,tst) = (132,32)
We perfomred the 5 fold cross validation

Accuracy Result(test set)
WX14 : 80.00%
Peng14 : 56.87%
```

***GSE103322***

```
Cancer type : Head and Neck(sinle cell)

Sample size
malignant : 2215
normal : 3687
total : 5902

Dataset splits(trn+dev,tst) = (4722,1180)
We perfomred the 5 fold cross validation

Accuracy Result(test set, 1180)
WX14 : 81.10%
Peng14 : 68.28%
```

