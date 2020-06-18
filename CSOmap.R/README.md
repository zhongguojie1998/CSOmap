# CSOmap
======== Installation ========

This is a R package. 

We have tested this package on these platforms:

Mac Os: R 3.6.0.

You don't have to install anything, but you need to run the code under folder CSOmap.R/

======== Demo ========

We have provided a demo dataset, to run it:

(1). change your R current working directory to CSOmap_R/ (on Code Ocean, the name is code/)

(2). type a single line of code:

    runme('demo');

It will take about 5 minutes to run.

(3). Outputs will be saved under results/demo/, the structure of this folder is below:
    
    coordinates.txt
    counts.txt
    pvalue.txt
    qvalue.txt
    
(4). For details, continue reading.

======== Overview ========

This package is the software for our paper "Reconstruction of cell spatial organization from single-cell RNA sequencing data based on ligand-receptor mediated self-assembly"

To run this package, you should read and follow the instructions as below:

(1). Data Preparation:

All the data should be saved in a single folder under 'data/' (All the tables mentioned below should be tab-seperated): 

 For TPM data, the first column should be gene names (ATTENTION! NOT GENE ID!); The first row should be cell names, starting with a leading tab (\t). The contants of this table should be the TPM (ATTENTION! NOT log(TPM+1)!) value. The file name should be TPM.txt

 For label data, the first columns should be cell names (we use that to match the cell names in TPM.txt, order is not necessary) and the second column should be the corresponding label. The file name should be label.txt.

 For ligand receptor data, each row is a pair of ligand-receptor and the first, second, and third columns are the ligand, receptor, and their interacting weight, respectively. The file name should be LR_pairs.txt.

Please see the folder demo/ to view an example of a prepared dataset. You can use the LR_pairs.txt for human in that folder by just copying the LR_pairs.txt file to your own directory.
 
(2). Run Main Function:

After data preparation, you just need to type a single line of code in R: 

    CSOmap('folder_name'); 

======== Outputs ========

All the outputs will be saved in a new directory 'results/folder_name'

    coordinates.txt: coordinates
    counts.txt: connection number
    pvalue.txt: statistical test result
    qvalue.txt: adjusted pvalue

and it will return a list() object which stores coords, counts, pvalue and qvalue
    


