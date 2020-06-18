# CSOmap
======== Installation ========

This is a matlab package. 

We have tested this package on these platforms:

Linux: Matlab R2016b, R2017b and R2018a.
Windows: Matlab R2016b.
Mac Os: Matlab R2018a.

You don't have to install anything, but you need to run the code under folder CSOmap.matlab/

======== Demo ========

We have provided a demo dataset, to run it:

1. change your matlab current working directory to CSOmap.matlab

2. type a single line of code:

	runme('demo');

It will take about 5 minutes to run.

3. Outputs will be saved under output/demo/, the structure of this folder is below:

output/demo/:
.
..
analyst.mat
data.mat
information.txt
result
workspace.mat

output/demo/result:
.
..
demo_3d_global_origin.pdf
demo_3dplot.fig
demo_3dplot.gif
demo_3dplot.pdf
demo_3d_views.pdf
demo_cellcounts.txt
demo_conclusion.txt
demo_connection_number_normalized.pdf
demo_connection_number.pdf
demo_connections.txt
demo_coordinate.txt
demo_density.pdf
demo_mainLR.txt
demo_qvalue.pdf
demo_sections_density.pdf
demo_sections_normal.pdf
demo_section_z=0.pdf
demo_statistics.txt
demo_totalcells.txt
rotate

output/demo/result/rotate:
.
..
demo_3dplot_100.pdf
demo_3dplot_10.pdf
demo_3dplot_110.pdf
demo_3dplot_120.pdf
demo_3dplot_130.pdf
demo_3dplot_140.pdf
demo_3dplot_150.pdf
demo_3dplot_160.pdf
demo_3dplot_170.pdf
demo_3dplot_180.pdf
demo_3dplot_190.pdf
demo_3dplot_200.pdf
demo_3dplot_20.pdf
demo_3dplot_210.pdf
demo_3dplot_220.pdf
demo_3dplot_230.pdf
demo_3dplot_240.pdf
demo_3dplot_250.pdf
demo_3dplot_260.pdf
demo_3dplot_270.pdf
demo_3dplot_280.pdf
demo_3dplot_290.pdf
demo_3dplot_300.pdf
demo_3dplot_30.pdf
demo_3dplot_310.pdf
demo_3dplot_320.pdf
demo_3dplot_330.pdf
demo_3dplot_340.pdf
demo_3dplot_350.pdf
demo_3dplot_360.pdf
demo_3dplot_40.pdf
demo_3dplot_50.pdf
demo_3dplot_60.pdf
demo_3dplot_70.pdf
demo_3dplot_80.pdf
demo_3dplot_90.pdf

4. For details, continue reading.

======== Overview ========

This package is the software for our paper "Reconstruction of cell spatial organization from single-cell RNA sequencing data based on ligand-receptor mediated self-assembly"

To run this package, you should read and follow the instructions as below:

1. Data Preparation:
All the data should be saved in a single folder under 'data/' (All the tables mentioned below should be tab-seperated): 

  For TPM data, the first column should be gene names (ATTENTION! NOT GENE ID!); The first row should be cell names, starting with a leading tab (\t). The contants of this table should be the TPM (ATTENTION! NOT log(TPM+1)!) value. The file name should be TPM.txt

  For label data, the first columns should be cell names (we use that to match the cell names in TPM.txt, order is not necessary) and the second column should be the corresponding label. The file name should be label.txt.

  For ligand receptor data, each row is a pair of ligand-receptor and the first, second, and third columns are the ligand, receptor, and their interacting weight, respectively. The file name should be LR_pairs.txt.

Please see the folder demo/ to view an example of a prepared dataset. You can use the LR_pairs.txt for human in that folder by just copying the LR_pairs.txt file to your own directory.
 
2. Run Main Function:
After data preparation, you just need to type a single line of code in Matlab: 

    runme('folder_name'); 

where folder_name is the name of the folder where your data are.

To run on demo set, use code runme('demo');

This code will preprocess the data, reconstruct 3d, do statistical analysis and draw pictures with default configurations, if you are not satisfied with default configurations, you can refer to the annotations in runme.m for details.

Those pictures only contains basic information for this dataset, including 3d plot, density plot, etc. If you want to study a certain gene or compare this dataset with another, like we did in our paper, you can refer to functions in draw_pictures/ for details, or you can see draw_all_pictures.m for examples. We also used Microsoft Excel and Powerpoint to draw some of the pictures based on extracted information stored in analyst. 

======== Outputs ========

All the outputs will be saved in a new directory 'output/folder_name', including preprocessed data, matlab workspace in 3D reconstruction and an object analyst.mat.

The file information.txt will show you the basic features of this dataset, like how many ligand-receptor pairs are found in it, if this number is critically low, we suggest you check your data preparation.

data.mat stores the processed data

workspace.mat stores all the variables in 3D reconstruction

analyst.mat is an object, it stores all the information needed for analysis, including TPM, 3d coordinate, connections, etc. It also has some in-built functions which can automatically output pictures, and those pictures should be stored in folder 'output/folder_name/result'. For details about those in-built functions, see analyst.m for details.

======== About pictures ========

The object analyst stores all the information needed for analysis and provides several built-in functions which can automatically output pictures. Each function has an annotation and you can refer to file analyst.m for details.

If you want to redraw those pictures in our paper, just run function draw_all_picuters.m , it will use the functions in draw_pictures/, for details of those functions, see the corresponding annotation.

======== Other functions ========

We described how to run CSOmap on an ordinary dataset above. In the paper, we talked about in-silico experiments and some analysis on certain genes, however, it is impossible for us to provide a single line of code to automatically do that, since you should manually decide which gene to study and which kind of in-silico experiment you want to do. 

We suggest you read the annotations in the following functions:

for in-silico gene changes, see annotations in changegenes.m

for in-silico knock out cells, see annotations in knockout cells.m

for in-silico adoptive transfer, see annotations in pseudo.m

======== About reproducibility ========

We stored our results in the folder data/analysts/ . It contains all the analyst objects, we seperate results according to cancer types, see annotation in draw_all_pictures.m for details. To reproduce those pictures in our paper, just run draw_all_pictures.m , pictures will be saved in folder pictures_in_paper/ . For result of adoptive transfer, however, it's impossible to store all the analysts due to limited storage, we stored only information for drawing pictures in workspace_HNC.mat and workspace_melanoma.mat.

We use random seed based on your system time in myoptimize.m . If you want to use a certain random seed, you can change that yourself. It is unnecessary because this algorithm is robust enough: although the coordinates might change a little with different random seeds, the statistical result and spatial pattern will not change.
