% this script will draw all the pictures in our paper, 
% each function's annotation could be found in that function's .m file
addpath('draw_pictures/');
addpath('support_analysis_functions/');
addpath('draw_pictures/gramm-master/');
mkdir('pictures_in_paper')
%% load data
try
load('data/analysts/lung.mat');
load('data/analysts/liver.mat');
load('data/analysts/CRC.mat');
load('data/analysts/HNC.mat');
load('data/analysts/melanoma.mat');
catch
disp('Please check you have folder data/ under CSOmap/');
end
%% CRC_NPT
mkdir('pictures_in_paper/CRC_NPT/');
draw_result3d_or_split_or_gif_with_gramm(CRC_NPT, 'pictures_in_paper/CRC_NPT/CRC_NPT_3d_global', 'normal');
draw_result3d_with_gramm(CRC_NPT, 'pictures_in_paper/CRC_NPT/CRC_NPT_3d_views', 0);
draw_sections_with_gramm(CRC_NPT, 'pictures_in_paper/CRC_NPT/CRC_NPT_sections_normal', 'normal');
draw_sections_with_gramm(CRC_NPT, 'pictures_in_paper/CRC_NPT/CRC_NPT_sections_density', 'density');
draw_bar_of_connection_number_with_gramm(CRC_NPT, 'pictures_in_paper/CRC_NPT/CRC_NPT_connection_number', 0);
draw_bar_of_connection_number_with_gramm(CRC_NPT, 'pictures_in_paper/CRC_NPT/CRC_NPT_connection_number_normalized', 1);
draw_density_with_gramm(CRC_NPT, 'pictures_in_paper/CRC_NPT/CRC_NPT_density',3,0.5);
draw_qvalue_with_gramm(CRC_NPT, 'pictures_in_paper/CRC_NPT/CRC_NPT_qvalue', [], 30, 5, 20, 0.6, 0.6, [0.1,0.1,0.72,0.72]);
%% liver_NPT
mkdir('pictures_in_paper/liver_NPT/');
draw_result3d_or_split_or_gif_with_gramm(liver_NPT, 'pictures_in_paper/liver_NPT/liver_NPT_3d_global', 'normal');
draw_result3d_with_gramm(liver_NPT, 'pictures_in_paper/liver_NPT/liver_NPT_3d_views', 0);
draw_sections_with_gramm(liver_NPT, 'pictures_in_paper/liver_NPT/liver_NPT_sections_normal', 'normal');
draw_sections_with_gramm(liver_NPT, 'pictures_in_paper/liver_NPT/liver_NPT_sections_density', 'density');
draw_bar_of_connection_number_with_gramm(liver_NPT, 'pictures_in_paper/liver_NPT/liver_NPT_connection_number', 0);
draw_bar_of_connection_number_with_gramm(liver_NPT, 'pictures_in_paper/liver_NPT/liver_NPT_connection_number_normalized', 1);
draw_density_with_gramm(liver_NPT, 'pictures_in_paper/liver_NPT/liver_NPT_density',3,0.5);
draw_qvalue_with_gramm(liver_NPT, 'pictures_in_paper/liver_NPT/liver_NPT_qvalue', [], 30, 5, 20, 0.6, 0.6, [0.1,0.1,0.72,0.72]);
%% lung_NPT
mkdir('pictures_in_paper/lung_NPT/');
draw_result3d_or_split_or_gif_with_gramm(lung_NPT, 'pictures_in_paper/lung_NPT/lung_NPT_3d_global', 'normal');
draw_result3d_with_gramm(lung_NPT, 'pictures_in_paper/lung_NPT/lung_NPT_3d_views', 0);
draw_sections_with_gramm(lung_NPT, 'pictures_in_paper/lung_NPT/lung_NPT_sections_normal', 'normal');
draw_sections_with_gramm(lung_NPT, 'pictures_in_paper/lung_NPT/lung_NPT_sections_density', 'density');
draw_bar_of_connection_number_with_gramm(lung_NPT, 'pictures_in_paper/lung_NPT/lung_NPT_connection_number', 0);
draw_bar_of_connection_number_with_gramm(lung_NPT, 'pictures_in_paper/lung_NPT/lung_NPT_connection_number_normalized', 1);
draw_density_with_gramm(lung_NPT, 'pictures_in_paper/lung_NPT/lung_NPT_density',2,0.5);
draw_qvalue_with_gramm(lung_NPT, 'pictures_in_paper/lung_NPT/lung_NPT_qvalue', [], 30, 5, 20, 0.6, 0.6, [0.1,0.1,0.72,0.72]);
%% CRC
mkdir('pictures_in_paper/CRC/');
draw_result3d_or_split_or_gif_with_gramm(CRC, 'pictures_in_paper/CRC/CRC_3d_global', 'normal');
draw_result3d_with_gramm(CRC, 'pictures_in_paper/CRC/CRC_3d_views', 0);
draw_sections_with_gramm(CRC, 'pictures_in_paper/CRC/CRC_sections_normal', 'normal');
draw_sections_with_gramm(CRC, 'pictures_in_paper/CRC/CRC_sections_density', 'density');
draw_bar_of_connection_number_with_gramm(CRC, 'pictures_in_paper/CRC/CRC_connection_number', 0);
draw_bar_of_connection_number_with_gramm(CRC, 'pictures_in_paper/CRC/CRC_connection_number_normalized', 1);
draw_qvalue_with_gramm(CRC, 'pictures_in_paper/CRC/CRC_qvalue', CRC.standards, 1.4, 0.14, 15, 0.25, 0.25);
draw_one_gene_with_gramm(CRC, 'CCR8', 'pictures_in_paper/CRC/CRC_CCR8_TPM', 10, 0.5);
draw_density_with_gramm(CRC, 'pictures_in_paper/CRC/CRC_density',10,2);
draw_result3d_with_section_with_gramm(CRC, 'pictures_in_paper/CRC/CRC_section_z=0', [], [0,90], [-inf,inf],[-inf,inf],[-5,5],'normal');
draw_one_gene_with_gramm(CRC, 'CCR8', 'pictures_in_paper/CRC/CRC_CCR8_TPM', 10, 2);
draw_one_gene_with_gramm(CRC, 'CCL4', 'pictures_in_paper/CRC/CRC_CCL4_TPM', 10, 2);
genes = {'CCR8', 'CCL4', 'MKI67', 'IFNG', 'TNF', 'PDCD1', 'CD274', 'CTLA4', 'CD80', 'CD86'};
for i = 1 : size(genes, 2)
    draw_result3d_with_genes_with_gramm(CRC, ['pictures_in_paper/CRC/CRC_', genes{i}], [30,30], [-inf,inf],[-inf,inf],[-inf,inf],genes(i));
    draw_result3d_with_genes_with_gramm(CRC, ['pictures_in_paper/CRC/CRC_', genes{i}, '_section'], [0,90], [-inf,inf],[-inf,inf],[-5,5],genes(i));
end
%% CRC_knockout_Treg: knockout CD4-CTLA4 cells in CRC
mkdir('pictures_in_paper/CRC_knockout_Treg/');
draw_result3d_or_split_or_gif_with_gramm(CRC_knockout_Treg, 'pictures_in_paper/CRC_knockout_Treg/CRC_knockout_Treg_3d_global', 'normal');
draw_result3d_with_gramm(CRC_knockout_Treg, 'pictures_in_paper/CRC_knockout_Treg/CRC_knockout_Treg_3d_views', 0);
draw_sections_with_gramm(CRC_knockout_Treg, 'pictures_in_paper/CRC_knockout_Treg/CRC_knockout_Treg_sections_normal', 'normal');
draw_sections_with_gramm(CRC_knockout_Treg, 'pictures_in_paper/CRC_knockout_Treg/CRC_knockout_Treg_sections_density', 'density');
draw_bar_of_connection_number_with_gramm(CRC_knockout_Treg, 'pictures_in_paper/CRC_knockout_Treg/CRC_knockout_Treg_connection_number', 0);
draw_bar_of_connection_number_with_gramm(CRC_knockout_Treg, 'pictures_in_paper/CRC_knockout_Treg/CRC_knockout_Treg_connection_number_normalized', 1);
draw_density_with_gramm(CRC_knockout_Treg, 'pictures_in_paper/CRC_knockout_Treg/CRC_knockout_Treg_density',3,0.5);
draw_qvalue_with_gramm(CRC_knockout_Treg, 'pictures_in_paper/CRC_knockout_Treg/CRC_knockout_Treg_qvalue', CRC_knockout_Treg.standards, 1.7, 0.17, 15, 0.25, 0.25);
%% liver
mkdir('pictures_in_paper/liver/');
draw_result3d_or_split_or_gif_with_gramm(liver, 'pictures_in_paper/liver/liver_3d_global', 'normal');
draw_result3d_with_gramm(liver, 'pictures_in_paper/liver/liver_3d_views', 0);
draw_sections_with_gramm(liver, 'pictures_in_paper/liver/liver_sections_normal', 'normal');
draw_sections_with_gramm(liver, 'pictures_in_paper/liver/liver_sections_density', 'density');
draw_bar_of_connection_number_with_gramm(liver, 'pictures_in_paper/liver/liver_connection_number', 0);
draw_bar_of_connection_number_with_gramm(liver, 'pictures_in_paper/liver/liver_connection_number_normalized', 1);
draw_qvalue_with_gramm(liver, 'pictures_in_paper/liver/liver_qvalue', liver.standards, 4, 0.4, 15, 0.2, 0.2);
draw_density_with_gramm(liver, 'pictures_in_paper/liver/liver_density',10,2);
draw_result3d_with_section_with_gramm(liver, 'pictures_in_paper/liver/liver_section_z=0', [], [0,90], [-inf,inf],[-inf,inf],[-5,5],'normal');
draw_one_gene_with_gramm(liver, 'CCR8', 'pictures_in_paper/liver/liver_CCR8_TPM', 10, 0.5);
genes = {'CCR8', 'CCL4', 'MKI67', 'IFNG', 'TNF', 'PDCD1', 'CD274', 'CTLA4', 'CD80', 'CD86'};
for i = 1 : size(genes, 2)
    draw_result3d_with_genes_with_gramm(liver, ['pictures_in_paper/liver/liver_', genes{i}], [30,30], [-inf,inf],[-inf,inf],[-inf,inf],genes(i));
    draw_result3d_with_genes_with_gramm(liver, ['pictures_in_paper/liver/liver_', genes{i}, '_section'], [0,90], [-inf,inf],[-inf,inf],[-5,5],genes(i));
end
%% lung
mkdir('pictures_in_paper/lung/');
draw_result3d_or_split_or_gif_with_gramm(lung, 'pictures_in_paper/lung/lung_3d_global', 'normal');
draw_result3d_with_gramm(lung, 'pictures_in_paper/lung/lung_3d_views', 0);
draw_sections_with_gramm(lung, 'pictures_in_paper/lung/lung_sections_normal', 'normal');
draw_sections_with_gramm(lung, 'pictures_in_paper/lung/lung_sections_density', 'density');
draw_bar_of_connection_number_with_gramm(lung, 'pictures_in_paper/lung/lung_connection_number', 0);
draw_bar_of_connection_number_with_gramm(lung, 'pictures_in_paper/lung/lung_connection_number_normalized', 1);
draw_qvalue_with_gramm(lung, 'pictures_in_paper/lung/lung_qvalue', lung.standards, 1.5, 0.15, 15, 0.2, 0.2);
draw_one_gene_with_gramm(lung, 'CCR8', 'pictures_in_paper/lung/lung_CCR8_TPM', 10, 1.5);
draw_density_with_gramm(lung, 'pictures_in_paper/lung/lung_density',10,2);
draw_result3d_with_section_with_gramm(lung, 'pictures_in_paper/lung/lung_section_z=0', [], [0,90], [-inf,inf],[-inf,inf],[-5,5],'normal');
genes = {'CCR8', 'CCL4', 'MKI67', 'IFNG', 'TNF', 'PDCD1', 'CD274', 'CTLA4', 'CD80', 'CD86'};
for i = 1 : size(genes, 2)
    draw_result3d_with_genes_with_gramm(lung, ['pictures_in_paper/lung/lung_', genes{i}], [30,30], [-inf,inf],[-inf,inf],[-inf,inf],genes(i));
    draw_result3d_with_genes_with_gramm(lung, ['pictures_in_paper/lung/lung_', genes{i}, '_section'], [0,90], [-inf,inf],[-inf,inf],[-5,5],genes(i));
end
%% 3D heatmaps of CCR8, CCL4, MKI67, IFNG, TNF, PDCD1, CD274, CTLA4, CD80, CD86
genes = {'CCR8', 'CCL4', 'MKI67', 'IFNG', 'TNF', 'PDCD1', 'CD274', 'CTLA4', 'CD80', 'CD86'};
for i = 1 : size(genes, 2)
    draw_result3d_with_genes_with_gramm(CRC, ['pictures_in_paper/CRC/CRC_', genes{i}], [30,30], [-inf,inf],[-inf,inf],[-inf,inf],genes(i));
    draw_result3d_with_genes_with_gramm(liver, ['pictures_in_paper/liver/liver_', genes{i}], [30,30], [-inf,inf],[-inf,inf],[-inf,inf],genes(i));
    draw_result3d_with_genes_with_gramm(lung, ['pictures_in_paper/lung/lung_', genes{i}], [30,30], [-inf,inf],[-inf,inf],[-inf,inf],genes(i));
end
%% sections of CCR8, CCL4, MKI67, IFNG, TNF, PDCD1, CD274, CTLA4, CD80, CD86
genes = {'CCR8', 'CCL4', 'MKI67', 'IFNG', 'TNF', 'PDCD1', 'CD274', 'CTLA4', 'CD80', 'CD86'};
for i = 1 : size(genes, 2)
    draw_result3d_with_genes_with_gramm(CRC, ['pictures_in_paper/CRC/CRC_', genes{i}, '_section'], [0,90], [-inf,inf],[-inf,inf],[-5,5],genes(i));
    draw_result3d_with_genes_with_gramm(liver, ['pictures_in_paper/liver/liver_', genes{i}, '_section'], [0,90], [-inf,inf],[-inf,inf],[-5,5],genes(i));
    draw_result3d_with_genes_with_gramm(lung, ['pictures_in_paper/lung/lung_', genes{i}, '_section'], [0,90], [-inf,inf],[-inf,inf],[-5,5],genes(i));
end

%% HNC
mkdir('pictures_in_paper/HNC/');
draw_result3d_or_split_or_gif_with_gramm(HNC, 'pictures_in_paper/HNC/HNC_3d_global', 'normal');
draw_result3d_with_gramm(HNC, 'pictures_in_paper/HNC/HNC_3d_views', 0);
draw_sections_with_gramm(HNC, 'pictures_in_paper/HNC/HNC_sections_normal', 'normal');
draw_sections_with_gramm(HNC, 'pictures_in_paper/HNC/HNC_sections_density', 'density');
draw_bar_of_connection_number_with_gramm(HNC, 'pictures_in_paper/HNC/HNC_connection_number', 0);
draw_bar_of_connection_number_with_gramm(HNC, 'pictures_in_paper/HNC/HNC_connection_number_normalized', 1);
draw_density_with_gramm(HNC, 'pictures_in_paper/HNC/HNC_density',3,0.5);
draw_qvalue_with_gramm(HNC, 'pictures_in_paper/HNC/HNC_qvalue', HNC.standards, 10, 1, 15, 0.2, 0.2);
draw_result3d_with_section_with_gramm(HNC, 'pictures_in_paper/HNC/HNC_section_z=0_normal', [], [0,90], [-inf,inf],[-inf,inf],[-5,5],'normal');
draw_result3d_with_section_with_gramm(HNC, 'pictures_in_paper/HNC/HNC_section_z=0_density', [], [0,90], [-inf,inf],[-inf,inf],[-5,5],'density');
genes={'CD63','TIMP1','PDPN','LAMC2','LAMB3'};
for i = 1 : size(genes, 2)
draw_result3d_with_genes_with_gramm(HNC, ['pictures_in_paper/HNC/HNC_', genes{i}], [-35,30], [-inf,inf],[-inf,inf],[-inf,inf],genes(i));
draw_result3d_with_genes_with_gramm(HNC, ['pictures_in_paper/HNC/HNC_', genes{i}, '_section'], [0,90], [-inf,inf],[-inf,inf],[-5,5],genes(i));
draw_compare_one_gene_with_gramm(HNC, melanoma, 'HNC',  'malignant', 'melanoma', 'malignant', genes{i}, ['pictures_in_paper/HNC/HNC_melanoma_compare_', genes{i}], 1, 0.2);
end
%% HNC_pEMT: identified pEMT cells, coordinates are the same as HNC
mkdir('pictures_in_paper/HNC_pEMT/');
draw_result3d_or_split_or_gif_with_gramm(HNC_pEMT, 'pictures_in_paper/HNC_pEMT/HNC_pEMT_3d_global', 'normal');
draw_result3d_with_gramm(HNC_pEMT, 'pictures_in_paper/HNC_pEMT/HNC_pEMT_3d_views', 0);
draw_result3d_with_section_with_gramm(HNC_pEMT, 'pictures_in_paper/HNC_pEMT/HNC_pEMT_section_z=0_normal', [], [0,90], [-inf,inf],[-inf,inf],[-5,5], 'normal');
draw_result3d_with_section_with_gramm(HNC_pEMT, 'pictures_in_paper/HNC_pEMT/HNC_pEMT_section_z=0_desity', [], [0,90], [-inf,inf],[-inf,inf],[-5,5], 'density');
draw_sections_with_gramm(HNC_pEMT, 'pictures_in_paper/HNC_pEMT/HNC_pEMT_sections_normal', 'normal');
draw_sections_with_gramm(HNC_pEMT, 'pictures_in_paper/HNC_pEMT/HNC_pEMT_sections_density', 'density');
draw_bar_of_connection_number_with_gramm(HNC_pEMT, 'pictures_in_paper/HNC_pEMT/HNC_pEMT_connection_number', 0);
draw_bar_of_connection_number_with_gramm(HNC_pEMT, 'pictures_in_paper/HNC_pEMT/HNC_pEMT_connection_number_normalized', 1);
draw_density_with_gramm(HNC_pEMT, 'pictures_in_paper/HNC_pEMT/HNC_pEMT_density',3,0.5);
draw_qvalue_with_gramm(HNC_pEMT, 'pictures_in_paper/HNC_pEMT/HNC_pEMT_qvalue', HNC_pEMT.standards, 8, 0.8, 20, 0.2, 0.2);
draw_qvalue_with_gramm(HNC_pEMT_origin, 'pictures_in_paper/HNC_pEMT/HNC_pEMT_origin_qvalue', HNC_pEMT_origin.standards, 8, 0.8, 20, 0.2, 0.2);
%% HNC_pEMT_new: identified p-EMT cells, combine T cells, B cells as 'others'
mkdir('pictures_in_paper/HNC_pEMT_new/');
draw_result3d_or_split_or_gif_with_gramm(HNC_pEMT_new, 'pictures_in_paper/HNC_pEMT_new/HNC_pEMT_new_3d_global', 'normal');
draw_result3d_with_gramm(HNC_pEMT_new, 'pictures_in_paper/HNC_pEMT_new/HNC_pEMT_new_3d_views', 0);
draw_result3d_with_section_with_gramm(HNC_pEMT_new, 'pictures_in_paper/HNC_pEMT_new/HNC_pEMT_new_section_z=0_normal', [], [0,90], [-inf,inf],[-inf,inf],[-2,2], 'normal');
draw_result3d_with_section_with_gramm(HNC_pEMT_new, 'pictures_in_paper/HNC_pEMT_new/HNC_pEMT_new_section_z=0_desity', [], [0,90], [-inf,inf],[-inf,inf],[-2,2], 'density');
draw_sections_with_gramm(HNC_pEMT_new, 'pictures_in_paper/HNC_pEMT_new/HNC_pEMT_new_sections_normal', 'normal');
draw_sections_with_gramm(HNC_pEMT_new, 'pictures_in_paper/HNC_pEMT_new/HNC_pEMT_new_sections_density', 'density');
draw_bar_of_connection_number_with_gramm(HNC_pEMT_new, 'pictures_in_paper/HNC_pEMT_new/HNC_pEMT_new_connection_number', 0);
draw_bar_of_connection_number_with_gramm(HNC_pEMT_new, 'pictures_in_paper/HNC_pEMT_new/HNC_pEMT_new_connection_number_normalized', 1);
draw_density_with_gramm(HNC_pEMT_new, 'pictures_in_paper/HNC_pEMT_new/HNC_pEMT_new_density',3,0.5);
draw_qvalue_with_gramm(HNC_pEMT_new, 'pictures_in_paper/HNC_pEMT_new/HNC_pEMT_new_qvalue', HNC_pEMT_new.standards, 10, 1, 15, 0.6, 0.6);
draw_compare_of_connection_number_with_gramm(HNC_pEMT_new, [], {'malignant', 'Fibroblast'}, {'p-EMT cells', 'Fibroblast'}, 'malignant', 'p-EMT cells', 'pictures_in_paper/HNC_pEMT_new/HNC_pEMT_malignant_p-EMT_compare', 'density');
%% HNC_pEMT_origin_new: use affinitymat to do prediction, identified p-EMT cells, combine T cells, B cells as 'others'
mkdir('pictures_in_paper/HNC_pEMT_origin_new/');
draw_bar_of_connection_number_with_gramm(HNC_pEMT_origin_new, 'pictures_in_paper/HNC_pEMT_origin_new/HNC_pEMT_origin_new_connection_number', 0);
draw_bar_of_connection_number_with_gramm(HNC_pEMT_origin_new, 'pictures_in_paper/HNC_pEMT_origin_new/HNC_pEMT_origin_new_connection_number_normalized', 1);
draw_density_with_gramm(HNC_pEMT_origin_new, 'pictures_in_paper/HNC_pEMT_origin_new/HNC_pEMT_origin_new_density',3,0.5);
draw_qvalue_with_gramm(HNC_pEMT_origin_new, 'pictures_in_paper/HNC_pEMT_origin_new/HNC_pEMT_origin_new_qvalue', HNC_pEMT_origin_new.standards, 10, 1, 15, 0.6, 0.6);
draw_compare_density_with_gramm(HNC_pEMT_origin_new, HNC_pEMT_new, {'p-EMT cells','Fibroblast', 'malignant'}, 'density(not embeded)', 'density(embeded)', 'pictures_in_paper/HNC_pEMT_origin_new/HNC_pEMT_before&after_embed_density_compare',1,0.1, 'pair');
draw_compare_of_connection_number_with_gramm(HNC_pEMT_origin_new, HNC_pEMT_new, {'p-EMT cells', 'p-EMT cells'}, {'p-EMT cells', 'p-EMT cells'}, 'not embeded', 'embeded', 'pictures_in_paper/HNC_pEMT_origin_new/HNC_pEMT_before&after_mapping_p-EMT_compare', 'density');
draw_compare_of_connection_number_with_gramm(HNC_pEMT_origin_new, HNC_pEMT_new, {'p-EMT cells', 'Fibroblast'}, {'p-EMT cells', 'Fibroblast'}, 'not embeded', 'embeded', 'pictures_in_paper/HNC_pEMT_origin_new/HNC_pEMT_before&after_mapping_p-EMT_Fibroblast_compare', 'density');
%% HNC_pEMT_CD63: up the expression of CD63 in malignant and p-EMT cells
mkdir('pictures_in_paper/HNC_pEMT_CD63/');
draw_result3d_or_split_or_gif_with_gramm(HNC_pEMT_CD63, 'pictures_in_paper/HNC_pEMT_CD63/HNC_pEMT_CD63_3d_global', 'normal');
draw_result3d_with_gramm(HNC_pEMT_CD63, 'pictures_in_paper/HNC_pEMT_CD63/HNC_pEMT_CD63_3d_views', 0);
draw_sections_with_gramm(HNC_pEMT_CD63, 'pictures_in_paper/HNC_pEMT_CD63/HNC_pEMT_CD63_sections_normal', 'normal');
draw_sections_with_gramm(HNC_pEMT_CD63, 'pictures_in_paper/HNC_pEMT_CD63/HNC_pEMT_CD63_sections_density', 'density');
draw_bar_of_connection_number_with_gramm(HNC_pEMT_CD63, 'pictures_in_paper/HNC_pEMT_CD63/HNC_pEMT_CD63_connection_number', 0);
draw_bar_of_connection_number_with_gramm(HNC_pEMT_CD63, 'pictures_in_paper/HNC_pEMT_CD63/HNC_pEMT_CD63_connection_number_normalized', 1);
draw_density_with_gramm(HNC_pEMT_CD63, 'pictures_in_paper/HNC_pEMT_CD63/HNC_pEMT_CD63_density',3,0.5);
draw_qvalue_with_gramm(HNC_pEMT_CD63, 'pictures_in_paper/HNC_pEMT_CD63/HNC_pEMT_CD63_qvalue', HNC_pEMT_CD63.standards, 10, 1, 15, 0.6, 0.6);
draw_compare_density_with_gramm(HNC_pEMT_new, HNC_pEMT_CD63, [], 'origin', 'CD63_over_expression', 'pictures_in_paper/HNC_pEMT_CD63/HNC_pEMT_CD63_density_compare',1,0.1);
%% melanoma
mkdir('pictures_in_paper/melanoma/');
draw_result3d_or_split_or_gif_with_gramm(melanoma, 'pictures_in_paper/melanoma/melanoma_3d_global', 'normal');
draw_result3d_with_gramm(melanoma, 'pictures_in_paper/melanoma/melanoma_3d_views', 0);
draw_sections_with_gramm(melanoma, 'pictures_in_paper/melanoma/melanoma_sections_normal', 'normal');
draw_sections_with_gramm(melanoma, 'pictures_in_paper/melanoma/melanoma_sections_density', 'density');
draw_bar_of_connection_number_with_gramm(melanoma, 'pictures_in_paper/melanoma/melanoma_connection_number', 0);
draw_bar_of_connection_number_with_gramm(melanoma, 'pictures_in_paper/melanoma/melanoma_connection_number_normalized', 1);
draw_density_with_gramm(melanoma, 'pictures_in_paper/melanoma/melanoma_density',8,1.5);
draw_qvalue_with_gramm(melanoma, 'pictures_in_paper/melanoma/melanoma_qvalue', melanoma.standards, 10, 1, 20, 0.2, 0.2);
draw_result3d_with_section_with_gramm(melanoma, 'pictures_in_paper/melanoma/melanoma_section_z=0', [], [0,90], [-inf,inf],[-inf,inf],[-2,2], 'normal');
draw_result3d_with_section_with_gramm(melanoma, 'pictures_in_paper/melanoma/melanoma_section_z=0', [], [0,90], [-inf,inf],[-inf,inf],[-2,2], 'density');
%%
genes={'CD63','TIMP1','PDPN','LAMC2','LAMB3'};
for i = 1 : size(genes, 2)
draw_result3d_with_genes_with_gramm(melanoma, ['pictures_in_paper/melanoma/melanoma_', genes{i}], [-35,30], [-inf,inf],[-inf,inf],[-inf,inf],genes(i));
draw_result3d_with_genes_with_gramm(melanoma, ['pictures_in_paper/melanoma/melanoma_', genes{i}, '_section'], [0,90], [-inf,inf],[-inf,inf],[-5,5],genes(i));
end
%% melanoma_down_CD63
mkdir('pictures_in_paper/melanoma_down_CD63/');
draw_result3d_or_split_or_gif_with_gramm(melanoma_down_CD63, 'pictures_in_paper/melanoma_down_CD63/melanoma_down_CD63_3d_global', 'normal');
draw_result3d_with_gramm(melanoma_down_CD63, 'pictures_in_paper/melanoma_down_CD63/melanoma_down_CD63_3d_views', 0);
draw_sections_with_gramm(melanoma_down_CD63, 'pictures_in_paper/melanoma_down_CD63/melanoma_down_CD63_sections_normal', 'normal');
draw_sections_with_gramm(melanoma_down_CD63, 'pictures_in_paper/melanoma_down_CD63/melanoma_down_CD63_sections_density', 'density');
draw_bar_of_connection_number_with_gramm(melanoma_down_CD63, 'pictures_in_paper/melanoma_down_CD63/melanoma_down_CD63_connection_number', 0);
draw_bar_of_connection_number_with_gramm(melanoma_down_CD63, 'pictures_in_paper/melanoma_down_CD63/melanoma_down_CD63_connection_number_normalized', 1);
draw_density_with_gramm(melanoma_down_CD63, 'pictures_in_paper/melanoma_down_CD63/melanoma_down_CD63_density',8,1.5);
draw_qvalue_with_gramm(melanoma_down_CD63, 'pictures_in_paper/melanoma_down_CD63/melanoma_down_CD63_qvalue', melanoma.standards, 10, 1, 15, 0.2, 0.2);
draw_compare_density_with_gramm(melanoma,melanoma_down_CD63,[], 'origin', 'CD63_down_expression', 'pictures_in_paper/melanoma_down_CD63/melanoma_down_CD63_density_compare',1,0.1);
draw_compare_of_connection_number_with_gramm(melanoma,melanoma_down_CD63, {'malignant', 'Macrophage'},{'malignant', 'Macrophage'}, 'origin', 'CD63_down_expression', 'pictures_in_paper/melanoma_down_CD63/melanoma_down_CD63_malignant_Macrophage_compare','normalized_number');
%% melanoma_cell_ICR: treated melanoma
mkdir('pictures_in_paper/melanoma_cell_ICR/');
draw_result3d_or_split_or_gif_with_gramm(melanoma_cell_ICR, 'pictures_in_paper/melanoma_cell_ICR/melanoma_cell_ICR_3d_global', 'normal');
draw_result3d_with_gramm(melanoma_cell_ICR, 'pictures_in_paper/melanoma_cell_ICR/melanoma_cell_ICR_3d_views', 0);
draw_sections_with_gramm(melanoma_cell_ICR, 'pictures_in_paper/melanoma_cell_ICR/melanoma_cell_ICR_sections_normal', 'normal');
draw_sections_with_gramm(melanoma_cell_ICR, 'pictures_in_paper/melanoma_cell_ICR/melanoma_cell_ICR_sections_density', 'density');
draw_bar_of_connection_number_with_gramm(melanoma_cell_ICR, 'pictures_in_paper/melanoma_cell_ICR/melanoma_cell_ICR_connection_number', 0);
draw_bar_of_connection_number_with_gramm(melanoma_cell_ICR, 'pictures_in_paper/melanoma_cell_ICR/melanoma_cell_ICR_connection_number_normalized', 1);
draw_density_with_gramm(melanoma_cell_ICR, 'pictures_in_paper/melanoma_cell_ICR/melanoma_cell_ICR_density',8,1.5);
draw_qvalue_with_gramm(melanoma_cell_ICR, 'pictures_in_paper/melanoma_cell_ICR/melanoma_cell_ICR_qvalue', melanoma_cell_ICR.standards, 8, 0.8, 15, 0.22, 0.22);
draw_compare_of_connection_number_with_gramm(melanoma_cell_ICR,melanoma, {'maligant-ICR', 'T.CD8-ICR'},{'malignant', 'T-cell'}, 'melanoma_ICR', 'melanoma', 'pictures_in_paper/melanoma_cell_ICR/melanoma_cell_ICR_melanoma_compare','normalized_number');
%% pseudo_HNC
% note this section sometimes could meet an error because of the fit
% function in matlab. This fit function uses random initial solution and 
% sometimes it could not give correct result if the initial solution is bad,
% just ignore it and run it again.
load('data/analysts/workspace_HNC.mat')
draw_pseudo_with_gramm(tumor_connection/2215,blood_connection/2215,'Tumor','Blood','pictures_in_paper/HNC/HNC_pseudo_connection','linear');
draw_pseudo_with_gramm(tumorT/2215,bloodT/2215,'Tumor','Blood','pictures_in_paper/HNC/HNC_pseudo_T_cellnumber','log');
draw_pseudo_with_gramm(tumormali/2215,bloodmali/2215,'Tumor','Blood','pictures_in_paper/HNC/HNC_pseudo_mali_cellnumber','log');
HNC_bloodT=bloodT;
HNC_tumorT=tumorT;
HNC_bloodmali=bloodmali;
HNC_tumormali=tumormali;
HNC_bloodconnection=blood_connection;
HNC_tumorconnection=tumor_connection;
%% pseudo_melanoma
load('data/analysts/workspace_melanoma.mat')
draw_pseudo_with_gramm(tumor_connection/1257,blood_connection/1257,'Tumor','Blood','pictures_in_paper/melanoma/melanoma_pseudo_connection','linear');
draw_pseudo_with_gramm(tumorT/1257,bloodT/1257,'Tumor','Blood','pictures_in_paper/melanoma/melanoma_pseudo_T_cellnumber','log');
draw_pseudo_with_gramm(tumormali/1257,bloodmali/1257,'Tumor','Blood','pictures_in_paper/melanoma/melanoma_pseudo_mali_cellnumber','log');
melanoma_bloodT=bloodT;
melanoma_tumorT=tumorT;
melanoma_bloodmali=bloodmali;
melanoma_tumormali=tumormali;
melanoma_bloodconnection=blood_connection;
melanoma_tumorconnection=tumor_connection;
%% compare
draw_pseudo_with_gramm(melanoma_tumorT(1:25,:)/1257,HNC_tumorT/2215,'melanoma','HNC', 'pictures_in_paper/HNC/HNC_pseudo_compare_melanoma_tumorT', 'log')
draw_pseudo_with_gramm(melanoma_tumormali(1:25,:)/1257,HNC_tumormali/2215,'melanoma','HNC', 'pictures_in_paper/HNC/HNC_pseudo_compare_melanoma_tumormali', 'log')
draw_pseudo_with_gramm(melanoma_bloodT(1:25,:)/1257,HNC_bloodT/2215,'melanoma','HNC', 'pictures_in_paper/HNC/HNC_pseudo_compare_melanoma_bloodT', 'log')
draw_pseudo_with_gramm(melanoma_bloodmali(1:25,:)/1257,HNC_bloodmali/2215,'melanoma','HNC', 'pictures_in_paper/HNC/HNC_pseudo_compare_melanoma_bloodmali', 'log')
draw_pseudo_with_gramm(melanoma_bloodconnection(1:25,:)/1257,HNC_bloodconnection/2215,'melanoma','HNC', 'pictures_in_paper/HNC/HNC_pseudo_compare_melanoma_bloodconnection', 'linear')
draw_pseudo_with_gramm(melanoma_tumorconnection(1:25,:)/1257,HNC_tumorconnection/2215,'melanoma','HNC', 'pictures_in_paper/HNC/HNC_pseudo_compare_melanoma_tumorconnection', 'linear')
%% correlation genes
%#ok<*ASGLU,*NOPTS>
% IFNG CCL4 in LAYN
[tbl1,tbl2,p1,p2]=correlation_two_genes(liver,'CD8-LAYN','CD4-CTLA4',{'IFNG','CCL4'},[1,4000], 'liver Tex', 'pictures_in_paper/liver/liver_LAYN_corr_IFNG_CCL4')
[tbl1,tbl2,p1,p2]=correlation_two_genes(lung,'CD8_C6-LAYN','CD4_C9-CTLA4',{'IFNG','CCL4'},[1,4000],'lung Tex', 'pictures_in_paper/lung/lung_LAYN_corr_IFNG_CCL4')
[tbl1,tbl2,p1,p2]=correlation_two_genes(CRC,'CD8_C07-LAYN','CD4_C12-CTLA4',{'IFNG','CCL4'},[1,4000], 'CRC Tex', 'pictures_in_paper/CRC/CRC_LAYN_corr_IFNG_CCL4')

% TNF CCL4 in LAYN
[tbl1,tbl2,p1,p2]=correlation_two_genes(liver,'CD8-LAYN','CD4-CTLA4',{'TNF','CCL4'}, [1,4000],'liver Tex','pictures_in_paper/liver/liver_LAYN_corr_TNF_CCL4') 
[tbl1,tbl2,p1,p2]=correlation_two_genes(lung,'CD8_C6-LAYN','CD4_C9-CTLA4',{'TNF','CCL4'},[1,4000],'lung Tex','pictures_in_paper/lung/lung_LAYN_corr_TNF_CCL4')
[tbl1,tbl2,p1,p2]=correlation_two_genes(CRC,'CD8_C07-LAYN','CD4_C12-CTLA4',{'TNF','CCL4'},[1,4000], 'CRC Tex', 'pictures_in_paper/CRC/CRC_LAYN_corr_TNF_CCL4')

% IFNG CCL4 in GZMK
[tbl1,tbl2,p1,p2]=correlation_two_genes(liver,'CD8-GZMK','CD4-CTLA4',{'IFNG','CCL4'}, [1,4000],'liver Tem', 'pictures_in_paper/liver/liver_GZMK_corr_IFNG_CCL4')
[tbl1,tbl2,p1,p2]=correlation_two_genes(lung,'CD8_C4-GZMK','CD4_C9-CTLA4',{'IFNG','CCL4'},[1,4000],'lung Tem','pictures_in_paper/lung/lung_GZMK_corr_IFNG_CCL4')
[tbl1,tbl2,p1,p2]=correlation_two_genes(CRC,'CD8_C04-GZMK','CD4_C12-CTLA4',{'IFNG','CCL4'}, [1,4000],'CRC Tem', 'pictures_in_paper/CRC/CRC_GZMK_corr_IFNG_CCL4')

% TNF CCL4 in GZMK
[tbl1,tbl2,p1,p2]=correlation_two_genes(liver,'CD8-GZMK','CD4-CTLA4',{'TNF','CCL4'}, [1,4000],'liver Tem', 'pictures_in_paper/liver/liver_GZMK_corr_TNF_CCL4')
[tbl1,tbl2,p1,p2]=correlation_two_genes(lung,'CD8_C4-GZMK','CD4_C9-CTLA4',{'TNF','CCL4'},[1,4000],'lung Tem','pictures_in_paper/lung/lung_GZMK_corr_TNF_CCL4')
[tbl1,tbl2,p1,p2]=correlation_two_genes(CRC,'CD8_C04-GZMK','CD4_C12-CTLA4',{'TNF','CCL4'}, [1,4000],'CRC Tem','pictures_in_paper/CRC/CRC_GZMK_corr_TNF_CCL4')

%% Human pancrea
load('data/analysts/human_pancrea.mat');
draw_for_one_dataset(a, 'human_pancrea', 'pictures_in_paper/human_pancrea/');

%% Mouse pancrea
load('data/analysts/mouse_pancrea.mat');
draw_for_one_dataset(a, 'mouse_pancrea', 'pictures_in_paper/mouse_pancrea/');

%% PMID30283141
analysis_for_PMID30283141;

%% GSE84498
analysis_for_GSE84498;

%% GSE108561
analysis_for_GSE108561;

%% human placenta
load('data/analysts/human_placenta.mat');
order = {'DC1';'Granulocytes';'ILC3';'Plasma_cells';'dM1';'dM2';'dNK1';'dNK2';'dNK3';'dNK_p';'dS1';'dS2';'dS3';'dT_cells';'T_cells';'DC2';'Monocytes';'NK';'NK_CD16+';'EVT';'EVT_p';'F1';'F2';'HB';'M3';'SCT';'VCT';'VCT_p';'un_annoted'};
draw_for_one_dataset(a, 'human_placenta', 'pictures_in_paper/human_placenta/');
draw_qvalue_with_gramm(a, ['pictures_in_paper/human_placenta/', 'human_placenta', '_qvalue'], order, 200/size(a.standards,1), 15/(size(a.standards,1)^2), 11, 0.15, 0.15);

%% random LR corr
load('data/qiming.IHC.random.LR.norm.corr.mat')
g(1,1) = gramm('x', allcorrs);
g(1,2) = gramm('x', allcorrs);
g(1,1).stat_bin('geom','stacked_bar')
g(1,2).stat_density();
g.set_names('x', 'correlations', 'y', 'counts');
g.set_text_options('base_size', 30);
f = figure();
set(f,'units','normalized','position',[0, 0, 0.6, 0.4]);
g.axe_property('TickDir', 'out');
g.draw();
g(1,2).update('x', ones(10,1)*realcorr, 'y', 0.2:0.2:2);
g(1,2).geom_line();
g(1,2).set_color_options('map', 'matlab');
g(1,2).set_line_options('styles', {'--'});
g(1,2).draw();
g.export('file_name', ['pictures_in_paper/compare.random.LR/compare.connection.number.random.LR.corr', '.pdf'], 'file_type', 'pdf');
%% random coord corr
load('data/qiming.IHC.random.coord.norm.corr.mat')
g(1,1) = gramm('x', allcorrs);
g(1,2) = gramm('x', allcorrs);
g(1,1).stat_bin('geom','stacked_bar')
g(1,2).stat_density();
g.set_names('x', 'correlations', 'y', 'counts');
g.set_text_options('base_size', 30);
f = figure();
set(f,'units','normalized','position',[0, 0, 0.6, 0.4]);
g.axe_property('TickDir', 'out');
g.draw();
g(1,2).update('x', ones(10,1)*realcorr, 'y', 0.2:0.2:2);
g(1,2).geom_line();
g(1,2).set_color_options('map', 'matlab');
g(1,2).set_line_options('styles', {'--'});
g(1,2).draw();
g.export('file_name', ['pictures_in_paper/compare.random.coord/compare.connection.number.random.coord.corr', '.pdf'], 'file_type', 'pdf');
%% corr sumLR and aff
correlation_aff_sumLR;

%% robust_dim_cutoff.m
robust_dim_cutoff;

%% robust_cluster_dp.m
robust_cluster_dp;

%% correlation_cellphoneDB.m
correlation_cellphoneDB;
