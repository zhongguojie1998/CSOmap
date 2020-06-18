% this script will draw all the pictures in our paper, 
% each function's annotation could be found in that function's m file

%% load data
load('/Users/zhongguojie/Desktop/temporary/analysts/CRC_all.mat');
CRC_NPT = a;
load('/Users/zhongguojie/Desktop/temporary/analysts/liver_all.mat');
liver_NPT = a;
load('/Users/zhongguojie/Desktop/temporary/analysts/lung_all.mat');
lung_NPT = a;
%% CRC_NPT
draw_result3d_or_split_or_gif_with_gramm(CRC_NPT, 'CRC_NPT/CRC_NPT_3d_global', 'normal');
draw_result3d_with_gramm(CRC_NPT, 'CRC_NPT/CRC_NPT_3d_views', 0);
draw_sections_with_gramm(CRC_NPT, 'CRC_NPT/CRC_NPT_sections_normal', 'normal');
draw_sections_with_gramm(CRC_NPT, 'CRC_NPT/CRC_NPT_sections_density', 'density');
draw_bar_of_connection_number_with_gramm(CRC_NPT, 'CRC_NPT/CRC_NPT_connection_number', 0);
draw_bar_of_connection_number_with_gramm(CRC_NPT, 'CRC_NPT/CRC_NPT_connection_number_normalized', 1);
draw_density_with_gramm(CRC_NPT, 'CRC_NPT/CRC_NPT_density',3,0.5);
draw_qvalue_with_gramm(CRC_NPT, 'CRC_NPT/CRC_NPT_qvalue', [], 30, 5, 20, 0.6, 0.6, [0.1,0.1,0.72,0.72]);
%% liver_NPT
draw_result3d_or_split_or_gif_with_gramm(liver_NPT, 'liver_NPT/liver_NPT_3d_global', 'normal');
draw_result3d_with_gramm(liver_NPT, 'liver_NPT/liver_NPT_3d_views', 0);
draw_sections_with_gramm(liver_NPT, 'liver_NPT/liver_NPT_sections_normal', 'normal');
draw_sections_with_gramm(liver_NPT, 'liver_NPT/liver_NPT_sections_density', 'density');
draw_bar_of_connection_number_with_gramm(liver_NPT, 'liver_NPT/liver_NPT_connection_number', 0);
draw_bar_of_connection_number_with_gramm(liver_NPT, 'liver_NPT/liver_NPT_connection_number_normalized', 1);
draw_density_with_gramm(liver_NPT, 'liver_NPT/liver_NPT_density',3,0.5);
draw_qvalue_with_gramm(liver_NPT, 'liver_NPT/liver_NPT_qvalue', [], 30, 5, 20, 0.6, 0.6, [0.1,0.1,0.72,0.72]);
%% lung_NPT
draw_result3d_or_split_or_gif_with_gramm(lung_NPT, 'lung_NPT/lung_NPT_3d_global', 'normal');
draw_result3d_with_gramm(lung_NPT, 'lung_NPT/lung_NPT_3d_views', 0);
draw_sections_with_gramm(lung_NPT, 'lung_NPT/lung_NPT_sections_normal', 'normal');
draw_sections_with_gramm(lung_NPT, 'lung_NPT/lung_NPT_sections_density', 'density');
draw_bar_of_connection_number_with_gramm(lung_NPT, 'lung_NPT/lung_NPT_connection_number', 0);
draw_bar_of_connection_number_with_gramm(lung_NPT, 'lung_NPT/lung_NPT_connection_number_normalized', 1);
draw_density_with_gramm(lung_NPT, 'lung_NPT/lung_NPT_density',2,0.5);
draw_qvalue_with_gramm(lung_NPT, 'lung_NPT/lung_NPT_qvalue', [], 30, 5, 20, 0.6, 0.6, [0.1,0.1,0.72,0.72]);
%% CRC
draw_result3d_or_split_or_gif_with_gramm(CRC, 'CRC/CRC_3d_global', 'normal');
draw_result3d_with_gramm(CRC, 'CRC/CRC_3d_views', 0);
draw_sections_with_gramm(CRC, 'CRC/CRC_sections_normal', 'normal');
draw_sections_with_gramm(CRC, 'CRC/CRC_sections_density', 'density');
draw_bar_of_connection_number_with_gramm(CRC, 'CRC/CRC_connection_number', 0);
draw_bar_of_connection_number_with_gramm(CRC, 'CRC/CRC_connection_number_normalized', 1);
draw_qvalue_with_gramm(CRC, 'CRC/CRC_qvalue', CRC.standards, 1.4, 0.14, 15, 0.25, 0.25);
draw_one_gene_with_gramm(CRC, 'CCR8', 'CRC/CRC_CCR8_TPM', 10, 0.5);
draw_density_with_gramm(CRC, 'CRC/CRC_density',10,2);
draw_result3d_with_section_with_gramm(CRC, 'CRC/CRC_section_z=0', [], [0,90], [-inf,inf],[-inf,inf],[-5,5],'normal');
draw_one_gene_with_gramm(CRC, 'CCR8', 'CRC/CRC_CCR8_TPM', 10, 2);
draw_one_gene_with_gramm(CRC, 'CCL4', 'CRC/CRC_CCL4_TPM', 10, 2);
genes = {'CCR8', 'CCL4', 'MKI67', 'IFNG', 'TNF', 'PDCD1', 'CD274', 'CTLA4', 'CD80', 'CD86'};
for i = 1 : size(genes, 2)
    draw_result3d_with_genes_with_gramm(CRC, ['CRC/CRC_', genes{i}], [30,30], [-inf,inf],[-inf,inf],[-inf,inf],genes(i));
    draw_result3d_with_genes_with_gramm(CRC, ['CRC/CRC_', genes{i}, '_section'], [0,90], [-inf,inf],[-inf,inf],[-5,5],genes(i));
end
%% CRC_knockout_Treg
draw_result3d_or_split_or_gif_with_gramm(CRC_knockout_Treg, 'CRC_knockout_Treg/CRC_knockout_Treg_3d_global', 'normal');
draw_result3d_with_gramm(CRC_knockout_Treg, 'CRC_knockout_Treg/CRC_knockout_Treg_3d_views', 0);
draw_sections_with_gramm(CRC_knockout_Treg, 'CRC_knockout_Treg/CRC_knockout_Treg_sections_normal', 'normal');
draw_sections_with_gramm(CRC_knockout_Treg, 'CRC_knockout_Treg/CRC_knockout_Treg_sections_density', 'density');
draw_bar_of_connection_number_with_gramm(CRC_knockout_Treg, 'CRC_knockout_Treg/CRC_knockout_Treg_connection_number', 0);
draw_bar_of_connection_number_with_gramm(CRC_knockout_Treg, 'CRC_knockout_Treg/CRC_knockout_Treg_connection_number_normalized', 1);
draw_density_with_gramm(CRC_knockout_Treg, 'CRC_knockout_Treg/CRC_knockout_Treg_density',3,0.5);
draw_qvalue_with_gramm(CRC_knockout_Treg, 'CRC_knockout_Treg/CRC_knockout_Treg_qvalue', CRC_knockout_Treg.standards, 1.7, 0.17, 15, 0.25, 0.25);
%% CRC_ziyi
draw_result3d_or_split_or_gif_with_gramm(CRC_ziyi, 'CRC_ziyi/CRC_ziyi_3d_global', 'normal');
draw_result3d_with_gramm(CRC_ziyi, 'CRC_ziyi/CRC_ziyi_3d_views', 0);
draw_sections_with_gramm(CRC_ziyi, 'CRC_ziyi/CRC_ziyi_sections_normal', 'normal');
draw_sections_with_gramm(CRC_ziyi, 'CRC_ziyi/CRC_ziyi_sections_density', 'density');
draw_bar_of_connection_number_with_gramm(CRC_ziyi, 'CRC_ziyi/CRC_ziyi_connection_number', 0);
draw_bar_of_connection_number_with_gramm(CRC_ziyi, 'CRC_ziyi/CRC_ziyi_connection_number_normalized', 1);
draw_qvalue_with_gramm(CRC_ziyi, 'CRC_ziyi/CRC_ziyi_qvalue', CRC_ziyi.standards, 1, 0.1, 15, 0.25, 0.25);
draw_density_with_gramm(CRC_ziyi, 'CRC_ziyi/CRC_ziyi_density',10,2);
genes = {'CCR8', 'CCL4', 'MKI67', 'IFNG', 'TNF', 'PDCD1', 'CD274', 'CTLA4', 'CD80', 'CD86'};
for i = 1 : size(genes, 2)
    draw_result3d_with_genes_with_gramm(CRC_ziyi, ['CRC_ziyi/CRC_ziyi_', genes{i}], [30,30], [-inf,inf],[-inf,inf],[-inf,inf],genes(i));
    draw_result3d_with_genes_with_gramm(CRC_ziyi, ['CRC_ziyi/CRC_ziyi_', genes{i}, '_section'], [0,90], [-inf,inf],[-inf,inf],[-5,5],genes(i));
end
%% liver
draw_result3d_or_split_or_gif_with_gramm(liver, 'liver/liver_3d_global', 'normal');
draw_result3d_with_gramm(liver, 'liver/liver_3d_views', 0);
draw_sections_with_gramm(liver, 'liver/liver_sections_normal', 'normal');
draw_sections_with_gramm(liver, 'liver/liver_sections_density', 'density');
draw_bar_of_connection_number_with_gramm(liver, 'liver/liver_connection_number', 0);
draw_bar_of_connection_number_with_gramm(liver, 'liver/liver_connection_number_normalized', 1);
draw_qvalue_with_gramm(liver, 'liver/liver_qvalue', liver.standards, 4, 0.4, 15, 0.2, 0.2);
draw_density_with_gramm(liver, 'liver/liver_density',10,2);
draw_result3d_with_section_with_gramm(liver, 'liver/liver_section_z=0', [], [0,90], [-inf,inf],[-inf,inf],[-5,5],'normal');
draw_one_gene_with_gramm(liver, 'CCR8', 'liver/liver_CCR8_TPM', 10, 0.5);
genes = {'CCR8', 'CCL4', 'MKI67', 'IFNG', 'TNF', 'PDCD1', 'CD274', 'CTLA4', 'CD80', 'CD86'};
for i = 1 : size(genes, 2)
    draw_result3d_with_genes_with_gramm(liver, ['liver/liver_', genes{i}], [30,30], [-inf,inf],[-inf,inf],[-inf,inf],genes(i));
    draw_result3d_with_genes_with_gramm(liver, ['liver/liver_', genes{i}, '_section'], [0,90], [-inf,inf],[-inf,inf],[-5,5],genes(i));
end
%% lung
draw_result3d_or_split_or_gif_with_gramm(lung, 'lung/lung_3d_global', 'normal');
draw_result3d_with_gramm(lung, 'lung/lung_3d_views', 0);
draw_sections_with_gramm(lung, 'lung/lung_sections_normal', 'normal');
draw_sections_with_gramm(lung, 'lung/lung_sections_density', 'density');
draw_bar_of_connection_number_with_gramm(lung, 'lung/lung_connection_number', 0);
draw_bar_of_connection_number_with_gramm(lung, 'lung/lung_connection_number_normalized', 1);
draw_qvalue_with_gramm(lung, 'lung/lung_qvalue', lung.standards, 1.5, 0.15, 15, 0.2, 0.2);
draw_one_gene_with_gramm(lung, 'CCR8', 'lung/lung_CCR8_TPM', 10, 1.5);
draw_density_with_gramm(lung, 'lung/lung_density',10,2);
draw_result3d_with_section_with_gramm(lung, 'lung/lung_section_z=0', [], [0,90], [-inf,inf],[-inf,inf],[-5,5],'normal');
genes = {'CCR8', 'CCL4', 'MKI67', 'IFNG', 'TNF', 'PDCD1', 'CD274', 'CTLA4', 'CD80', 'CD86'};
for i = 1 : size(genes, 2)
    draw_result3d_with_genes_with_gramm(lung, ['lung/lung_', genes{i}], [30,30], [-inf,inf],[-inf,inf],[-inf,inf],genes(i));
    draw_result3d_with_genes_with_gramm(lung, ['lung/lung_', genes{i}, '_section'], [0,90], [-inf,inf],[-inf,inf],[-5,5],genes(i));
end
%% 3D heatmaps of CCR8, CCL4, MKI67, IFNG, TNF, PDCD1, CD274, CTLA4, CD80, CD86
genes = {'CCR8', 'CCL4', 'MKI67', 'IFNG', 'TNF', 'PDCD1', 'CD274', 'CTLA4', 'CD80', 'CD86'};
for i = 1 : size(genes, 2)
    draw_result3d_with_genes_with_gramm(CRC, ['CRC/CRC_', genes{i}], [30,30], [-inf,inf],[-inf,inf],[-inf,inf],genes(i));
    draw_result3d_with_genes_with_gramm(liver, ['liver/liver_', genes{i}], [30,30], [-inf,inf],[-inf,inf],[-inf,inf],genes(i));
    draw_result3d_with_genes_with_gramm(lung, ['lung/lung_', genes{i}], [30,30], [-inf,inf],[-inf,inf],[-inf,inf],genes(i));
end
%% sections of CCR8, CCL4, MKI67, IFNG, TNF, PDCD1, CD274, CTLA4, CD80, CD86
genes = {'CCR8', 'CCL4', 'MKI67', 'IFNG', 'TNF', 'PDCD1', 'CD274', 'CTLA4', 'CD80', 'CD86'};
for i = 1 : size(genes, 2)
    draw_result3d_with_genes_with_gramm(CRC, ['CRC/CRC_', genes{i}, '_section'], [0,90], [-inf,inf],[-inf,inf],[-5,5],genes(i));
    draw_result3d_with_genes_with_gramm(liver, ['liver/liver_', genes{i}, '_section'], [0,90], [-inf,inf],[-inf,inf],[-5,5],genes(i));
    draw_result3d_with_genes_with_gramm(lung, ['lung/lung_', genes{i}, '_section'], [0,90], [-inf,inf],[-inf,inf],[-5,5],genes(i));
end
%% liver knock out CCR8
draw_result3d_or_split_or_gif_with_gramm(liver_knockout_CCR8, 'liver_knockout_CCR8/liver_knockout_CCR8_3d_global', 'normal');
draw_result3d_with_gramm(liver_knockout_CCR8, 'liver_knockout_CCR8/liver_knockout_CCR8_3d_views', 0);
draw_sections_with_gramm(liver_knockout_CCR8, 'liver_knockout_CCR8/liver_knockout_CCR8_sections_normal', 'normal');
draw_sections_with_gramm(liver_knockout_CCR8, 'liver_knockout_CCR8/liver_knockout_CCR8_sections_density', 'density');
draw_bar_of_connection_number_with_gramm(liver_knockout_CCR8, 'liver_knockout_CCR8/liver_knockout_CCR8_connection_number', 0);
draw_bar_of_connection_number_with_gramm(liver_knockout_CCR8, 'liver_knockout_CCR8/liver_knockout_CCR8_connection_number_normalized', 1);
draw_qvalue_with_gramm(liver_knockout_CCR8, 'liver_knockout_CCR8/liver_knockout_CCR8_qvalue', liver_knockout_CCR8.standards, 4, 0.4, 15, 0.2, 0.2);
draw_density_with_gramm(liver_knockout_CCR8, 'liver_knockout_CCR8/liver_knockout_CCR8_density',10,2);
genes = {'CCR8', 'CCL4', 'MKI67', 'IFNG', 'TNF', 'PDCD1', 'CD274', 'CTLA4', 'CD80', 'CD86'};
for i = 1 : size(genes, 2)
    draw_result3d_with_genes_with_gramm(liver_knockout_CCR8, ['liver_knockout_CCR8/liver_knockout_CCR8_', genes{i}], [30,30], [-inf,inf],[-inf,inf],[-inf,inf],genes(i));
    draw_result3d_with_genes_with_gramm(liver_knockout_CCR8, ['liver_knockout_CCR8/liver_knockout_CCR8_', genes{i}, '_section'], [0,90], [-inf,inf],[-inf,inf],[-5,5],genes(i));
end

%% HNC
draw_result3d_or_split_or_gif_with_gramm(HNC, 'HNC/HNC_3d_global', 'normal');
draw_result3d_with_gramm(HNC, 'HNC/HNC_3d_views', 0);
draw_sections_with_gramm(HNC, 'HNC/HNC_sections_normal', 'normal');
draw_sections_with_gramm(HNC, 'HNC/HNC_sections_density', 'density');
draw_bar_of_connection_number_with_gramm(HNC, 'HNC/HNC_connection_number', 0);
draw_bar_of_connection_number_with_gramm(HNC, 'HNC/HNC_connection_number_normalized', 1);
draw_density_with_gramm(HNC, 'HNC/HNC_density',3,0.5);
draw_qvalue_with_gramm(HNC, 'HNC/HNC_qvalue', HNC.standards, 10, 1, 15, 0.2, 0.2);
draw_result3d_with_section_with_gramm(HNC, 'HNC/HNC_section_z=0', [], [0,90], [-inf,inf],[-inf,inf],[-5,5]);
genes={'CD63','TIMP1','PDPN','LAMC2','LAMB3'};
%%
for i = 1 : size(genes, 2)
draw_result3d_with_genes_with_gramm(HNC, ['HNC/HNC_', genes{i}], [-35,30], [-inf,inf],[-inf,inf],[-inf,inf],genes(i));
draw_result3d_with_genes_with_gramm(HNC, ['HNC/HNC_', genes{i}, '_section'], [0,90], [-inf,inf],[-inf,inf],[-5,5],genes(i));
draw_compare_one_gene_with_gramm(HNC, melanoma, 'HNC',  'malignant', 'melanoma', 'malignant', genes{i}, ['HNC/HNC_melanoma_compare_', genes{i}], 1, 0.2);
end
%% HNC_pEMT
draw_result3d_or_split_or_gif_with_gramm(HNC_pEMT, 'HNC_pEMT/HNC_pEMT_3d_global', 'normal');
draw_result3d_with_gramm(HNC_pEMT, 'HNC_pEMT/HNC_pEMT_3d_views', 0);
draw_result3d_with_section_with_gramm(HNC_pEMT, 'HNC_pEMT/HNC_pEMT_section_z=0_normal', [], [0,90], [-inf,inf],[-inf,inf],[-5,5], 'normal');
draw_result3d_with_section_with_gramm(HNC_pEMT, 'HNC_pEMT/HNC_pEMT_section_z=0_desity', [], [0,90], [-inf,inf],[-inf,inf],[-5,5], 'density');
draw_sections_with_gramm(HNC_pEMT, 'HNC_pEMT/HNC_pEMT_sections_normal', 'normal');
draw_sections_with_gramm(HNC_pEMT, 'HNC_pEMT/HNC_pEMT_sections_density', 'density');
draw_bar_of_connection_number_with_gramm(HNC_pEMT, 'HNC_pEMT/HNC_pEMT_connection_number', 0);
draw_bar_of_connection_number_with_gramm(HNC_pEMT, 'HNC_pEMT/HNC_pEMT_connection_number_normalized', 1);
draw_density_with_gramm(HNC_pEMT, 'HNC_pEMT/HNC_pEMT_density',3,0.5);
draw_qvalue_with_gramm(HNC_pEMT, 'HNC_pEMT/HNC_pEMT_qvalue', HNC_pEMT.standards, 8, 0.8, 20, 0.2, 0.2);
draw_qvalue_with_gramm(HNC_pEMT_origin, 'HNC_pEMT/HNC_pEMT_origin_qvalue', HNC_pEMT.standards, 8, 0.8, 20, 0.2, 0.2);
%% HNC_pEMT_new
draw_result3d_or_split_or_gif_with_gramm(HNC_pEMT_new, 'HNC_pEMT_new/HNC_pEMT_new_3d_global', 'normal');
draw_result3d_with_gramm(HNC_pEMT_new, 'HNC_pEMT_new/HNC_pEMT_new_3d_views', 0);
draw_result3d_with_section_with_gramm(HNC_pEMT_new, 'HNC_pEMT_new/HNC_pEMT_new_section_z=0_normal', [], [0,90], [-inf,inf],[-inf,inf],[-2,2], 'normal');
draw_result3d_with_section_with_gramm(HNC_pEMT_new, 'HNC_pEMT_new/HNC_pEMT_new_section_z=0_desity', [], [0,90], [-inf,inf],[-inf,inf],[-2,2], 'density');
draw_sections_with_gramm(HNC_pEMT_new, 'HNC_pEMT_new/HNC_pEMT_new_sections_normal', 'normal');
draw_sections_with_gramm(HNC_pEMT_new, 'HNC_pEMT_new/HNC_pEMT_new_sections_density', 'density');
draw_bar_of_connection_number_with_gramm(HNC_pEMT_new, 'HNC_pEMT_new/HNC_pEMT_new_connection_number', 0);
draw_bar_of_connection_number_with_gramm(HNC_pEMT_new, 'HNC_pEMT_new/HNC_pEMT_new_connection_number_normalized', 1);
draw_density_with_gramm(HNC_pEMT_new, 'HNC_pEMT_new/HNC_pEMT_new_density',3,0.5);
draw_qvalue_with_gramm(HNC_pEMT_new, 'HNC_pEMT_new/HNC_pEMT_new_qvalue', HNC_pEMT_new.standards, 10, 1, 15, 0.6, 0.6);
draw_compare_of_connection_number_with_gramm(HNC_pEMT_new, [], {'malignant', 'Fibroblast'}, {'p-EMT cells', 'Fibroblast'}, 'malignant', 'p-EMT cells', 'HNC_pEMT_new/HNC_pEMT_malignant_p-EMT_compare', 'density');
%% HNC_pEMT_origin_new
draw_bar_of_connection_number_with_gramm(HNC_pEMT_origin_new, 'HNC_pEMT_origin_new/HNC_pEMT_origin_new_connection_number', 0);
draw_bar_of_connection_number_with_gramm(HNC_pEMT_origin_new, 'HNC_pEMT_origin_new/HNC_pEMT_origin_new_connection_number_normalized', 1);
draw_density_with_gramm(HNC_pEMT_origin_new, 'HNC_pEMT_origin_new/HNC_pEMT_origin_new_density',3,0.5);
draw_qvalue_with_gramm(HNC_pEMT_origin_new, 'HNC_pEMT_origin_new/HNC_pEMT_origin_new_qvalue', HNC_pEMT_origin_new.standards, 10, 1, 15, 0.6, 0.6);
draw_compare_density_with_gramm(HNC_pEMT_origin_new, HNC_pEMT_new, {'p-EMT cells','Fibroblast', 'malignant'}, 'density(not embeded)', 'density(embeded)', 'HNC_pEMT_origin_new/HNC_pEMT_before&after_embed_density_compare',1,0.1, 'pair');
draw_compare_of_connection_number_with_gramm(HNC_pEMT_origin_new, HNC_pEMT_new, {'p-EMT cells', 'p-EMT cells'}, {'p-EMT cells', 'p-EMT cells'}, 'not embeded', 'embeded', 'HNC_pEMT_origin_new/HNC_pEMT_before&after_mapping_p-EMT_compare', 'density');
draw_compare_of_connection_number_with_gramm(HNC_pEMT_origin_new, HNC_pEMT_new, {'p-EMT cells', 'Fibroblast'}, {'p-EMT cells', 'Fibroblast'}, 'not embeded', 'embeded', 'HNC_pEMT_origin_new/HNC_pEMT_before&after_mapping_p-EMT_Fibroblast_compare', 'density');
%% HNC_pEMT_CD63
draw_result3d_or_split_or_gif_with_gramm(HNC_pEMT_CD63, 'HNC_pEMT_CD63/HNC_pEMT_CD63_3d_global', 'normal');
draw_result3d_with_gramm(HNC_pEMT_CD63, 'HNC_pEMT_CD63/HNC_pEMT_CD63_3d_views', 0);
draw_sections_with_gramm(HNC_pEMT_CD63, 'HNC_pEMT_CD63/HNC_pEMT_CD63_sections_normal', 'normal');
draw_sections_with_gramm(HNC_pEMT_CD63, 'HNC_pEMT_CD63/HNC_pEMT_CD63_sections_density', 'density');
draw_bar_of_connection_number_with_gramm(HNC_pEMT_CD63, 'HNC_pEMT_CD63/HNC_pEMT_CD63_connection_number', 0);
draw_bar_of_connection_number_with_gramm(HNC_pEMT_CD63, 'HNC_pEMT_CD63/HNC_pEMT_CD63_connection_number_normalized', 1);
draw_density_with_gramm(HNC_pEMT_CD63, 'HNC_pEMT_CD63/HNC_pEMT_CD63_density',3,0.5);
draw_qvalue_with_gramm(HNC_pEMT_CD63, 'HNC_pEMT_CD63/HNC_pEMT_CD63_qvalue', HNC_pEMT_CD63.standards, 10, 1, 15, 0.6, 0.6);
draw_compare_density_with_gramm(HNC_pEMT_new, HNC_pEMT_CD63, [], 'origin', 'CD63_over_expression', 'HNC_pEMT_CD63/HNC_pEMT_CD63_density_compare',1,0.1);
%% melanoma
draw_result3d_or_split_or_gif_with_gramm(melanoma, 'melanoma/melanoma_3d_global', 'normal');
draw_result3d_with_gramm(melanoma, 'melanoma/melanoma_3d_views', 0);
draw_sections_with_gramm(melanoma, 'melanoma/melanoma_sections_normal', 'normal');
draw_sections_with_gramm(melanoma, 'melanoma/melanoma_sections_density', 'density');
draw_bar_of_connection_number_with_gramm(melanoma, 'melanoma/melanoma_connection_number', 0);
draw_bar_of_connection_number_with_gramm(melanoma, 'melanoma/melanoma_connection_number_normalized', 1);
draw_density_with_gramm(melanoma, 'melanoma/melanoma_density',8,1.5);
draw_qvalue_with_gramm(melanoma, 'melanoma/melanoma_qvalue', melanoma.standards, 10, 1, 20, 0.2, 0.2);
draw_result3d_with_section_with_gramm(melanoma, 'melanoma/melanoma_section_z=0', [], [0,90], [-inf,inf],[-inf,inf],[-2,2], 'normal');
draw_result3d_with_section_with_gramm(melanoma, 'melanoma/melanoma_section_z=0', [], [0,90], [-inf,inf],[-inf,inf],[-2,2], 'density');
%%
genes={'CD63','TIMP1','PDPN','LAMC2','LAMB3'};
for i = 1 : size(genes, 2)
draw_result3d_with_genes_with_gramm(melanoma, ['melanoma/melanoma_', genes{i}], [-35,30], [-inf,inf],[-inf,inf],[-inf,inf],genes(i));
draw_result3d_with_genes_with_gramm(melanoma, ['melanoma/melanoma_', genes{i}, '_section'], [0,90], [-inf,inf],[-inf,inf],[-5,5],genes(i));
end
%% melanoma_down_CD63
draw_result3d_or_split_or_gif_with_gramm(melanoma_down_CD63, 'melanoma_down_CD63/melanoma_down_CD63_3d_global', 'normal');
draw_result3d_with_gramm(melanoma_down_CD63, 'melanoma_down_CD63/melanoma_down_CD63_3d_views', 0);
draw_sections_with_gramm(melanoma_down_CD63, 'melanoma_down_CD63/melanoma_down_CD63_sections_normal', 'normal');
draw_sections_with_gramm(melanoma_down_CD63, 'melanoma_down_CD63/melanoma_down_CD63_sections_density', 'density');
draw_bar_of_connection_number_with_gramm(melanoma_down_CD63, 'melanoma_down_CD63/melanoma_down_CD63_connection_number', 0);
draw_bar_of_connection_number_with_gramm(melanoma_down_CD63, 'melanoma_down_CD63/melanoma_down_CD63_connection_number_normalized', 1);
draw_density_with_gramm(melanoma_down_CD63, 'melanoma_down_CD63/melanoma_down_CD63_density',8,1.5);
draw_qvalue_with_gramm(melanoma_down_CD63, 'melanoma_down_CD63/melanoma_down_CD63_qvalue', melanoma.standards, 10, 1, 15, 0.2, 0.2);
draw_compare_density_with_gramm(melanoma,melanoma_down_CD63,[], 'origin', 'CD63_down_expression', 'melanoma_down_CD63/melanoma_down_CD63_density_compare',1,0.1);
draw_compare_of_connection_number_with_gramm(melanoma,melanoma_down_CD63, {'malignant', 'Macrophage'},{'malignant', 'Macrophage'}, 'origin', 'CD63_down_expression', 'melanoma_down_CD63/melanoma_down_CD63_malignant_Macrophage_compare','normalized_number');
%% melanoma_cell
draw_result3d_or_split_or_gif_with_gramm(melanoma_cell, 'melanoma_cell/melanoma_cell_3d_global', 'normal');
draw_result3d_with_gramm(melanoma_cell, 'melanoma_cell/melanoma_cell_3d_views', 0);
draw_sections_with_gramm(melanoma_cell, 'melanoma_cell/melanoma_cell_sections_normal', 'normal');
draw_sections_with_gramm(melanoma_cell, 'melanoma_cell/melanoma_cell_sections_density', 'density');
draw_bar_of_connection_number_with_gramm(melanoma_cell, 'melanoma_cell/melanoma_cell_connection_number', 0);
draw_bar_of_connection_number_with_gramm(melanoma_cell, 'melanoma_cell/melanoma_cell_connection_number_normalized', 1);
draw_density_with_gramm(melanoma_cell, 'melanoma_cell/melanoma_cell_density',8,1.5);
draw_qvalue_with_gramm(melanoma_cell, 'melanoma_cell/melanoma_cell_qvalue', melanoma_cell.standards, 8, 0.8, 15, 0.2, 0.2);
%% melanoma_cell_2
draw_result3d_or_split_or_gif_with_gramm(melanoma_cell_2, 'melanoma_cell_2/melanoma_cell_2_3d_global', 'normal');
draw_result3d_with_gramm(melanoma_cell_2, 'melanoma_cell_2/melanoma_cell_2_3d_views', 0);
draw_sections_with_gramm(melanoma_cell_2, 'melanoma_cell_2/melanoma_cell_2_sections_normal', 'normal');
draw_sections_with_gramm(melanoma_cell_2, 'melanoma_cell_2/melanoma_cell_2_sections_density', 'density');
draw_bar_of_connection_number_with_gramm(melanoma_cell_2, 'melanoma_cell_2/melanoma_cell_2_connection_number', 0);
draw_bar_of_connection_number_with_gramm(melanoma_cell_2, 'melanoma_cell_2/melanoma_cell_2_connection_number_normalized', 1);
draw_density_with_gramm(melanoma_cell_2, 'melanoma_cell_2/melanoma_cell_2_density',8,1.5);
draw_qvalue_with_gramm(melanoma_cell_2, 'melanoma_cell_2/melanoma_cell_2_qvalue', melanoma_cell_2.standards, 10, 2, 15, 0.2, 0.2);
%% melanoma_cell_TN
draw_result3d_or_split_or_gif_with_gramm(melanoma_cell_TN, 'melanoma_cell_TN/melanoma_cell_TN_3d_global', 'normal');
draw_result3d_with_gramm(melanoma_cell_TN, 'melanoma_cell_TN/melanoma_cell_TN_3d_views', 0);
draw_sections_with_gramm(melanoma_cell_TN, 'melanoma_cell_TN/melanoma_cell_TN_sections_normal', 'normal');
draw_sections_with_gramm(melanoma_cell_TN, 'melanoma_cell_TN/melanoma_cell_TN_sections_density', 'density');
draw_bar_of_connection_number_with_gramm(melanoma_cell_TN, 'melanoma_cell_TN/melanoma_cell_TN_connection_number', 0);
draw_bar_of_connection_number_with_gramm(melanoma_cell_TN, 'melanoma_cell_TN/melanoma_cell_TN_connection_number_normalized', 1);
draw_density_with_gramm(melanoma_cell_TN, 'melanoma_cell_TN/melanoma_cell_TN_density',8,1.5);
draw_qvalue_with_gramm(melanoma_cell_TN, 'melanoma_cell_TN/melanoma_cell_TN_qvalue', melanoma_cell_TN.standards, 8,0.8, 15, 0.2, 0.2);
%% melanoma_cell_ICR
draw_result3d_or_split_or_gif_with_gramm(melanoma_cell_ICR, 'melanoma_cell_ICR/melanoma_cell_ICR_3d_global', 'normal');
draw_result3d_with_gramm(melanoma_cell_ICR, 'melanoma_cell_ICR/melanoma_cell_ICR_3d_views', 0);
draw_sections_with_gramm(melanoma_cell_ICR, 'melanoma_cell_ICR/melanoma_cell_ICR_sections_normal', 'normal');
draw_sections_with_gramm(melanoma_cell_ICR, 'melanoma_cell_ICR/melanoma_cell_ICR_sections_density', 'density');
draw_bar_of_connection_number_with_gramm(melanoma_cell_ICR, 'melanoma_cell_ICR/melanoma_cell_ICR_connection_number', 0);
draw_bar_of_connection_number_with_gramm(melanoma_cell_ICR, 'melanoma_cell_ICR/melanoma_cell_ICR_connection_number_normalized', 1);
draw_density_with_gramm(melanoma_cell_ICR, 'melanoma_cell_ICR/melanoma_cell_ICR_density',8,1.5);
draw_qvalue_with_gramm(melanoma_cell_ICR, 'melanoma_cell_ICR/melanoma_cell_ICR_qvalue', melanoma_cell_ICR.standards, 8, 0.8, 15, 0.22, 0.22);
draw_compare_of_connection_number_with_gramm(melanoma_cell_ICR,melanoma, {'maligant-ICR', 'T.CD8-ICR'},{'malignant', 'T-cell'}, 'melanoma_ICR', 'melanoma', 'melanoma_cell_ICR/melanoma_cell_ICR_melanoma_compare','normalized_number');
%% pseudo_HNC
load('/Users/zhongguojie/Desktop/temporary/analysts/workspace_HNC.mat')
draw_pseudo_with_gramm(tumor_connection/2215,blood_connection/2215,'Tumor','Blood','HNC/HNC_pseudo_connection','linear');
draw_pseudo_with_gramm(tumorT/2215,bloodT/2215,'Tumor','Blood','HNC/HNC_pseudo_T_cellnumber','log');
draw_pseudo_with_gramm(tumormali/2215,bloodmali/2215,'Tumor','Blood','HNC/HNC_pseudo_mali_cellnumber','log');
HNC_bloodT=bloodT;
HNC_tumorT=tumorT;
HNC_bloodmali=bloodmali;
HNC_tumormali=tumormali;
HNC_bloodconnection=blood_connection;
HNC_tumorconnection=tumor_connection;
%% pseudo_melanoma
load('/Users/zhongguojie/Desktop/temporary/analysts/workspace_melanoma.mat')
draw_pseudo_with_gramm(tumor_connection/1257,blood_connection/1257,'Tumor','Blood','melanoma/melanoma_pseudo_connection','linear');
draw_pseudo_with_gramm(tumorT/1257,bloodT/1257,'Tumor','Blood','melanoma/melanoma_pseudo_T_cellnumber','log');
draw_pseudo_with_gramm(tumormali/1257,bloodmali/1257,'Tumor','Blood','melanoma/melanoma_pseudo_mali_cellnumber','log');
melanoma_bloodT=bloodT;
melanoma_tumorT=tumorT;
melanoma_bloodmali=bloodmali;
melanoma_tumormali=tumormali;
melanoma_bloodconnection=blood_connection;
melanoma_tumorconnection=tumor_connection;
%% compare
draw_pseudo_with_gramm(melanoma_tumorT(1:25,:)/1257,HNC_tumorT/2215,'melanoma','HNC', 'HNC/HNC_pseudo_compare_melanoma_tumorT', 'log')
draw_pseudo_with_gramm(melanoma_tumormali(1:25,:)/1257,HNC_tumormali/2215,'melanoma','HNC', 'HNC/HNC_pseudo_compare_melanoma_tumormali', 'log')
draw_pseudo_with_gramm(melanoma_bloodT(1:25,:)/1257,HNC_bloodT/2215,'melanoma','HNC', 'HNC/HNC_pseudo_compare_melanoma_bloodT', 'log')
draw_pseudo_with_gramm(melanoma_bloodmali(1:25,:)/1257,HNC_bloodmali/2215,'melanoma','HNC', 'HNC/HNC_pseudo_compare_melanoma_bloodmali', 'log')
draw_pseudo_with_gramm(melanoma_bloodconnection(1:25,:)/1257,HNC_bloodconnection/2215,'melanoma','HNC', 'HNC/HNC_pseudo_compare_melanoma_bloodconnection', 'linear')
draw_pseudo_with_gramm(melanoma_tumorconnection(1:25,:)/1257,HNC_tumorconnection/2215,'melanoma','HNC', 'HNC/HNC_pseudo_compare_melanoma_tumorconnection', 'linear')

%% correlation genes
%#ok<*ASGLU,*NOPTS>
% IFNG CCL4 in LAYN
[tbl1,tbl2,p1,p2]=correlation_two_genes(liver,'CD8-LAYN','CD4-CTLA4',{'IFNG','CCL4'},[1,4000], 'liver Tex', 'liver/liver_LAYN_corr_IFNG_CCL4')
[tbl1,tbl2,p1,p2]=correlation_two_genes(lung,'CD8_C6-LAYN','CD4_C9-CTLA4',{'IFNG','CCL4'},[1,4000],'lung Tex', 'lung/lung_LAYN_corr_IFNG_CCL4')
[tbl1,tbl2,p1,p2]=correlation_two_genes(CRC,'CD8_C07-LAYN','CD4_C12-CTLA4',{'IFNG','CCL4'},[1,4000], 'CRC Tex', 'CRC/CRC_LAYN_corr_IFNG_CCL4')

% TNF CCL4 in LAYN
[tbl1,tbl2,p1,p2]=correlation_two_genes(liver,'CD8-LAYN','CD4-CTLA4',{'TNF','CCL4'}, [1,4000],'liver Tex','liver/liver_LAYN_corr_TNF_CCL4') 
[tbl1,tbl2,p1,p2]=correlation_two_genes(lung,'CD8_C6-LAYN','CD4_C9-CTLA4',{'TNF','CCL4'},[1,4000],'lung Tex','lung/lung_LAYN_corr_TNF_CCL4')
[tbl1,tbl2,p1,p2]=correlation_two_genes(CRC,'CD8_C07-LAYN','CD4_C12-CTLA4',{'TNF','CCL4'},[1,4000], 'CRC Tex', 'CRC/CRC_LAYN_corr_TNF_CCL4')

% IFNG CCL4 in GZMK
[tbl1,tbl2,p1,p2]=correlation_two_genes(liver,'CD8-GZMK','CD4-CTLA4',{'IFNG','CCL4'}, [1,4000],'liver Tem', 'liver/liver_GZMK_corr_IFNG_CCL4')
[tbl1,tbl2,p1,p2]=correlation_two_genes(lung,'CD8_C4-GZMK','CD4_C9-CTLA4',{'IFNG','CCL4'},[1,4000],'lung Tem','lung/lung_GZMK_corr_IFNG_CCL4')
[tbl1,tbl2,p1,p2]=correlation_two_genes(CRC,'CD8_C04-GZMK','CD4_C12-CTLA4',{'IFNG','CCL4'}, [1,4000],'CRC Tem', 'CRC/CRC_GZMK_corr_IFNG_CCL4')

% TNF CCL4 in GZMK
[tbl1,tbl2,p1,p2]=correlation_two_genes(liver,'CD8-GZMK','CD4-CTLA4',{'TNF','CCL4'}, [1,4000],'liver Tem', 'liver/liver_GZMK_corr_TNF_CCL4')
[tbl1,tbl2,p1,p2]=correlation_two_genes(lung,'CD8_C4-GZMK','CD4_C9-CTLA4',{'TNF','CCL4'},[1,4000],'lung Tem','lung/lung_GZMK_corr_TNF_CCL4')
[tbl1,tbl2,p1,p2]=correlation_two_genes(CRC,'CD8_C04-GZMK','CD4_C12-CTLA4',{'TNF','CCL4'}, [1,4000],'CRC Tem','CRC/CRC_GZMK_corr_TNF_CCL4')

