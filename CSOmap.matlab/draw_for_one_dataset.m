function ignored = draw_for_one_dataset(a, cancertype, outputpath) %#ok<STOUT>
% this function draw the basic pictures of dataset a, to use advance
% functions to study a certain gene, you need to turn to other functions.
% ATTENTION: those pictures might not be 'perfect', you need to manually
% set the points' size, picutre's size, etc, please refer to each
% function's annotation for details.
addpath('draw_pictures/');
addpath('draw_pictures/gramm-master/');
a.outputpath=outputpath;
a.writeresult3d([cancertype, '_coordinate']);
a.writecounts(cancertype);
a.writestatistics([cancertype, '_statistics']);
a.drawconclusion(0.05, [cancertype, '_conclusion']);
a.savegif([cancertype, '_3dplot'], [cancertype, '_3dplot']);
a.mainLR([cancertype, '_mainLR']);
draw_result3d_or_split_or_gif_with_gramm(a, [outputpath, cancertype, '_3d_global'], 'normal');
draw_result3d_with_gramm(a, [outputpath, cancertype, '_3d_views'], 0);
draw_sections_with_gramm(a, [outputpath, cancertype, '_sections_normal'], 'normal');
draw_sections_with_gramm(a, [outputpath, cancertype, '_sections_density'], 'density');
draw_bar_of_connection_number_with_gramm(a, [outputpath, cancertype, '_connection_number'], 0);
draw_bar_of_connection_number_with_gramm(a, [outputpath, cancertype, '_connection_number_normalized'], 1);
draw_qvalue_with_gramm(a, [outputpath, cancertype, '_qvalue'], a.standards, 150/size(a.standards,1), 15/(size(a.standards,1)^2), 15, 0.2, 0.2);
draw_density_with_gramm(a, [outputpath, cancertype, '_density'],10,2);
draw_result3d_with_section_with_gramm(a, [outputpath, cancertype, '_section_z=0'], [], [0,90], [-inf,inf],[-inf,inf],[-5,5],'normal');
end

