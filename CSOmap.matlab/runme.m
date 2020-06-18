function ignored = runme(cancertype, condition) %#ok<STOUT>
    % All the functions for a cancertype, data should be saved in
    % "data/cancertype/" directory, outputs would be saved in
    % "output/cancertype/" directory. condition can be loose or tight, or
    % you can just leave that variable empty and it will automatically be
    % determined.
    % ATTENTION: We set default configurations for the point size, text
    % size, etc.
    % for our datasets, you can run draw_all_pictures.m to
    % reproduce our pictures. For a new dataset, if you are not satisfied
    % with those default configurations, you can refer to each function's
    % annotation for details.
    % These pictures only has the basic information about this dataset, if
    % you want to study specific genes, or compare this dataset with another
    % see draw_all_pictures for an example.
    addpath('draw_pictures/');
    addpath('draw_pictures/gramm-master/');
    if ~exist('condition' ,'var')
        condition = 'tobedetermined';
    end
    % set data path and output path
    datapath = ['data/', cancertype, '/'];
    outputpath = ['output/', cancertype, '/'];
    % pre-process data, save pre-processed data in data.mat
    if ~exist([outputpath, 'data.mat'], 'file')
    	preprocess(datapath, outputpath);
    end
    % reconstruct 3d, save result in workspace.mat
    if ~exist([outputpath, 'workspace.mat'], 'file')
    	reconstruct_3d(outputpath, outputpath, 2, 0, 50, condition);
    end
    % build an object analyst, it will store all the information for analysis
    a = analyst(outputpath, datapath, [outputpath, 'result/'], 1);
    % save that object in a file analyst.mat
    save([outputpath, 'analyst.mat'], 'a', '-v7.3');
    % once analyst is saved, other information is not needed
    delete([outputpath, 'data.mat']);
    delete([outputpath, 'workspace.mat']);
    % use in-built functions to draw pictures and write result files. see
    % each function's annotation for details
    a.writeresult3d([cancertype, '_coordinate']);
    a.writecounts(cancertype);
    a.writestatistics([cancertype, '_statistics']);
    a.drawconclusion(0.05, [cancertype, '_conclusion']);
    a.savegif([cancertype, '_3dplot'], [cancertype, '_3dplot']);
    a.mainLR([cancertype, '_mainLR']);
    % draw pictures using functions in folder draw_pictures/
    outputpath = ['output/', cancertype, '/result/'];
    draw_result3d_or_split_or_gif_with_gramm(a, [outputpath, cancertype, '_3d_global'], 'normal');
    draw_result3d_with_gramm(a, [outputpath, cancertype, '_3d_views'], 0);
    draw_sections_with_gramm(a, [outputpath, cancertype, '_sections_normal'], 'normal');
    draw_sections_with_gramm(a, [outputpath, cancertype, '_sections_density'], 'density');
    draw_bar_of_connection_number_with_gramm(a, [outputpath, cancertype, '_connection_number'], 0);
    draw_bar_of_connection_number_with_gramm(a, [outputpath, cancertype, '_connection_number_normalized'], 1);
    draw_qvalue_with_gramm(a, [outputpath, cancertype, '_qvalue'], a.standards, 150/size(a.standards,1), 15/(size(a.standards,1)^2), 15, 0.2, 0.2);
    draw_density_with_gramm(a, [outputpath, cancertype, '_density'],10,2);
    draw_result3d_with_section_with_gramm(a, [outputpath, cancertype, '_section_z=0'], [], [0,90], [-inf,inf],[-inf,inf],[-5,5],'normal');
    close all;
end
