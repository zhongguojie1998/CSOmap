classdef analyst
    % This class save all the data and provide several methods for
    % statistical analysis and save pictures.
    properties
        TPM                 % TPM data, rows are genes and columns are cells
        cells               % cell names
        genes               % gene names
        affinitymat         % after denoise
        labels              % cell labels (numeric, represent indexs in standards)
        standards           % string, unique labels
        ligands             % ligand names
        receptors           % receptor names
        scores              % weights for each pair of ligand-receptor
        ligandindex         % each ligand's index in this dataset's gene list
        receptorindex       % each receptor's index in this dataset's gene list
        result2d            % 2d coordinates, if possible (default is empty)
        result3d            % 3d coordinates
        process             % coordinates during the optimizing iterations
        connection          % q-val, if connection <= 0.05, then the two clusters have enriched interaction
        reverseconnection   % reverse q-val, if reverseconnection <= 0.05, then the two clusters have depleted interaction
        neighbor            % each row represents one cell and each column represents one cluster, each element means all the interactions between this cell and this cluster
        degree              % each cell's connection number
        counts              % each element includes all the connections between the row cluster and column cluster, use a pair of indices to represent a connection
        clustercounts       % each cluster's number of cells
        outputpath          % default image output path
    end
    methods
        %% main statistical functions
        function obj = analyst(workspace_path, labelpath, finaloutputpath, stat)
            % load data from workspace.mat, and generate an object analyst
            % if stat, then do statistial analysis.
            % ATTENTION, If you don't do stat, some functions can't be used.
            if ~exist(workspace_path, 'dir')
                load(workspace_path); %#ok<LOAD>
            else
                load([workspace_path, 'workspace.mat']); %#ok<LOAD>
            end
            if ~exist('stat', 'var')
                stat = 1;
            end
            [labels, standards] = analyst.identify_label(cells, labelpath);
            % assign attributes
            obj.TPM = TPM;
            obj.genes = genes;
            obj.affinitymat = affinitymat;
            obj.cells = cells;
            obj.labels = labels;
            obj.standards = standards;
            obj.ligands = ligands;
            obj.receptors = receptors;
            obj.scores = scores;
            obj.ligandindex = ligandindex;
            obj.receptorindex = receptorindex;
            obj.result3d = result3d;
            obj.process = process;
            obj.outputpath = finaloutputpath;
            if ~exist(finaloutputpath, 'dir')
               mkdir(finaloutputpath);
            else
               disp('Warning! directory already exists, this program might change the files in it');
            end
            if dim ~= 3
               obj.result2d = result2d;
            end
            if stat
                obj = obj.getconnection(3);
            end
        end
        
        function obj = getconnection(obj, k, dim, useaffinitymat, method, porq)
            % identify the connections between cells, k is the median
            % number of connections, dim indicates the dimension, default
            % is 3, if dim==2, then obj.result2d MUST NOT be empty. You can
            % also use affinitymat instead of coordinates, in that case,
            % you should input affinitymat and assign it to variable
            % 'useaffinitymat'; the default method to calculate q-value is
            % 'hpgdistri', we offer another option, 'permutation'
            if ~exist('dim', 'var')
                dim = 3;
            end
            if ~exist('useaffinitymat','var')
                useaffinitymat = 0;
            end
            if ~exist('method', 'var')
                method = 'hpgdistri';
            end
            if ~exist('porq', 'var')
                % 0 for q, 1 for p.
                porq = 0;
            end
            if dim ~= 3
                tsneresult = obj.result2d;
            else
                tsneresult = obj.result3d;
            end
            
            n = size(tsneresult, 1);
            labelsize = max(obj.labels);
            if labelsize <= 1
                disp('Error: there is only one kind of label in your dataset, unable to do statistical analysis, please check your label.txt!!!')
                return
            end
            if n <= k
                disp('obj not changed, because k is smaller than cell number');
                return
            end
            
            if ~useaffinitymat
                try
                    distancemat = analyst.getdistancemat(tsneresult);            
                    distancemat(1:n+1:end) = Inf;
                    cutoffs = zeros(size(tsneresult, 1), 1);
                    parfor i = 1 : n
                        row = distancemat(i, :);
                        [B, ~] = sort(row);
                        cutoffs(i) = B(k);
                    end
                    cutoff = median(cutoffs);  
                    [row, col, ~] = find(distancemat <= cutoff);
                catch
                    distancemat = sparse(n, n);
                    parfor i = 1 : n
                        point = tsneresult(i, :);
                        row = zeros(n, 1);
                        for j = 1 : n
                            if i == k
                                row(j) = Inf;
                            end
                            row(j) = sqrt(sum((point - tsneresult(j, :)) .^ 2));
                        end
                        [B, ~] = sort(row);
                        cutoffs(i) = B(k);
                        row = 1 ./ row;
                        row(row < 1/B(k)) = 0;
                        row(i) = 0;
                        row = sparse(row);
                        distancemat(:, i) = row;
                    end
                    cutoff = median(cutoffs); 
                    [row, col, ~] = find(distancemat >= 1/cutoff);
                end        
            else
                distancemat = useaffinitymat;
                distancemat(1:n+1:end) = 0;
                cutoffs = zeros(size(tsneresult, 1), 1);
                parfor i = 1 : n
                    row = distancemat(i, :);
                    [B, ~] = sort(row, 'descend');
                    cutoffs(i) = B(k);
                end
                cutoff = median(cutoffs);
                [row, col, ~] = find(distancemat >= cutoff);
            end
            
            % initialize obj counts
            objcounts = cell(labelsize, labelsize);
            parfor i = 1 : labelsize
                for j = 1 : labelsize
                    objcounts{i, j} = zeros([0, 2]);
                end
            end
            % get obj cluster counts
            objclustercounts = tabulate(obj.labels);
            objclustercounts = objclustercounts(:, 2);
            % initialize obj degree
            objneighbor = cell(size(obj.labels, 1), size(obj.standards, 1));
            for i = 1 : size(obj.labels, 1)
                for j = 1 : size(obj.standards, 1)
                    objneighbor{i, j} = zeros([0, 1]);
                end
            end
            % initialize realconnection, realconnection represents connection number
            realconnection = zeros(labelsize, labelsize);
            % update realconnection, obj counts and obj degree
            for i = 1 : size(row, 1)
                realconnection(obj.labels(row(i)), obj.labels(col(i))) = realconnection(obj.labels(row(i)), obj.labels(col(i))) + 1;
                objcounts{obj.labels(row(i)), obj.labels(col(i))} = [objcounts{obj.labels(row(i)), obj.labels(col(i))}; row(i), col(i)];
                objneighbor{row(i), obj.labels(col(i))} = [objneighbor{row(i), obj.labels(col(i))}; col(i)];
            end
            realconnection = realconnection - 1/2 * diag(diag(realconnection));
            if strcmp(method, 'permutation')
                realconnection = realconnection(:);
                randomtotal = zeros([size(realconnection,1),1000]);
                parfor j = 1 : 1000
                    randomconnection = zeros(labelsize, labelsize);
                    randomlabels = obj.labels(randperm(size(obj.labels, 1)), :);
                    for i = 1 : size(row, 1)
                        randomconnection(randomlabels(row(i)), randomlabels(col(i))) = randomconnection(randomlabels(row(i)), randomlabels(col(i))) + 1;
                    end
                    randomconnection = randomconnection - 1/2 * diag(diag(randomconnection));
                    randomtotal(:, j) = randomconnection(:);
                end
                p_value = zeros(size(randomtotal, 1), 1);
                reverse_p_value = zeros(size(randomtotal, 1), 1);
                parfor i = 1 : size(randomtotal, 1)
                    [muhat, sigmahat] = normfit(randomtotal(i, :));
                    if realconnection(i) == 0
                        p_value(i, 1) = 1;
                    else
                        p_value(i, 1) = normcdf(realconnection(i), muhat, sigmahat, 'upper');
                        reverse_p_value(i, 1) = normcdf(realconnection(i), muhat, sigmahat);
                    end
                end
                [~, q_value] = mafdr(p_value);
                [~, reverse_q_value] = mafdr(reverse_p_value);
                q_value = reshape(q_value, [labelsize, labelsize]);
                reverse_q_value = reshape(reverse_q_value, [labelsize, labelsize]);
            elseif strcmp(method, 'hpgdistri')
                pop = size(obj.labels, 1) * (size(obj.labels, 1) - 1) / 2;   % statistical population
                samplenumber = floor(size(row, 1) / 2);
                p_value = zeros([labelsize, labelsize]);
                reverse_p_value = zeros([labelsize, labelsize]);
                for i = 1 : labelsize
                    for j = i : labelsize
                        if i == j
                            propertynumber = objclustercounts(i) * (objclustercounts(i) - 1) / 2;
                        else
                            propertynumber = objclustercounts(i) * objclustercounts(j);
                        end
                        eventnumber = realconnection(i, j);
                        if eventnumber == 0
                            p_value(i, j) = 1;
                            reverse_p_value(i, j) = 0;
                        else
                            p_value(i, j) = 1 - hygecdf(eventnumber, pop, propertynumber, samplenumber);
                            reverse_p_value(i, j) = hygecdf(eventnumber, pop, propertynumber, samplenumber);
                        end
                        p_value(j, i) = p_value(i, j);
                        reverse_p_value(j, i) = reverse_p_value(i, j);
                    end
                end
                [~, q_value] = mafdr(p_value(:));
                [~, reverse_q_value] = mafdr(reverse_p_value(:));
                q_value = reshape(q_value, [labelsize, labelsize]);
                reverse_q_value = reshape(reverse_q_value, [labelsize, labelsize]);
            end
            if porq
                obj.connection = p_value;
                obj.reverseconnection = reverse_p_value;
            else
                obj.connection = q_value;
                obj.reverseconnection = reverse_q_value;
            end
            obj.neighbor = objneighbor;
            % process objcounts, when i==j, has to be corrected, becaues in
            % that case one connection is counted twice
            for i = 1 : size(objcounts, 1)
                count = objcounts{i, i};
                correctedcount = zeros(0, 2);
                countID = {};
                for k = 1 : size(count, 1)
                    if ~ismember([num2str(min(count(k, :))), '-', num2str(max(count(k, :)))], countID)
                        correctedcount = [correctedcount; count(k, :)];
                        countID = [countID; [num2str(min(count(k, :))), '-', num2str(max(count(k, :)))]];
                    end
                end
                objcounts{i, i} = correctedcount;
            end
            obj.counts = objcounts;
            obj.clustercounts = objclustercounts;
            obj.degree = zeros(size(obj.neighbor));
            for i = 1 : size(obj.neighbor, 1)
                for j = 1 : size(obj.neighbor, 2)
                    obj.degree(i, j) = size(obj.neighbor{i, j}, 1);
                end
            end
        end
        
        %% functions about affinitymat
        function f = affinitymathistogram(obj, filename, Title)
            % show affinity distribution, default is not save. 
            % If you want to save this histogram, please specify filename
            % and title, it will be saved in obj.outputpath
            if ~exist('filename', 'var')
                % no filename, don't save
                save = 0;
            else
                save = 1;
            end
            if ~exist('Title', 'var')
                % no filename, don't save
                Title = 'histogram of affinity';
            end
            f = histogram(obj.affinitymat(:));
            title(Title);
            if save
                saveas(f, [obj.outputpath, filename, '.jpg']);
            end
        end
        
        function f = affinitymatshow(obj, normalize, fontsize, filename)
            % show affinitymat, this is a rough version, if you want to
            % draw the picture in our paper, please refer to folder
            % draw_pictures.
            if ~exist('filename', 'var')
                % no filename, don't save
                save = 0;
            else
                save = 1;
            end
            [sorted, rank] = sort(obj.labels);
            Paffinitymat = obj.affinitymat(:, rank);
            Paffinitymat = Paffinitymat(rank, :);
            if normalize
                normalizefactor = diag(sqrt(1 ./ mean(Paffinitymat, 2)));
                Paffinitymat = normalizefactor * Paffinitymat * normalizefactor;
            end
            f = figure;
            imagesc(Paffinitymat);
            colorbar;
            axis off;
            set(f,'units','normalized','position',[0.1 0.1 0.8 0.8]);
            for i = 1 : size(obj.standards, 1)
                [position1, ~] = find(sorted == i);
                position1 = position1(1) - 0.5;
                if i ~= size(obj.standards, 1)
                    [position2, ~] = find(sorted == i+1);
                    position2 = position2(1) - 0.5;
                    position = (position1 + position2)/2;
                else
                    position = (position1 + size(Paffinitymat, 1) + 0.5) / 2;
                end
                line([0, size(Paffinitymat, 1)+0.5], [position1, position1], 'Color', 'red', 'LineWidth', 0.05)
                text(0.5, position, strrep([num2str(i), ':', obj.standards{i}, '-'], '_', '\_'), 'Color', 'red', 'FontSize', fontsize, 'HorizontalAlignment', 'right');
                line([position1, position1], [0, size(Paffinitymat, 1)+0.5], 'Color', 'red', 'LineWidth', 0.05)
                text(position, 0.5, strrep(['-', num2str(i), ':', obj.standards{i}], '_', '\_'), 'Color', 'red', 'FontSize', fontsize, 'HorizontalAlignment', 'left', 'Rotation', 30);
            end
            if save
                saveas(f, [obj.outputpath, filename, '.jpg']);
            end
        end
        
        function ignored = writeaffinitymat(obj, filename) %#ok<STOUT>
            %write affinitymat into txt file
            file = fopen([obj.outputpath, filename, '.txt'], 'w');
            fprintf(file, '%s\t', 'affinity');
            for i = 1 : size(obj.cells, 1)
                fprintf(file, '%s\t', obj.cells{i});
            end
            fprintf(file, '\n');
            for i = 1 : size(obj.affinitymat, 1)
                fprintf(file, '%s\t', obj.cells{i});
                for j = 1 : size(obj.affinitymat, 2)
                    fprintf(file, '%.4f\t', obj.affinitymat(i, j));
                end
                fprintf(file, '\n');
            end
        end
        
        %% functions about result 3d coordinate
        function ignored = writeresult3d(obj, filename) %#ok<STOUT>
            %write coordinates into txt file
            file = fopen([obj.outputpath, filename, '.txt'], 'w');
            fprintf(file, '%s\t', 'ID');
            fprintf(file, '%s\t', 'x');
            fprintf(file, '%s\t', 'y');
            fprintf(file, '%s\t', 'z');
            fprintf(file, '%s\n', 'cluster');
            for i = 1 : size(obj.cells, 1)
                fprintf(file, '%s\t', obj.cells{i});
                fprintf(file, '%f\t', obj.result3d(i, 1));
                fprintf(file, '%f\t', obj.result3d(i, 2));
                fprintf(file, '%f\t', obj.result3d(i, 3));
                fprintf(file, '%f\n', obj.labels(i));
            end
            fclose(file);
        end
        
        %% functions about show process
        function ignored = processshow(obj, filename, show) %#ok<STOUT>
            % show process
            if ~exist('filename', 'var')
                % no filename, then don't save
                save = 0;
                filename = 'process show';
            else
                save = 1;
            end
            if ~exist('show', 'var')
                % if save, then don't show, if not save, then show
                show = ~save;
            end
            iters = size(obj.process, 2)/size(obj.result3d, 2);
            dim = 3;
            for iter = 0 : iters-1
                current_coord = obj.process(:,iter*dim+1:iter*dim+3);
                picture = tdplot(current_coord, obj.labels, obj.standards, filename);
                axis on;
                axis tight;
                if save
                    M=getframe(picture);
                    nn=frame2im(M);
                    [nn,cm]=rgb2ind(nn,256);
                    if iter==0
                        imwrite(nn,cm,[obj.outputpath, filename, '.gif'],'gif','LoopCount',inf,'DelayTime',0.02);%loopcount i==1
                    else
                        imwrite(nn,cm,[obj.outputpath, filename, '.gif'],'gif','WriteMode','append','DelayTime',0.02)%i>=2 loopcount
                    end
                end
                if ~show
                    close all
                end
            end
            function f = tdplot(result3d, labels, standards, Title)
                % draw picture according to result and labels
                f = figure;
                notempty = size(labels, 1);
                Title = strrep(Title, '_', '\_');
                if notempty
                    standardcolors = analyst.getcolors(max(labels), 'jet');
                    for i = 1 : max(labels)
                        toplot = result3d(labels == i, :);
                        [n, ~] = size(toplot);
                        color = ones(n, 1) * standardcolors(i, :);
                        scatter3(toplot(:, 1), toplot(:, 2), toplot(:, 3), 20, color, 'filled');
                        hold on;
                    end
                    hold off;
                    legend(strrep(standards, '_', '\_'), 'Location', 'eastoutside');
                    title(Title);
                else
                    title(Title);
                end
            end
        end
        
        %% functions about counts
        function f = countsshow(obj, fontsize, filename)
            % show counts in a picture
            if ~exist('fontsize', 'var')
                fontsize = 12;
            end
            if ~exist('filename', 'var')
                % no filename, don't save
                save = 0;
            else
                save = 1;
            end
            if ~size(obj.counts, 1)
                disp('counts not found in object, now it will run "getconnection" first');
                obj = obj.getconnection(3);
            end
            alphadata = ones(size(obj.connection, 1), size(obj.connection, 2));
            for i = 1 : size(obj.connection, 1)
                for j = 1 : size(obj.connection, 2)
                    if j > i
                        alphadata(i, j) = 0;
                    end
                end
            end
            f = figure;
            imagesc(obj.connection, 'AlphaData', alphadata);
            colorbar('southoutside');
            axis off;
            set(f,'units','normalized','position',[0.1 0.1 0.8 0.8]);
            for i = 1 : size(obj.standards, 1)
                position1 = i - 0.5;
                position = i;
                line([0.5, position+0.5], [position1, position1], 'Color', 'white', 'LineWidth', 0.05)
                text(0.5, position, strrep([num2str(i), ':', obj.standards{i}, '-'], '_', '\_'), 'Color', 'red', 'FontSize', fontsize, 'HorizontalAlignment', 'right');
                line([position+0.5, position+0.5], [position1, size(obj.standards, 1)+0.5], 'Color', 'white', 'LineWidth', 0.05)
                text(position, position - 0.5, strrep(['-', num2str(i), ':', obj.standards{i}], '_', '\_'), 'Color', 'red', 'FontSize', fontsize, 'HorizontalAlignment', 'left', 'Rotation', 20);
                for j = 1 : i
                    count = obj.counts{i, j};
                    countnumber = size(count, 1);
                    if countnumber
                        if j ~= i
                            countfromA = size(unique(count(:, 1)), 1);
                            countfromB = size(unique(count(:, 2)), 1);
                        else
                            countfromA = size(unique([count(:, 1);count(:, 2)]), 1);
                            countfromB = countfromA;
                        end
                    else
                        countfromA = 0;
                        countfromB = 0;
                    end
                    text(j, i, sprintf([num2str(countnumber), '\n', num2str(countfromA), '/', num2str(countfromB), '\n', num2str(obj.clustercounts(i)), '/', num2str(obj.clustercounts(j))]), 'Color', 'white', 'FontSize', fontsize, 'HorizontalAlignment', 'center');
                end
            end
            if save
                saveas(f, [obj.outputpath, filename, '.jpg']);
            end
        end
        
        function ignored = writecounts(obj, filename) %#ok<STOUT>
            % write count number into 3 txt files
            % filename_connections.txt contains the exact number of connections
            % filename_cellcounts.txt counts the number of unique cells from one cluster
            % filename_totalcells.txt contains the exact number of cells in each cluster
            file = fopen([obj.outputpath, filename, '_connections.txt'], 'w');
            fprintf(file, '%s\t', 'number');
            for i = 1 : size(obj.standards, 1)
                fprintf(file, '%s\t', obj.standards{i});
            end
            fprintf(file, '\n');
            for i = 1 : size(obj.counts, 1)
                fprintf(file, '%s\t', obj.standards{i});
                for j = 1 : size(obj.counts, 2)
                    fprintf(file, '%d\t', size(obj.counts{i, j}, 1));
                end
                fprintf(file, '\n');
            end
            fclose(file);
            %write countfromA and countfromB into txt file
            file = fopen([obj.outputpath, filename, '_cellcounts.txt'], 'w');
            fprintf(file, '%s\t', 'clusterA/clusterB');
            for i = 1 : size(obj.standards, 1)
                fprintf(file, '%s\t', obj.standards{i});
            end
            fprintf(file, '\n');
            for i = 1 : size(obj.counts, 1)
                fprintf(file, '%s\t', obj.standards{i});
                for j = 1 : size(obj.counts, 2)
                    if size(obj.counts{i, j}, 1)
                        if i == j
                            fprintf(file, '%d\t', size(unique([obj.counts{i, j}(:, 1); obj.counts{i, j}(:, 2)]), 1));
                        else
                        fprintf(file, '%d\t', size(unique(obj.counts{i, j}(:, 1)), 1));
                        end
                    end
                end
                fprintf(file, '\n');
            end
            fclose(file);
            %write cluster count from A and cluster count from B into txt file
            file = fopen([obj.outputpath, filename, '_totalcells.txt'], 'w');
            fprintf(file, '%s\t', 'clusterA/clusterB');
            for i = 1 : size(obj.standards, 1)
                fprintf(file, '%s\t', obj.standards{i});
            end
            fprintf(file, '\n');
            for i = 1 : size(obj.counts, 1)
                fprintf(file, '%s\t', obj.standards{i});
                for j = 1 : size(obj.counts, 2)
                    fprintf(file, '%d\t', obj.clustercounts(i));
                end
                fprintf(file, '\n');
            end
            fclose(file);
        end
        
        %% functions about degree
        function ignored = writedegree(obj, filename) %#ok<STOUT>
            % write degree into a txt file
            file = fopen([obj.outputpath, filename, '.txt'], 'w');
            fprintf(file, '%s\t', 'degree');
            for i = 1 : size(obj.degree, 2)
                fprintf(file, '%s\t', obj.standards{i});
            end
            fprintf(file, '\n');
            for i = 1 : size(obj.degree, 1)
                fprintf(file, '%s\t', obj.standards{i});
                for j = 1 : size(obj.connection, 2)
                    fprintf(file, '%d\t', obj.degree(i, j));
                end
                fprintf(file, '\n');
            end
            fclose(file);
        end
             
        %% function about differential genes
        function [genes, genes_foldchange, ABcellsTPM, AnotBcellsTPM] = differential_genes(obj, clusterA, clusterB, pcutoff, fcutoff, testtype, draw, filename)
            % find the differentially expressed genes bewteen cells in clusterA that
            % interacts with clusterB and those cells in clusterA that
            % doesn't, genes_foldchange are the cooresponding genes' fold
            % change, ABcellsTPM is the TPM of all the cells in clusterA
            % that interacts with clusterB, AnotBcellsTPM is the TOM of all
            % the cells in clusterA that doesn't interact with clusterB.
            close all;
            if ~exist('testtype', 'var')
                testtype = 'ttest';
            end
            if ~exist('draw', 'var')
                draw = 0;
                save = 0;
            end
            if ~exist('filename', 'var')
                save = 0;
            end
            if isa(clusterA, 'char')
                clusterA = find(strcmp(obj.standards, clusterA));
            end
            if isa(clusterB, 'char')
                clusterB = find(strcmp(obj.standards, clusterB));
            end
            allindexes = (1:size(obj.labels, 1))';
            allAcells = allindexes(obj.labels == clusterA);
            if clusterA ~= clusterB
                ABcells = unique(obj.counts{clusterA, clusterB}(:, 1));
            else
                ABcells = unique([obj.counts{clusterA, clusterB}(:, 1); obj.counts{clusterA, clusterB}(:, 2)]);
            end
            ABcellsTPM = obj.TPM(:, ABcells);
            AnotBcells = allAcells(~ismember(allAcells, ABcells));
            AnotBcellsTPM = obj.TPM(:, AnotBcells);
            if strcmp(testtype, 'ttest')
                p = mattest(log2(ABcellsTPM+1), log2(AnotBcellsTPM+1));
            else
                p = zeros(size(ABcellsTPM, 1), 1);
                for i = 1 : size(ABcellsTPM, 1)
                    p(i, 1) = ranksum(log2(ABcellsTPM(i,:)+1), log2(AnotBcellsTPM(i,:)+1));
                end
            end
            [~, p] = mafdr(p);
            fold_change = mean(ABcellsTPM, 2) ./ mean(AnotBcellsTPM, 2);
            if fcutoff > 0
                genes = obj.genes(p<pcutoff&fold_change>fcutoff);
                genes_foldchange=fold_change(p<pcutoff&fold_change>fcutoff);
            elseif fcutoff < 0
                genes = obj.genes(p<pcutoff&fold_change<-1/fcutoff);
                genes_foldchange=fold_change(p<pcutoff&fold_change<-1/fcutoff);
                fcutoff=-fcutoff;
            else
                genes = obj.genes(p<pcutoff);
                genes_foldchange=fold_change(p<pcutoff);
            end
            [~, geneindex] = ismember(genes, obj.genes);
            allindexes = 1:size(obj.genes, 1);
            otherindex = allindexes(~ismember(allindexes, geneindex));
            if draw
                heatmat = log2(obj.TPM(geneindex, [ABcells; AnotBcells])/10+1);
                heatmat = zscore(heatmat, 0, 2);
                f = figure;        
                imagesc(heatmat);
                colormap parula;
                colorbar;
                line([size(ABcells, 1), size(ABcells, 1)], [0, size(allAcells, 1)], 'Color', 'red', 'LineWidth', 0.05)
                set(gca, 'ytick', 1:size(genes, 1), 'yticklabel', genes);
                g = figure;
                scatter(log2(fold_change(otherindex)), -log10(p(otherindex)), 'filled');
                xlabel('log_2 (fold\_change)');
                ylabel('-log_{10} (p)');
                hold on;
                scatter(log2(fold_change(geneindex)), -log10(p(geneindex)), 'filled', 'r');
                line([-10,20], [-log10(pcutoff), -log10(pcutoff)], 'Color', 'red', 'LineStyle', ':');
                hold on;
                line([log2(fcutoff), log2(fcutoff)], [0,25], 'Color', 'red', 'LineStyle', ':');
                hold on;
                line([-log2(fcutoff), -log2(fcutoff)], [0,25], 'Color', 'red', 'LineStyle', ':');
                for i = 1 : size(geneindex)
                    text(log2(fold_change(geneindex(i)))+0.2, -log10(p(geneindex(i))), genes(i));
                end
                axis tight;
                if save
                    print(f, [obj.outputpath, filename, '_heatmap.pdf'], '-dpdf', '-bestfit');
                    print(g, [obj.outputpath, filename, '_vocano.pdf'], '-dpdf', '-bestfit');
                end
            end
        end
        
        %% functions about statistical results and reverse statis(p-value matrix)
        function f = statisticsshow(obj, fontsize, filename)
            % show the q-value matrix as a picture
            if ~exist('fontsize', 'var')
                fontsize = 12;
            end
            if ~exist('filename', 'var')
                % no filename, don't save
                save = 0;
            else
                save = 1;
            end
            alphadata = ones(size(obj.connection, 1), size(obj.connection, 2));
            for i = 1 : size(obj.connection, 1)
                for j = 1 : size(obj.connection, 2)
                    if j > i
                        alphadata(i, j) = 0;
                    end
                end
            end
            f = figure;
            imagesc(obj.connection, 'AlphaData', alphadata);
            colorbar('southoutside');
            axis off;
            set(f,'units','normalized','position',[0.1 0.1 0.8 0.8]);
            for i = 1 : size(obj.standards, 1)
                position1 = i - 0.5;
                position = i;
                line([0.5, position+0.5], [position1, position1], 'Color', 'white', 'LineWidth', 0.05)
                text(0.5, position, strrep([num2str(i), ':', obj.standards{i}, '-'], '_', '\_'), 'Color', 'red', 'FontSize', fontsize, 'HorizontalAlignment', 'right');
                line([position+0.5, position+0.5], [position1, size(obj.standards, 1)+0.5], 'Color', 'white', 'LineWidth', 0.05)
                text(position, position - 0.5, strrep(['-', num2str(i), ':', obj.standards{i}], '_', '\_'), 'Color', 'red', 'FontSize', fontsize, 'HorizontalAlignment', 'left', 'Rotation', 20);
                for j = 1 : i
                    text(j, i, num2str(obj.connection(i, j), '%.4f'), 'Color', 'white', 'FontSize', fontsize, 'HorizontalAlignment', 'center');
                end
            end
            if save
                saveas(f, [obj.outputpath, filename, '.jpg']);
            end
        end
        
        function f = reversestatisticsshow(obj, fontsize, filename)
            % show reverse q-value as a picture
            if ~exist('fontsize', 'var')
                fontsize = 12;
            end
            if ~exist('filename', 'var')
                % no filename, don't save
                save = 0;
            else
                save = 1;
            end
            alphadata = ones(size(obj.connection, 1), size(obj.connection, 2));
            for i = 1 : size(obj.reverseconnection, 1)
                for j = 1 : size(obj.reverseconnection, 2)
                    if j > i
                        alphadata(i, j) = 0;
                    end
                end
            end
            f = figure;
            imagesc(obj.reverseconnection, 'AlphaData', alphadata);
            colorbar('southoutside');
            axis off;
            set(f,'units','normalized','position',[0.1 0.1 0.8 0.8]);
            for i = 1 : size(obj.standards, 1)
                position1 = i - 0.5;
                position = i;
                line([0.5, position+0.5], [position1, position1], 'Color', 'white', 'LineWidth', 0.05)
                text(0.5, position, strrep([num2str(i), ':', obj.standards{i}, '-'], '_', '\_'), 'Color', 'red', 'FontSize', fontsize, 'HorizontalAlignment', 'right');
                line([position+0.5, position+0.5], [position1, size(obj.standards, 1)+0.5], 'Color', 'white', 'LineWidth', 0.05)
                text(position, position - 0.5, strrep(['-', num2str(i), ':', obj.standards{i}], '_', '\_'), 'Color', 'red', 'FontSize', fontsize, 'HorizontalAlignment', 'left', 'Rotation', 20);
                for j = 1 : i
                    text(j, i, num2str(obj.reverseconnection(i, j), '%.4f'), 'Color', 'white', 'FontSize', fontsize, 'HorizontalAlignment', 'center');
                end
            end
            if save
                saveas(f, [obj.outputpath, filename, '.jpg']);
            end
        end
        
        function ignored = writestatistics(obj, filename) %#ok<STOUT>
            %write q-value matrix into txt file
            file = fopen([obj.outputpath, filename, '.txt'], 'w');
            fprintf(file, '%s\t', 'p-value');
            for i = 1 : size(obj.standards, 1)
                fprintf(file, '%s\t', obj.standards{i});
            end
            fprintf(file, '\n');
            for i = 1 : size(obj.connection, 1)
                fprintf(file, '%s\t', obj.standards{i});
                for j = 1 : size(obj.connection, 2)
                    fprintf(file, '%.4f\t', obj.connection(i, j));
                end
                fprintf(file, '\n');
            end
            fclose(file);
        end
        
        function conclusion = drawconclusion(obj, pcutoff, filename)
            % Save enriched interactions to a txt file, based on q-value and cutoff.
            conclusion = {};
            for i = 1 : size(obj.connection, 1)
                for j = 1 : i
                    if obj.connection(i, j) <= pcutoff
                        a = [cell2mat(obj.standards(i)), '---', cell2mat(obj.standards(j))];
                        conclusion = [conclusion; a]; %#ok<AGROW>
                    end
                end
            end
            conclusion = sort(conclusion);
            file = fopen([obj.outputpath, filename, '.txt'], 'w');
            for i = 1 : size(conclusion, 1)
                fprintf(file, '%s\n', conclusion{i});
            end
            fclose(file);
        end
        
        %% functions about presenting the cells' coordinates
        function f = scatter2d(obj, Title, filename)
            % draw picture according to result and labels
            % ATTENTION, obj.result2d must not be empty
            if ~exist('Title', 'var')
                Title = 'scatter 3d';
            end
            if ~exist('filename', 'var')
                % no filename, don't save
                save = 0;
            else
                save = 1;
            end
            Title = strrep(Title, '_', '\_');
            if ~size(obj.result2d, 1)
                Presult = tsne(obj.result3d, [], 2);
            else
                Presult = obj.result2d;
            end
            f = figure;
            notempty = size(obj.labels, 1);
            if notempty
                standardcolors = analyst.getcolors(max(obj.labels), 'jet');
                for i = 1 : max(obj.labels)
                    toplot = Presult(obj.labels == i, :);
                    [n, ~] = size(toplot);
                    color = ones(n, 1) * standardcolors(i, :);
                    scatter(toplot(:, 1), toplot(:, 2), 20, color, 'filled');
                    hold on;
                end
                hold off;
                legend(strrep(obj.standards, '_', '\_'), 'Location', 'eastoutside');
                title(Title);
            else
                title(Title);
            end
            set(f,'units','normalized','position',[0.1 0.1 0.8 0.8]);
            if save
                saveas(f, [obj.outputpath, filename, '.jpg']);
            end
        end
        
        function f = scatter3d(obj, Title, filename, iter)
            % draw picture according to result and labels
            if ~exist('iter', 'var')
                coordinates = obj.result3d;
            else
                coordinates = obj.process(:,iter:iter+2);
            end
            if ~exist('Title', 'var')
                Title = 'scatter 3d';
            end
            if ~exist('filename', 'var')
                % no filename, don't save
                save = 0;
            else
                save = 1;
            end
            f = figure;
            notempty = size(obj.labels, 1);
            Title = strrep(Title, '_', '\_');
            if notempty
                standardcolors = analyst.getcolors(max(obj.labels), 'parula');
                for i = 1 : max(obj.labels)
                    toplot = coordinates(obj.labels == i, :);
                    [n, ~] = size(toplot);
                    color = ones(n, 1) * standardcolors(i, :);
                    scatter3(toplot(:, 1), toplot(:, 2), toplot(:, 3), 50, color, 'filled');
                    hold on;
                end
                hold off;
                legend(strrep(obj.standards, '_', '\_'), 'Location', 'eastoutside');
                title(Title);
            else
                title(Title);
            end
            set(f,'units','normalized','position',[0.1 0.2 0.4 0.6]);
            axis tight;
            xlabel('x');
            ylabel('y');
            zlabel('z');
            if save
                print(f, [obj.outputpath, filename, '.pdf'], '-dpdf', '-bestfit');
                saveas(f, [obj.outputpath, filename, '.fig']);
            end
            xticks(-60:20:60);
            yticks(-60:20:60);
            zticks(-60:20:60);
        end
        
        function f = getsections(obj, Axis, Title, filename)
            % draw several sections from z axis's view
            if ~exist('filename', 'var')
                save = 0;
            else
                save = 1;
            end
            if ~exist('Title', 'var')
                Title = [Axis, ' view'];
            end
            f = obj.scatter3d();
            for i = -50:20:30
                if strcmp(Axis, 'x')
                    % view from x axis
                    axis([i, i+20, -inf, inf, -inf, inf]);
                    view([90, 0, 0]);
                    title([Title, ' range: [', num2str(i), ', ', num2str(i+20), ']']);
                elseif strcmp(Axis, 'y')
                    % view from y axis
                    axis([-inf, inf, i, i+20, -inf, inf]);
                    view([0, 90, 0]);
                    title([Title, ' range: [', num2str(i), ', ', num2str(i+20), ']']);
                elseif strcmp(Axis, 'z')
                    % view from z axis
                    axis([-inf, inf, -inf, inf, i, i+20]);
                    view([0, 0, 90]);
                    title([Title, ' range: [', num2str(i), ', ', num2str(i+20), ']']);
                end
                if save
                    print(f, [obj.outputpath, filename, Title, ' range[', num2str(i), ',', num2str(i+20), '].pdf'], '-dpdf', '-bestfit');
                end
            end
            %back to origin
            if save
                axis([-inf, inf, -inf, inf, -inf, inf]);
                view([90, 90, 90]);
                close all;
            end
        end
        
        function f = scattercluster(obj, cluster, Title, filename)
            % scatter only one cluster
            if ~exist('filename', 'var')
                save = 0;
            else
                save = 1;
            end
            if ~exist('Title', 'var')
                Title = ['scatter of ', cluster];
            end
            f = figure;
            notempty = size(obj.labels, 1);
            Title = strrep(Title, '_', '\_');
            i = find(strcmp(obj.standards, cluster));
            if notempty
                standardcolors = analyst.getcolors(max(obj.labels), 'jet');
                toplot = obj.result3d(obj.labels == i, :);
                [n, ~] = size(toplot);
                color = ones(n, 1) * standardcolors(i, :);
                scatter3(toplot(:, 1), toplot(:, 2), toplot(:, 3), 20, color, 'filled');
                legend(strrep(obj.standards{i}, '_', '\_'), 'Location', 'eastoutside');
                title(Title);
            else
                title(Title);
            end
            set(f,'units','normalized','position',[0.1 0.1 0.8 0.8]);
            xlabel('x');
            ylabel('y');
            zlabel('z');
            if save
                print(f, [obj.outputpath, filename, '.pdf'], 'dpdf', '-bestfit');
            end
        end
        
        function ignored = savegif(obj, Title, filename) %#ok<STOUT>
            % view 3d plot as gif
            picture = obj.scatter3d(Title, filename);
            axis vis3d
            shading interp
            mkdir([obj.outputpath, 'rotate/']);
            for i=1:36
                camorbit(10,0,'data',[0 0 1])
                print(picture, [obj.outputpath, 'rotate/', filename, '_', num2str(i*10), '.pdf'], '-dpdf', '-bestfit');
                M=getframe(picture);
                nn=frame2im(M);
                [nn,cm]=rgb2ind(nn,256);
                if i==1
                    imwrite(nn,cm,[obj.outputpath, filename, '.gif'],'gif','LoopCount',inf,'DelayTime',0.25);%loopcount i==1
                    else
                    imwrite(nn,cm,[obj.outputpath, filename, '.gif'],'gif','WriteMode','append','DelayTime',0.25)%i>=2 loopcount
                end
            end
        end
        
        %% function about a certain gene's expression
        function f = expressionshow(obj, gene, filename, Title)
            % show one gene's expression pattern
            if ~exist('Title', 'var')
                Title = ['Expression of ', gene];
            end
            if ~exist('filename', 'var')
                % no filename, don't save
                save = 0;
            else
                save = 1;
            end
            geneTPM = log2(obj.TPM(strcmp(obj.genes, gene), :)'+1);
            f = figure;
            scatter3(obj.result3d(:, 1), obj.result3d(:, 2), obj.result3d(:, 3), 'filled', 'cdata', geneTPM);
            title(Title);
            colorbar;
            if save
                saveas(f, [obj.outputpath, filename]);
            end
        end
        
        %% functions about main contributed Ligand and Receptor
        function ignored = mainLR(obj, filename) %#ok<STOUT>
            % calculates the contribution of each pair of Ligand and Receptor
            % save it to a txt file.
            file = fopen([obj.outputpath, filename, '.txt'], 'w');
            if ~size(obj.counts, 1)
                disp('counts not found in object, now it will run "getconnection" first');
                obj = obj.getconnection(3);
            end
            Atotake = obj.ligandindex;
            Btotake = obj.receptorindex;
            LRpairs = [obj.ligands, obj.receptors];
            allscores = obj.scores;
            for i = 1 : size(obj.ligandindex, 1)
                if obj.ligandindex(i) ~= obj.receptorindex(i)
                    Atotake = [Atotake; obj.receptorindex(i)];
                    Btotake = [Btotake; obj.ligandindex(i)];
                    LRpairs = [LRpairs; obj.receptors(i), obj.ligands(i)];
                    allscores = [allscores; obj.scores(i)];
                end
            end
            for i = 1 : size(obj.counts, 1)
                for j = 1 : i
                    count = obj.counts{i, j};
                    allconclusion = cell(0, 2);
                    if size(count, 1)
                        fprintf(file, '%s\t%s\n\n', [obj.standards{i}, '---', obj.standards{j}], ['q-value: ', num2str(obj.connection(i, j))]);
                        fprintf(file, '%s\t%s\n', 'Ligand---Receptor', 'contribution');
                        for k = 1 : size(count, 1)
                            A = obj.TPM(Atotake, count(k, 1));
                            B = obj.TPM(Btotake, count(k, 2));
                            affinity = sum(A.*allscores.*B);
                            contributes = A.*allscores.*B ./ affinity;
                            [contributes, I] = sort(contributes, 'descend');
                            top10 = LRpairs(I(1:10), :);
                            contributes = contributes(1:10, :);
                            for l = 1 : size(top10, 1)
                                [~, I] = ismember([top10{l, 1}, '---', top10{l, 2}], allconclusion(:, 1));
                                if I ~= 0
                                    allconclusion{I, 2} = allconclusion{I, 2} + contributes(l);
                                else
                                    if size(allconclusion, 1) == 0
                                        allconclusion = {[top10{l, 1}, '---', top10{l, 2}], contributes(l)};
                                    else
                                        allconclusion = [allconclusion; {[top10{l, 1}, '---', top10{l, 2}], contributes(l)}]; %#ok<AGROW>
                                    end
                                end
                            end
                        end
                    [~, I] = sort(cell2mat(allconclusion(:, 2)), 'descend');
                    allconclusion = allconclusion(I, :);
                    for k = 1 : size(allconclusion, 1)
                        countnumber = size(count, 1);
                        fprintf(file, '%s\t%s\n', allconclusion{k, 1}, num2str(allconclusion{k, 2}/countnumber, '%.6f'));
                    end
                    fprintf(file, '\n');
                    else
                        fprintf(file, '%s\t%s\n\n', [obj.standards{i}, '---', obj.standards{j}], ['p-value: ', num2str(obj.connection(i, j))]);
                    end
                end
            end
            fclose(file);
        end
        
        %% functions about change labels
        function obj = setnewlabels(obj, newlabelpath, restat)
            %change labels, if not restat, some functions cannot be used.
            [newlabels, newstandards] = analyst.identify_label(obj.cells, newlabelpath);
            obj.labels = newlabels;
            obj.standards = newstandards;
            if restat
                obj = obj.getconnection(3);
            else
                % clear old statistical results
                obj.connection = [];
                obj.reverseconnection = [];
                obj.counts = [];
                obj.clustercounts = [];
            end
        end
        
        %% functions about cluster_dp
        function [cl, halo, distribution] = cluster_dp(obj, k, mds, type, show)
            % perform cluster dp, this function is referenced from CDP
            % algorithm, published on Science in 2014.
            % cl and halo are the origin output in that function.
            % distribution calculates each cluster's composition
            close all;
            if ~exist('k', 'var')
                k = 3;
            end
            if ~exist('mds', 'var')
                mds = 0;
            end
            if ~exist('type', 'var')
                type = 'cutoff';
            end
            if ~exist('show','var')
                show = 0;
            end
            xx = obj.getplaindistance(obj.result3d);
            dist = obj.getdistancemat(obj.result3d);
            ND=max(xx(:,2));
            fprintf('median number of neighbours (hard coded): %5.6f\n', k);
            % get cutoff or gaussian kernel
            dist(1:ND+1:end) = Inf;
            cutoffs = zeros(ND, 1);
            for i = 1 : ND
                row = dist(i, :);
                [B, ~] = sort(row);
                cutoffs(i) = B(k);
            end
            dc = median(cutoffs);
            dist(1:ND+1:end) = 0;
%             N=size(xx,1);
%             percent=2.0;
%             fprintf('average percent of neighbours (hard coded): %5.6f\n', percent);
%             position=round(N*percent/100);
%             sda=sort(xx(:,3));
%             dc=sda(position);            
            fprintf('Computing Rho with gaussian kernel of radius: %12.6f\n', dc);
            
            rho = zeros(ND, 1);
            
            if strcmp(type, 'gaussian')
                % Gaussian kernel
                for i=1:ND-1
                    for j=i+1:ND
                        rho(i)=rho(i)+exp(-(dist(i,j)/dc)*(dist(i,j)/dc));
                        rho(j)=rho(j)+exp(-(dist(i,j)/dc)*(dist(i,j)/dc));
                    end
                end
            elseif strcmp(type, 'cutoff')
                % "Cut off" kernel
                for i=1:ND-1
                    for j=i+1:ND
                        if (dist(i,j)<dc)
                            rho(i)=rho(i)+1.;
                            rho(j)=rho(j)+1.;
                        end
                    end
                end
            end

            maxd=max(max(dist));

            [~,ordrho]=sort(rho,'descend');
            delta(ordrho(1))=-1.;
            nneigh(ordrho(1))=0;

            for ii=2:ND
               delta(ordrho(ii))=maxd;
               for jj=1:ii-1
                 if(dist(ordrho(ii),ordrho(jj))<delta(ordrho(ii)))
                    delta(ordrho(ii))=dist(ordrho(ii),ordrho(jj));
                    nneigh(ordrho(ii))=ordrho(jj);
                 end
               end
            end
            delta(ordrho(1))=max(delta(:));
%             disp('Generated file:DECISION GRAPH')
%             disp('column 1:Density')
%             disp('column 2:Delta')
% 
%             fid = fopen('DECISION_GRAPH', 'w');
%             for i=1:ND
%                 fprintf(fid, '%6.2f %6.2f\n', rho(i),delta(i));
%             end
%             disp('Select a rectangle enclosing cluster centers')
            scrsz = get(0,'ScreenSize');
            if mds
                figure('Position',[6 12 scrsz(3)/4. scrsz(4)/1.3]);
            else
                figure('Position',[6 12 scrsz(3)/2.5 scrsz(4)/2.5]);
            end
            for i=1:ND
                ind(i)=i;
                gamma(i)=rho(i)*delta(i);
            end
            if mds
                subplot(2,1,1)
            end
            tt=plot(rho(:),delta(:),'o','MarkerSize',5,'MarkerFaceColor','k','MarkerEdgeColor','k');
            title ('Decision Graph','FontSize',15.0)
            xlabel ('\rho')
            ylabel ('\delta')

            if mds
                subplot(2,1,1)
            end
            if show
                disp('Please Enter min rho and min delta')
                rhomin=input('Min Rho:\n');
                deltamin=input('Min Delta:\n');
            else
                rhomin=prctile(rho, 80);
                deltamin=prctile(delta, 80);
            end
            NCLUST=0;
            for i=1:ND
              cl(i)=-1;
            end
            for i=1:ND
              if ( (rho(i)>rhomin) && (delta(i)>deltamin))
                 NCLUST=NCLUST+1;
                 cl(i)=NCLUST;
                 icl(NCLUST)=i;
              end
            end
            fprintf('NUMBER OF CLUSTERS: %i \n', NCLUST);
            disp('Performing assignation')

            %assignation
            for i=1:ND
              if (cl(ordrho(i))==-1)
                cl(ordrho(i))=cl(nneigh(ordrho(i)));
              end
            end
            %halo
            for i=1:ND
              halo(i)=cl(i);
            end
            if (NCLUST>1)
              for i=1:NCLUST
                bord_rho(i)=0.;
              end
              for i=1:ND-1
                for j=i+1:ND
                  if ((cl(i)~=cl(j))&& (dist(i,j)<=dc))
                    rho_aver=(rho(i)+rho(j))/2.;
                    if (rho_aver>bord_rho(cl(i))) 
                      bord_rho(cl(i))=rho_aver;
                    end
                    if (rho_aver>bord_rho(cl(j))) 
                      bord_rho(cl(j))=rho_aver;
                    end
                  end
                end
              end
              for i=1:ND
                if (rho(i)<bord_rho(cl(i)))
                  halo(i)=0;
                end
              end
            end
            for i=1:NCLUST
              nc=0;
              nh=0;
              for j=1:ND
                if (cl(j)==i) 
                  nc=nc+1;
                end
                if (halo(j)==i) 
                  nh=nh+1;
                end
              end
              fprintf('CLUSTER: %i CENTER: %i ELEMENTS: %i CORE: %i HALO: %i \n', i,icl(i),nc,nh,nc-nh);
            end
            cmap=colormap;
            for i=1:NCLUST
               ic=int8((i*64.)/(NCLUST*1.));
               if mds
                subplot(2,1,1)
               end
               hold on
               plot(rho(icl(i)),delta(icl(i)),'o','MarkerSize',8,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
            end
            saveas(gca, [obj.outputpath, 'cluster_dp_rho&delta.jpg']);
            halo = halo';
            cl = cl';
            distribution = cell(NCLUST+2, max(obj.labels)+2);
            distribution(1, 2:max(obj.labels)+1) = obj.standards';
            distribution{1, max(obj.labels)+2} = 'total';
            for i = 2 : NCLUST+2
                distribution{i, 1} = num2str(i-2);
                for j = 2 : max(obj.labels)+1
                    distribution{i, j} = 0;
                end
            end
            total = zeros(NCLUST+1, 1);
            for i = 1:size(halo, 1)
                total(halo(i)+1) = total(halo(i)+1) + 1;
                distribution{halo(i)+2, obj.labels(i)+1} = distribution{halo(i)+2, obj.labels(i)+1} + 1;
            end
            for i = 2 : NCLUST+2
                for j = 2 : max(obj.labels)+1
                    distribution{i, j} = distribution{i, j}/total(i-1);
                    distribution{i, max(obj.labels)+2} = total(i-1);
                end
            end
            if mds
                % plot mds result
                subplot(2,1,2)
                disp('Performing 2D nonclassical multidimensional scaling')
                Y1 = mdscale(dist, 2, 'criterion','metricstress');
                plot(Y1(:,1),Y1(:,2),'o','MarkerSize',2,'MarkerFaceColor','k','MarkerEdgeColor','k');
                title ('2D Nonclassical multidimensional scaling','FontSize',15.0)
                xlabel ('X')
                ylabel ('Y')
                A = zeros(ND, 2);
                for i=1:NCLUST
                  nn=0;
                  ic=int8((i*64.)/(NCLUST*1.));
                  for j=1:ND
                    if (halo(j)==i)
                      nn=nn+1;
                      A(nn,1)=Y1(j,1);
                      A(nn,2)=Y1(j,2);
                    end
                  end
                  hold on
                  plot(A(1:nn,1),A(1:nn,2),'o','MarkerSize',2,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
                end
            else
                % plot 3d result
                f = figure();
                set(f,'units','normalized','position',[0 0 1 0.35]);
                disp('Print 3d plot')
                standardcolors = analyst.getcolors(NCLUST+1, 'jet');
                legends = {};
                for j = 1 : 5
                    subplot(1, 5, j);
                    for i = 0 : NCLUST
                        toplot = obj.result3d(halo == i, :);
                        [n, ~] = size(toplot);
                        color = ones(n, 1) * standardcolors(i+1, :);
                        scatter3(toplot(:, 1), toplot(:, 2), toplot(:, 3), 10, color, 'filled');
                        legends = [legends; num2str(i)];
                        hold on;
                    end
                    axis([-40, 40, -40, 40, -70+j*20+7, -50+j*20-7])
                    view([0, 0, 90]);
                    hold off;
                end
                legend(legends, 'Location', 'eastoutside');
                title('3d plot');
            end
            saveas(gcf, [obj.outputpath, 'cluster_dp_sections.jpg']);
            save([obj.outputpath, 'cluster_dp_distribution.mat'], 'distribution', 'halo', 'cl', '-v7.3');
            %for i=1:ND
            %   if (halo(i)>0)
            %      ic=int8((halo(i)*64.)/(NCLUST*1.));
            %      hold on
            %      plot(Y1(i,1),Y1(i,2),'o','MarkerSize',2,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
            %   end
            %end
        %     faa = fopen('CLUSTER_ASSIGNATION', 'w');
        %     disp('Generated file:CLUSTER_ASSIGNATION')
        %     disp('column 1:element id')
        %     disp('column 2:cluster assignation without halo control')
        %     disp('column 3:cluster assignation with halo control')
        %     for i=1:ND
        %        fprintf(faa, '%i %i %i\n',i,cl(i),halo(i));
        %     end
        end
        
        %% functions about spatial non-randomly distributed genes
        function [p1, p2, p3, p4, p5, p6] = spatial_nonrandom(obj, gene, allrandomcells, dim, expressioncutoff, numhist)
            % check whether a gene is spatially non-randomly distributed
            if ~exist('allrandomcells', 'var')
                n = size(obj.cells, 1);
                allrandomcells = zeros(n, 1000);
                parfor i = 1:1000
                    allrandomcells(:, i) = randperm(n)';
                end
            end
            
            if ~exist('dim', 'var')
                dim = 3;
            end
            if ~exist('expressioncutoff', 'var')
                expressioncutoff = 0;
            end
            if ~exist('numhist', 'var')
                numhist = floor(size(obj.result3d,1) / 10);
            end
            
            if dim ==3
                coordinates = obj.result3d;
            elseif dim == 2
                coordinates = obj.result2d;
            end
            mindist = 0;
            maxdist = sqrt(sum((max(coordinates) - min(coordinates)) .^ 2));
            step = (maxdist - mindist) ./ numhist;
            histrange = mindist+step/2 : step : maxdist-step/2;
            
            indexestotake = obj.TPM(strcmp(obj.genes, gene) , :)' > expressioncutoff;
            realcells = coordinates(indexestotake, :);
            realdist = hist(pdist(realcells), histrange) ./ (size(realcells,1)*(size(realcells,1)-1)/2);
            realmeandist = mean(pdist(realcells));
            randomdist = zeros(1000, size(realdist, 2));
            randommeandist = zeros(1000, 1);
            parfor i = 1 : 1000
                randomcoordinates = coordinates(allrandomcells(:, i), :);
                randomcells = randomcoordinates(indexestotake, :);
                randomdist(i,:) = hist(pdist(randomcells), histrange) ./ (size(realcells,1)*(size(realcells,1)-1)/2);
                randommeandist(i,:) = mean(pdist(randomcells));
            end
            
            meandist = mean(randomdist);
            randombackground = sum(abs(bsxfun(@minus, randomdist, meandist)), 2);
            
            realstat = sum(abs(realdist - meandist));
            
            [muhat, sigmahat] = normfit(randombackground);
            [mu2, sig2] = normfit(randommeandist);
            p1 = normcdf(realstat, muhat, sigmahat, 'upper');
            p2 = sum(randombackground >= realstat) ./ 1000;
            p3 = normcdf(realmeandist, mu2, sig2, 'upper');
            p4 = normcdf(realmeandist, mu2, sig2);
            p5 = sum(randommeandist >= realmeandist) ./ 1000;
            p6 = sum(randommeandist <= realmeandist) ./ 1000;
            disp([gene, ' ', num2str(p1), ' ', num2str(p2), ' ', num2str(p3), ' ' , num2str(p4), ' ', num2str(p5), ' ', num2str(p6)]);
        end
        
        function [p1, p2, p3, p4, p5, p6, allrandomcells, allgenes] = all_spatial_nonrandom(obj)
            n = size(obj.cells, 1);
            allrandomcells = zeros(n, 1000);
            parfor i = 1:1000
                allrandomcells(:, i) = randperm(n)';
            end
            num = size(obj.genes,1);
            allgenes = obj.genes;
            p1 = zeros(num,1);
            p2 = zeros(num,1);
            p3 = zeros(num,1);
            p4 = zeros(num,1);
            p5 = zeros(num,1);
            p6 = zeros(num,1);
            parfor i = 1:num
                try
                    [p_1,p_2,p_3,p_4,p_5,p_6] = obj.spatial_nonrandom(allgenes{i}, allrandomcells);
                    p1(i, 1) = p_1;
                    p2(i, 1) = p_2;
                    p3(i, 1) = p_3;
                    p4(i, 1) = p_4;
                    p5(i, 1) = p_5;
                    p6(i, 1) = p_6;
                catch
                    disp(['error in ', obj.genes{i}]);
                end
            end
            disp('end of this data set');
        end

        %% other supported functions, only for developer.
        function affinitymat = calculate_affinity_mat(obj)
        % This function calculate the affinity mat from ligand and receptor TPM.
            Atotake = obj.ligandindex;
            Btotake = obj.receptorindex;
            allscores = obj.scores;
            for i = 1:size(obj.ligandindex, 1)
                if obj.ligandindex(i) ~= obj.receptorindex(i)
                    Atotake = [Atotake; obj.receptorindex(i)];
                    Btotake = [Btotake; obj.ligandindex(i)];
                    allscores = [allscores; obj.scores(i)];
                end
            end
            A = obj.TPM(Atotake, :);
            B = obj.TPM(Btotake, :);
            affinitymat = (diag(allscores) * A)' * B;
        end

        function result = discretization(obj, k)
            % this function discretize the affinitymat to denoise
            result = obj.affinitymat;
            n = max(size(obj.affinitymat));
            for i = 1 : n
                row = obj.affinitymat(i, :);
                [~, I] = sort(row, 'descend');
                for j = k+1 : n
                    result(i, I(j)) = 0;
                end
            end
            result = (result + result') / 2;
        end
        
        %% manually set the attributes
        function obj = set.TPM(obj, newTPM)
	    obj.TPM = newTPM;
	end

        function obj = set.labels(obj, newlabels)
            obj.labels = newlabels;
        end
        
        function obj = set.standards(obj, newstandards)
            obj.standards = newstandards;
        end
        
        function obj = set.ligandindex(obj, newligandindex)
            obj.ligandindex = newligandindex;
        end
        
        function obj = set.receptorindex(obj, newreceptorindex)
            obj.receptorindex = newreceptorindex;
        end
        
        function obj = set.result2d(obj, newresult2d)
            obj.result2d = newresult2d;
        end
        
        function obj = set.result3d(obj, newresult3d)
            obj.result3d = newresult3d;
        end
        
        function obj = set.process(obj, newprocess)
            obj.process = newprocess;
        end
        
        function obj = set.outputpath(obj, newpath)
            obj.outputpath = newpath;
        end
        
        function obj = set.connection(obj, newconnection)
            obj.connection = newconnection;
        end
        
        function obj = set.counts(obj, newcounts)
            obj.counts = newcounts;
        end
        
        function obj = set.clustercounts(obj, newclustercounts)
            obj.clustercounts = newclustercounts;
        end
        
        function obj = set.affinitymat(obj, newaffinitymat)
            obj.affinitymat = newaffinitymat;
        end

    end
    methods (Static)
        function [ilabels, standards] = identify_label(cells, labelpath)
        % Get the label of cells from labelset, return labels and standards.
        if exist(labelpath, 'dir')
            labelfile = fopen([labelpath, 'label.txt']);
        elseif exist(labelpath, 'file')
            labelfile = fopen(labelpath);
        else
            disp('Cannot find label, check your labelpath');
        end
        labelset = textscan(labelfile, '%s%s');
        fclose(labelfile);
        [~, cellindex] = ismember(cells, labelset{1});
        slabels = cell(size(cells, 1), 1); %labels, string
        for i = 1 : size(cells, 1)
            if cellindex(i)
                slabels{i, 1} = labelset{2}{cellindex(i), 1};
            else
                slabels{i, 1} = 'unlabeled';
            end
        end
        [standards, ~, ilabels] = unique(slabels); % labels, number and standards
        end
        
        function colors = getcolors(x, name) %#ok<STOUT>
            % return the standard color R, G, B, according to number x
            command = ['colors = ', name, '(', num2str(x), ');'];
            eval(command);
        end
        
        function distancemat = getdistancemat(tsneresult)
            % calculate the distances between all cells
            n = size(tsneresult, 1);
            sum_tsneresult = sum(tsneresult .^ 2, 2);
            distancemat = bsxfun(@plus, sum_tsneresult, bsxfun(@plus, sum_tsneresult', -2 * (tsneresult * tsneresult')));
            distancemat = sqrt(distancemat);
            distancemat(1:n+1:end) = Inf;
        end
        
        function plaindistance = getplaindistance(tsneresult)
            % same as distancemat, only represent in an array-like way.
            n = size(tsneresult, 1);
            sum_tsneresult = sum(tsneresult .^ 2, 2);
            distancemat = bsxfun(@plus, sum_tsneresult, bsxfun(@plus, sum_tsneresult', -2 * (tsneresult * tsneresult')));
            distancemat = sqrt(distancemat);
            plaindistance = zeros(n*(n-1), 3);
            k = 1;
            for i = 1 : n
                for j = 1 : n
                    if j ~= i
                        plaindistance(k, 1) = i;
                        plaindistance(k, 2) = j;
                        plaindistance(k, 3) = distancemat(i, j);
                        k = k + 1;
                    end
                end
            end
        end
   end
end
