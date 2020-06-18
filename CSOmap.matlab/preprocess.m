function ignored = preprocess(rawdatapath, processed_datapath, pseudo, LRpairstoadd, genestochange, TPMtochange, newcells, newcellsTPM, cellstoknock) %#ok<STOUT>
    % This function import raw data and process them into a *.mat file.
    % Also, you can use it in pseudo mode, which is provided in
    % pseudo.m file.
    %% set variables properly
    if ~exist('pseudo', 'var') || isempty(pseudo)
        pseudo = 0;
        LRpairstoadd = cell(0, 3);
        genestochange = cell(0, 1);
    end
    if ~exist('LRpairstoadd', 'var') || isempty(LRpairstoadd)
        LRpairstoadd = cell(0, 3);
    end
    if ~exist('genestochange', 'var') || ~exist('TPMtochange', 'var') || isempty(genestochange) || isempty(TPMtochange)
        genestochange = cell(0, 1);
    end
    if ~exist('newcells', 'var') || ~exist('newcellsTPM', 'var') || isempty(newcells) || isempty(newcellsTPM)
        newcells = cell(0, 1);
    end
    if ~exist('cellstoknock', 'var') || isempty(cellstoknock)
        cellstoknock = cell(0, 1);
    end
    %% import data as cell structure
    file = fopen([rawdatapath, 'LR_pairs.txt']);
    genepairs = textscan(file, '%s%s%s');
    fclose(file);
    data = importdata([rawdatapath, 'TPM.txt']);
    if ~exist(processed_datapath, 'dir')
        mkdir(processed_datapath);
    else
       disp('Warning! directory already exists, this program might change the files in it')
    end
    
    ligands = genepairs{1, 1};
    receptors = genepairs{1, 2};
    scores = str2double(genepairs{1, 3});
    
    %% get the names of ligands, receptors, genes and cells
    % filter ligand-receptor pairs, only keep literature supported L-R pairs
    
    if pseudo
        ligands = [ligands; LRpairstoadd(:, 1)];
        receptors = [receptors; LRpairstoadd(:, 2)];
        scores = [scores; cell2mat(LRpairstoadd(:, 3))];
    end

    % get the TPM data(attention)
    TPM = data.data;

    %get the names and number of all the genes
    genes = data.textdata(2 : end, 1);
    genenumber = size(genes, 1);
    if pseudo && size(genestochange, 1)
        [~, loc] = ismember(genestochange, genes);
        for i = 1:size(loc, 1)
            if loc(i) == 0
                genes = [genes; genestochange(i)];
                genenumber = genenumber + 1;
                TPM = [TPM; TPMtochange(i, :)];
            else
                TPM(loc(i), :) = TPMtochange(i,:);
            end
        end
    end

    %get the names and number of all the cells
    cells = data.textdata(1, 2 : end)';
    cellnumber = size(cells, 1);
    if pseudo && size(newcells, 1)
        cells = [cells; newcells];
        TPM = [TPM, newcellsTPM];
    end
    
    if pseudo && size(cellstoknock, 1)
        cells(cellstoknock, :) = [];
        TPM(:, cellstoknock) = [];
    end

    %calculate the ligands' and receptors' indexes in the TPM matrix
    [~, ligandindex] = ismember(ligands, genes);
    [~, receptorindex] = ismember(receptors, genes);
    found = logical((ligandindex~=0) .* (receptorindex~=0));
    before = max(size(ligandindex));
    ligandindex = ligandindex(found);
    receptorindex = receptorindex(found); %#ok<*NASGU>
    ligands = ligands(found); % names of ligands found in the dataset
    receptors = receptors(found); % names of receptors found in the dataset
    scores = scores(found);
    after = max(size(ligandindex)); % the L-R pairs found in the dataset
    %% save processed data
    clear data;
    save([processed_datapath, 'data'], '-v7.3');
end
