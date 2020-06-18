function mappedX = my_umap(X)
    % Here X is affinity_mat
    X = double(X);
    % Write affinity matrix, added by Guojie Zhong
    h = fopen('affinitymat.dat', 'wb');
    fwrite(h, X, 'double');
    fclose(h);
    disp('Data written');
    my_umap_path='./';
    cmd = sprintf('python %s %d',fullfile(my_umap_path,'/my_umap.py'), size(X, 1));
    [flag, cmdout] = system(cmd, '-echo');
    if(flag~=0)
        error(cmdout);
    end
    mappedX = read_data('result.dat', size(X, 1), 3);   
    delete('affinitymat.dat');
    delete('result.dat');
end

% Reads the result file from umap implementation
function X = read_data(file_name, cellnumber, dimension)
    h = fopen(file_name, 'rb');
    X = fread(h, cellnumber * dimension, 'double');
    X = reshape(X, [dimension cellnumber])';
    fclose(h);
end
