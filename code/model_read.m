function [numVertices, numFaces, flist, vlist] = model_read(model)

[~, ~, fileExt] = fileparts(model);

fileExt = lower(fileExt);

if strcmp(fileExt, '.ply')
    [v, f] = read_ply(model);
    vlist = v';
    flist = f';
    numVertices = length(vlist(1, :));
    numFaces = length(flist(1, :));

elseif strcmp(fileExt, '.off')
    fid = fopen(model, 'r');
    if fid == -1
        error('无法打开文件：%s', model);
    end

    try
        [numVertices, numFaces, flist, vlist] = offopen(fid);
    catch ME
        fclose(fid);
        rethrow(ME);
    end

    fclose(fid);
else
    error('不支持的文件格式：%s', fileExt);
end

vlist(1, :) = vlist(1, :) - min(vlist(1, :));
vlist(2, :) = vlist(2, :) - min(vlist(2, :));
vlist(3, :) = vlist(3, :) - min(vlist(3, :));
max_vertex = max(max(vlist));
while max_vertex >= 1
    vlist = vlist/10;
    max_vertex = max(max(vlist));
end
end
