function [neighbor, num_neighbor_element]=get_neighbors(numFaces, flist, numVertices)
neighbor = zeros(200, numVertices);
num_neighbor_element = zeros(1, numVertices);

for v = 1: numVertices
    index = 1;
    for f = 1:numFaces

        if (flist(1, f) == v || flist(2, f) == v || flist(3, f) == v)
            if flist(1, f) == v
                neighbor(index, v) = flist(2, f);
                neighbor(index+1, v) = flist(3, f);
                index = index + 2;
            elseif flist(2, f) == v
                neighbor(index, v) = flist(1, f);
                neighbor(index+1, v) = flist(3, f);
                index = index + 2;
            else
                neighbor(index, v) = flist(1, f);
                neighbor(index+1, v) = flist(2, f);
                index = index + 2;
            end

        end
    end
    neighbors = neighbor(:, v);
    neighbors(neighbors == 0)=[];
    neighbors = unique(neighbors);
    neighbor(:, v) = 0;
    neighbor(1:length(neighbors), v) = neighbors;
    num_neighbor_element(v) = nnz(neighbors);
end
end

