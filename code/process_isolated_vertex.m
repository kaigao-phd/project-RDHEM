function [pred_vertex, used_vertex, isolated_vertices] = process_isolated_vertex(sorted_num_neighbor_element, I, num_neighbor_element, pred_vertex, used_vertex, numVertices, vlist)

zero_indices = find(sorted_num_neighbor_element == 0);
isolated_vertices = I(zero_indices);
for i = 1:length(isolated_vertices)
    iso_pred = isolated_vertices(i);
    search_radius = 1;

    while num_neighbor_element(iso_pred) == 0
        % 检查当前半径的左右两个位置
        left_idx = iso_pred - search_radius;
        right_idx = iso_pred + search_radius;

        % 检查右侧点
        if right_idx <= numVertices && num_neighbor_element(right_idx) ~= 0
            iso_pred = right_idx;
            pred_vertex(:, isolated_vertices(i)) = vlist(:, iso_pred);
            used_vertex(isolated_vertices(i)) = 1;
            break;
        end

        % 检查左侧点
        if left_idx >= 1 && num_neighbor_element(left_idx) ~= 0
            iso_pred = left_idx;
            pred_vertex(:, isolated_vertices(i)) = vlist(:, iso_pred);
            used_vertex(isolated_vertices(i)) = 1;
            break;
        end

        search_radius = search_radius + 1;

    end

end
end