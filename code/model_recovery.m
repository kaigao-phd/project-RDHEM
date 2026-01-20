function [reconstruct_model, reconstruct_length] = model_recovery(flist, numFaces, numVertices, vlist, L, decrypted_model, m,  dict, sum_length_table)

%get neighbor
[neighbor, num_neighbor_element]=get_neighbors( numFaces, flist, numVertices);
[sorted_num_neighbor_element,I] = sort(num_neighbor_element, 'descend');

% Initialize...
% Record reference vertex
Reference = zeros(1, numVertices);
Reference(I(1)) = 1;



reconstruct_length = zeros(3, numVertices);
pointer_bit = L;
pointer_xyz = 1;
pointer_num = 1;
code_length_param = 7;

%extract_huffmandict
dict_r = cell(L+1, 2);


for i = 1: L+1
    dict_r{i, 1} = i;
    current_code_length_bin = [];
    for j = 1: code_length_param
        current_code_length_bin = [current_code_length_bin, bitget(decrypted_model(pointer_xyz, pointer_num), pointer_bit)];
        pointer_bit = pointer_bit - 1;
        if pointer_bit < 1
            pointer_bit = L;
            pointer_xyz = pointer_xyz + 1;
            if pointer_xyz > 3
                pointer_xyz = 1;
                pointer_num = pointer_num + 1;
            end
        end
    end
    current_code_length = bin2dec(num2str(current_code_length_bin));
    code = [];
    for k = 1:current_code_length
        code = [code, bitget(decrypted_model(pointer_xyz, pointer_num), pointer_bit)];
        pointer_bit = pointer_bit - 1;
        if pointer_bit < 1
            pointer_bit = L;
            pointer_xyz = pointer_xyz + 1;
            if pointer_xyz > 3
                pointer_xyz = 1;
                pointer_num = pointer_num + 1;
            end
        end
    end
    dict_r{i, 2} = code;
end
isequal(dict_r, dict)
disp('dict')


%extract label
for i = 1: numVertices
    if Reference(i) ~= 1
        for xyz = 1:3
            codeword = [];
            inner_flag = 0;
            while (inner_flag == 0)
                codeword = [codeword, bitget(decrypted_model(pointer_xyz, pointer_num), pointer_bit)];

                pointer_bit = pointer_bit - 1;
                if pointer_bit < 1
                    pointer_bit = L;
                    pointer_xyz = pointer_xyz + 1;
                    if pointer_xyz > 3
                        pointer_xyz = 1;
                        pointer_num = pointer_num + 1;
                    end
                end

                for k = 1:(L+1)
                    if isequal(dict_r{k, 2}, codeword)
                        recover_length = dict_r{k, 1};
                        reconstruct_length(xyz, i) = recover_length;
                        inner_flag = 1;
                    end

                end
            end

        end
    end

end

isequal(reconstruct_length, sum_length_table)
sum(Reference)
disp('length')

reconstruct_model = zeros(3, numVertices);
valid_res_bit = zeros(3, numVertices);


%recover res_bit

for i = 1:numVertices
    if Reference(i) ~= 1
        for xyz = 1:3
            if reconstruct_length(xyz, i) <= (L-2)
                for bit = L - reconstruct_length(xyz, i) - 1 : -1 : 1
                    res_bit = bitget(decrypted_model(pointer_xyz, pointer_num), pointer_bit);
                    pointer_bit = pointer_bit - 1;
                    if pointer_bit < 1
                        pointer_bit = L;
                        pointer_xyz = pointer_xyz +1;
                        if pointer_xyz > 3
                            pointer_xyz = 1;
                            pointer_num = pointer_num + 1;
                        end
                    end
                    reconstruct_model(xyz, i) = bitset(reconstruct_model(xyz, i), bit, res_bit);
                    %valid
                    valid_res_bit(xyz, i) = bitset(valid_res_bit(xyz, i), bit, bitget(vlist(xyz, i), bit));
                end
            end
        end
    else
        for xyz = 1:3
            for bit = L : -1 : 1
                res_bit = bitget(decrypted_model(pointer_xyz, pointer_num), pointer_bit);
                pointer_bit = pointer_bit - 1;
                if pointer_bit < 1
                    pointer_bit = L;
                    pointer_xyz = pointer_xyz +1;
                    if pointer_xyz > 3
                        pointer_xyz = 1;
                        pointer_num = pointer_num + 1;
                    end
                end
                reconstruct_model(xyz, i) = bitset(reconstruct_model(xyz, i), bit, res_bit);
                %valid
                valid_res_bit(xyz, i) = bitset(valid_res_bit(xyz, i), bit, bitget(vlist(xyz, i), bit));
            end
        end
    end
end

pointer_num - 1
isequal(round(reconstruct_model), round(valid_res_bit))


% extract LM
Recovered_LM = [];
Recovered_LM_index = 1;
while pointer_num < numVertices
    Recovered_LM(Recovered_LM_index) = bitget(decrypted_model(pointer_xyz, pointer_num), pointer_bit);
    Recovered_LM_index = Recovered_LM_index + 1;
    pointer_bit = pointer_bit - 1;
    if pointer_bit < 1
        pointer_bit = L;
        pointer_xyz = pointer_xyz +1;
        if pointer_xyz > 3
            pointer_xyz = 1;
            pointer_num = pointer_num + 1;
        end
    end
end


% recover MSB

recovered_table = zeros(1, numVertices);
recovered_table(I(1)) = 1;
recovered_table(neighbor(1, I(1))) = 1;
LM_index = 1;

for xyz = 1:3
    length_MSB = reconstruct_length(xyz, neighbor(1, I(1)));
    if length_MSB <= (L - 1)
        for bit = L: -1: (L - length_MSB + 1)
            reconstruct_model(xyz, neighbor(1, I(1))) = bitset(reconstruct_model(xyz,neighbor(1, I(1))), bit, bitget(reconstruct_model(xyz, I(1)), bit));
        end
        % diff_bit
        diff_bit = 1 - bitget(reconstruct_model(xyz, I(1)), L-length_MSB);
        reconstruct_model(xyz, neighbor(1, I(1))) = bitset(reconstruct_model(xyz, neighbor(1, I(1))), L-length_MSB, diff_bit);
    else
        for bit = L: -1: 1
            reconstruct_model(xyz, neighbor(1, I(1))) = bitset(reconstruct_model(xyz,neighbor(1, I(1))), bit, bitget(reconstruct_model(xyz, I(1)), bit));
        end
    end
end



rev_pred_vertex = zeros(3, numVertices);

rev_pred_vertex(:, neighbor(1, I(1))) = reconstruct_model(:, I(1));
queue = neighbor(:, I(1));
queue = queue(queue ~= 0);

num_available_neighbor = zeros(1, numVertices);

test_index = 1;
% 检查是否所有邻居都已标记为1
all_vertices = 1:length(recovered_table);
% 找出非孤立点的索引
zero_indices = find(sorted_num_neighbor_element == 0);
isolated_vertices = I(zero_indices);
non_isolated_vertices = setdiff(all_vertices, isolated_vertices);


while ~all(recovered_table(non_isolated_vertices) == 1)
    if isempty(queue)
        % 找到第一个未处理的顶点
        start_vertex = find(recovered_table == 0, 1);
        if ~isempty(start_vertex)
            % 初始化这个新的起始点
            search_radius = 1;

            while recovered_table(start_vertex) == 0
                % 检查当前半径的左右两个位置
                left_idx = start_vertex - search_radius;
                right_idx = start_vertex + search_radius;

                % 检查右侧点
                if right_idx <= numVertices && recovered_table(right_idx) == 1 && ~ismember(right_idx, isolated_vertices)
                    rev_pred_vertex(:, start_vertex) = reconstruct_model(:, right_idx);
                    recovered_table(start_vertex) = 1;
                    break;
                end

                % 检查左侧点
                if left_idx >= 1 && recovered_table(left_idx) == 1 && ~ismember(left_idx, isolated_vertices)
                    rev_pred_vertex(:, start_vertex) = reconstruct_model(:, left_idx);
                    recovered_table(start_vertex) = 1;
                    break;
                end

                search_radius = search_radius + 1;
            end
            
            pred_value = rev_pred_vertex(:, start_vertex);
            for xyz = 1:3
                length_MSB = reconstruct_length(xyz, start_vertex);
                if length_MSB <= (L - 1)
                    for bit = L: -1: (L - length_MSB + 1)
                        reconstruct_model(xyz, start_vertex) = bitset(reconstruct_model(xyz,start_vertex), bit, bitget(pred_value(xyz), bit));
                    end
                    % diff_bit
                    diff_bit = 1 - bitget(pred_value(xyz), L-length_MSB);
                    reconstruct_model(xyz, start_vertex) = bitset(reconstruct_model(xyz, start_vertex), L-length_MSB, diff_bit);
                else
                    for bit = L: -1: 1
                        reconstruct_model(xyz, start_vertex) = bitset(reconstruct_model(xyz, start_vertex), bit, bitget(pred_value(xyz), bit));
                    end
                end
            end


            curr_neighbors = neighbor(:, start_vertex);
            curr_neighbors = curr_neighbors(curr_neighbors ~= 0);

            % 标记其邻居并设置预测值
            if ~isempty(curr_neighbors)
                recovered_table(curr_neighbors(1)) = 1;
                rev_pred_vertex(:, curr_neighbors(1)) = reconstruct_model(:, start_vertex);

                pred_value = rev_pred_vertex(:, curr_neighbors(1));
                tamper_vertex = curr_neighbors(1);
                for xyz = 1:3
                    length_MSB = reconstruct_length(xyz, tamper_vertex);
                    if length_MSB <= (L - 1)
                        for bit = L: -1: (L - length_MSB + 1)
                            reconstruct_model(xyz, tamper_vertex) = bitset(reconstruct_model(xyz,tamper_vertex), bit, bitget(pred_value(xyz), bit));
                        end
                        % diff_bit
                        diff_bit = 1 - bitget(pred_value(xyz), L-length_MSB);
                        reconstruct_model(xyz, tamper_vertex) = bitset(reconstruct_model(xyz, tamper_vertex), L-length_MSB, diff_bit);
                    else
                        for bit = L: -1: 1
                            reconstruct_model(xyz, tamper_vertex) = bitset(reconstruct_model(xyz,tamper_vertex), bit, bitget(pred_value(xyz), bit));
                        end
                    end
                end
            end

            % 将新起始点的邻居加入队列
            queue = curr_neighbors;
            queue = queue(queue ~= 0);
            continue;  % 继续主循环
        end
    end

    new_queue = [];
    % fprintf('当前队列长度: %d\n', length(queue));

    for i = 1:length(queue)
        current_vertex = queue(i);

        new_queue = [new_queue, neighbor(:, current_vertex)];
        if recovered_table(current_vertex) ~= 1
            mark = 0;
            % 获取current_vertex的一度邻居
            curr_neighbors = neighbor(:, current_vertex);
            curr_neighbors = curr_neighbors(curr_neighbors ~= 0);

            % 找到满足条件的vertex_a
            valid_a = curr_neighbors(recovered_table(curr_neighbors) == 1);

            num_available_neighbor(current_vertex) = length(valid_a);

            if ~isempty(valid_a)
                mirror_vertex = vlist(:, valid_a);
                if size(mirror_vertex, 2) > 2
                    % 计算所有点对之间的距离
                    num_points = size(mirror_vertex, 2);
                    distances = inf(num_points); % 改为 inf 初始化
                    % 构建距离矩阵
                    for p1 = 1:num_points
                        for p2 = p1+1:num_points
                            point1 = round(mirror_vertex(:,p1) * 10^m) / 10^m;
                            point2 = round(mirror_vertex(:,p2) * 10^m) / 10^m;
                            distances(p1,p2) = norm(point1 - point2);
                            distances(p2,p1) = distances(p1,p2);
                        end
                        distances(p1,p1) = -inf; % 设置对角线为 -inf
                    end
                    % 找到最大距离对应的点对
                    [max_dist, idx] = max(distances(:));
                    [point1_idx, point2_idx] = ind2sub(size(distances), idx);

                    % 计算直线方向向量
                    line_dir = mirror_vertex(:,point2_idx) - mirror_vertex(:,point1_idx);
                    if norm(line_dir) == 0  % 如果所有点都重合
                        mirror_point = round(mirror_vertex(:,point1_idx) * 10^m) / 10^m;
                    else
                        % 计算最远两点的中点
                        midpoint = (mirror_vertex(:,point1_idx) + mirror_vertex(:,point2_idx)) / 2;
                        midpoint = round(midpoint * 10^m) / 10^m;

                        % 创建一个索引数组，排除最远的两个点
                        other_points_idx = setdiff(1:num_points, [point1_idx, point2_idx]);

                        % 只对其他点(除了最远的两点)进行镜像变换
                        mirror_points = zeros(3, length(other_points_idx));
                        for ii = 1:length(other_points_idx)
                            p = other_points_idx(ii);
                            point = mirror_vertex(:,p);
                            % 镜像变换：2*中点 - 原点 = 镜像点
                            mirror_points(:,ii) = 2 * midpoint - point;
                            mirror_points(:,ii) = round(mirror_points(:,ii) * 10^m) / 10^m;
                        end

                        if ~isempty(other_points_idx)
                            mirror_point = mean(mirror_points, 2);
                            mirror_point = round(mirror_point * 10^m) / 10^m;
                        else
                            % 如果只有两个最远点，直接使用中点
                            mirror_point = midpoint;
                            mirror_point = round(mirror_point * 10^m) / 10^m;
                        end
                    end

                    if Recovered_LM(LM_index) == 1
                        rev_pred_vertex(:, current_vertex) = max(round(mirror_point), 0);
                        pred_value = rev_pred_vertex(:, current_vertex);
                    else
                        temp_mean = mean(vlist(:, valid_a), 2);
                        temp_mean = round(temp_mean * 10^m) / 10^m;
                        rev_pred_vertex(:, current_vertex) = round(temp_mean);
                        pred_value = rev_pred_vertex(:, current_vertex);
                    end
                    LM_index = LM_index + 1;
                    mark = 1;
                else
                    % 只有一个点的情况
                    mean_point = mean(mirror_vertex, 2);
                    mirror_point = round(mean_point * 10^m) / 10^m;
                end
                if mark == 0
                    rev_pred_vertex(:, current_vertex) = max(round(mirror_point), 0);
                    pred_value = rev_pred_vertex(:, current_vertex);
                end

                recovered_table(current_vertex) = 1;

                for xyz = 1:3
                    length_MSB = reconstruct_length(xyz, current_vertex);
                    if length_MSB <= (L - 1)

                        for bit = L: -1: (L - length_MSB + 1)
                            reconstruct_model(xyz, current_vertex) = bitset(reconstruct_model(xyz,current_vertex), bit, bitget(pred_value(xyz), bit));
                        end
                        % diff_bit
                        diff_bit = 1 - bitget(rev_pred_vertex(xyz, current_vertex), L-length_MSB);
                        reconstruct_model(xyz, current_vertex) = bitset(reconstruct_model(xyz, current_vertex), L-length_MSB, diff_bit);
                    else
                        for bit = L: -1: 1
                            reconstruct_model(xyz, current_vertex) = bitset(reconstruct_model(xyz,current_vertex), bit, bitget(pred_value(xyz), bit));
                        end
                    end
                end

            end
        end

    end

    new_queue = new_queue(new_queue ~= 0);
    new_queue = unique(new_queue);  % 去重
    new_queue = new_queue(recovered_table(new_queue) ~= 1);  % 去除已处理的点
    queue = new_queue;
    test_index = test_index + 1;
end

% recover isolated_vertex
for i = 1:length(isolated_vertices)
    current_vertex = isolated_vertices(i);
    rev_iso_pred = isolated_vertices(i);
    search_radius = 1;

    while num_neighbor_element(rev_iso_pred) == 0
        left_idx = rev_iso_pred - search_radius;
        right_idx = rev_iso_pred + search_radius;

        % 检查右侧点
        if right_idx <= numVertices && num_neighbor_element(right_idx) ~= 0
            rev_iso_pred = right_idx;
            rev_pred_vertex(:, current_vertex) = reconstruct_model(:, rev_iso_pred);
            recovered_table(current_vertex) = 1;
            break;
        end

        % 检查左侧点
        if left_idx >= 1 && num_neighbor_element(left_idx) ~= 0
            rev_iso_pred = left_idx;
            rev_pred_vertex(:, current_vertex) = reconstruct_model(:, rev_iso_pred);
            recovered_table(current_vertex) = 1;
            break;
        end

        search_radius = search_radius + 1;

    end

    pred_value = rev_pred_vertex(:, current_vertex);
    for xyz = 1:3
        length_MSB = reconstruct_length(xyz, current_vertex);
        if length_MSB <= (L - 1)
            for bit = L: -1: (L - length_MSB + 1)
                reconstruct_model(xyz, current_vertex) = bitset(reconstruct_model(xyz,current_vertex), bit, bitget(pred_value(xyz), bit));
            end
            % diff_bit
            diff_bit = 1 - bitget(pred_value(xyz), L-length_MSB);
            reconstruct_model(xyz, current_vertex) = bitset(reconstruct_model(xyz, current_vertex), L-length_MSB, diff_bit);
        else
            for bit = L: -1: 1
                reconstruct_model(xyz, current_vertex) = bitset(reconstruct_model(xyz,current_vertex), bit, bitget(pred_value(xyz), bit));
            end
        end
    end

end

end