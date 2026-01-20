function [pred_vertex, used_vertex, Reference, LM, T_number, proportion, distance_mean, distance_vp, distance_hyb, LM_index] = process_vertex_v4(isolated_vertices, pred_vertex, used_vertex, I, neighbor, Reference, numVertices, vlist, m)
% 初始化队列
queue = neighbor(:, I(1));
queue = queue(queue ~= 0);

num_available_neighbor = zeros(1, numVertices);
test_index = 1;

T_number = 0;
P_number = 0;
distance_vp = 0;
distance_mean = 0;
distance_hyb = 0;
LM = [];
LM_index = 1;

% 主循环，处理所有顶点
while ~all(used_vertex(:) == 1)
    % 如果当前队列为空，说明一个子模型处理完成，需要寻找新的起始点
    if isempty(queue)
        % 找到第一个未处理的顶点
        start_vertex = find(used_vertex == 0, 1);
        if ~isempty(start_vertex)
            % 初始化这个新的起始点

            search_radius = 1;

            while used_vertex(start_vertex) == 0
                % 检查当前半径的左右两个位置
                left_idx = start_vertex - search_radius;
                right_idx = start_vertex + search_radius;

                % 检查右侧点
                if right_idx <= numVertices && used_vertex(right_idx) == 1 && ~ismember(right_idx, isolated_vertices)
                    pred_vertex(:, start_vertex) = vlist(:, right_idx);
                    used_vertex(start_vertex) = 1;
                    break;
                end

                % 检查左侧点
                if left_idx >= 1 && used_vertex(left_idx) == 1 && ~ismember(left_idx, isolated_vertices)
                    pred_vertex(:, start_vertex) = vlist(:, left_idx);
                    used_vertex(start_vertex) = 1;
                    break;
                end

                search_radius = search_radius + 1;
            end

            curr_neighbors = neighbor(:, start_vertex);
            curr_neighbors = curr_neighbors(curr_neighbors ~= 0);

            % 标记其邻居并设置预测值
            if ~isempty(curr_neighbors)
                used_vertex(curr_neighbors(1)) = 1;
                pred_vertex(:, curr_neighbors(1)) = vlist(:, start_vertex);
            end

            % 将新起始点的邻居加入队列
            queue = curr_neighbors;
            queue = queue(queue ~= 0);
            continue;  % 继续主循环
        end
    end

    new_queue = [];
    %     fprintf('当前队列长度: %d\n', length(queue));

    for i = 1:length(queue)
        current_vertex = queue(i);
        new_queue = [new_queue, neighbor(:, current_vertex)];

        if used_vertex(current_vertex) ~= 1
            mark = 0;
            % 获取current_vertex的一度邻居
            curr_neighbors = neighbor(:, current_vertex);
            curr_neighbors = curr_neighbors(curr_neighbors ~= 0);

            % 找到满足条件的vertex_a
            valid_a = curr_neighbors(used_vertex(curr_neighbors) == 1);
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
                    T_number = T_number + 1;
                    temp_mean = mean(vlist(:, valid_a), 2);
                    temp_mean = round(temp_mean * 10^m) / 10^m;
                    sum_1 = abs(temp_mean(1) - vlist(1, current_vertex)) + abs(temp_mean(2) - vlist(2, current_vertex)) + abs(temp_mean(3) - vlist(3, current_vertex));
                    temp_vp = round(mirror_point);
                    temp_vp = max(temp_vp(:), 0);
                    sum_2 = abs(temp_vp(1) - vlist(1, current_vertex)) + abs(temp_vp(2) - vlist(2, current_vertex)) + abs(temp_vp(3) - vlist(3, current_vertex));
                    distance_mean = distance_mean + sum_1;
                    distance_vp = distance_vp + sum_2;
                    if sum_1 > sum_2
                        P_number = P_number + 1;
                        pred_vertex(:, current_vertex) = round(mirror_point);
                        pred_vertex(:, current_vertex) = max(pred_vertex(:, current_vertex), 0);
                        LM(LM_index) = 1;
                        LM_index = LM_index + 1;

                    else
                        pred_vertex(:, current_vertex) = round(temp_mean);
                        LM(LM_index) = 0;
                        LM_index = LM_index + 1;

                    end
                    sum_3 = abs(pred_vertex(1, current_vertex) - vlist(1, current_vertex)) + abs(pred_vertex(2, current_vertex) - vlist(2, current_vertex)) + abs(pred_vertex(3, current_vertex) - vlist(3, current_vertex));
                    distance_hyb = distance_hyb + sum_3;
                    mark = 1;

                else
                    % 只有一个点的情况
                    mean_point = mean(mirror_vertex, 2);
                    mirror_point = round(mean_point * 10^m) / 10^m;
                end
                if mark == 0
                    pred_vertex(:, current_vertex) = round(mirror_point);
                    pred_vertex(:, current_vertex) = max(pred_vertex(:, current_vertex), 0);
                end
            end
            used_vertex(current_vertex) = 1;
        end
    end

    new_queue = new_queue(new_queue ~= 0);
    new_queue = unique(new_queue);  % 去重
    new_queue = new_queue(used_vertex(new_queue) ~= 1);  % 去除已处理的点
    queue = new_queue;
    test_index = test_index + 1;
end
proportion = P_number / T_number;
distance_mean = distance_mean / T_number;
distance_vp = distance_vp / T_number;
distance_hyb = distance_hyb / T_number;
LM_index = LM_index - 1;

end