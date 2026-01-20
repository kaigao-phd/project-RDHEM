function [reconstruct_model, reconstruct_length] = model_recovery(flist, numFaces, numVertices, vlist, L, decrypted_model, m,  dict, sum_length_table)


%依次取出就行
reconstruct_length = zeros(3, numVertices);
pointer_bit = L;
pointer_xyz = 1;
pointer_num = 1;

%这里设置每个code的长度用7bits来记录，可以更短
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


%recover res_bit reference vertices的全部bit会被直接复原

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

% pointer_num - 1
% isequal(round(reconstruct_model), round(valid_res_bit))





% recover MSB

recovered_table = zeros(1, numVertices);
recovered_table(I(1)) = 1;
recovered_table(neighbor(1, I(1))) = 1;



for 你的方法

    %重新得到每个顶点的预测值

    pred_value = rev_pred_vertex(:, current_vertex);

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