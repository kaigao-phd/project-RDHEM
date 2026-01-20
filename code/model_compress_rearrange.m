function [rearrange_model, final_EC_bpv, total_EC, acc, sum_length_table, test_table, length_huffman_table] = model_compress_rearrange(L, numVertices, vlist, pred_vertex, dict, Reference, LM)

msb_vnew = round(vlist);
msb_vpred = round(pred_vertex);


msb_hiding = msb_vnew;
sum_length_table = zeros(3, numVertices);
rearrange_model = zeros(3, numVertices);
pointer_bit = L;
pointer_xyz = 1;
pointer_num = 1;

% Hide Huffman table
bit_stream = [];
code_length_param = 7;
for i = 1: L+1
    current_code_length = length(dict{i, 2});
    current_code_length_bin = double(dec2bin(current_code_length, code_length_param)) - '0';
    bit_stream = [bit_stream, current_code_length_bin, dict{i, 2}];
end

length_huffman_table = length(bit_stream)

for i = 1: length(bit_stream)
    rearrange_model(pointer_xyz, pointer_num) = bitset(rearrange_model(pointer_xyz, pointer_num), pointer_bit, bit_stream(i));
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

fprintf('huffman_table_position = %d\n', pointer_num)

% Hide label
test_table = zeros(1, L+1);
for i = 1 : numVertices
    if Reference(i) ~= 1
        for xyz = 1:3
            bit_plane = L;
            E2 = 0;
            while (bit_plane > 0) && (bitget(msb_vnew(xyz, i), bit_plane) ==  bitget(msb_vpred(xyz, i), bit_plane))

                E2 = E2 + 1;
                bit_plane = bit_plane - 1;

            end
            sum_length_table(xyz, i) = E2;
            if E2 == 0
                bits_to_set = length(dict{L+1, 2});
                test_table(L+1) = test_table(L+1) + 1;
            else
                bits_to_set = length(dict{E2, 2});
                test_table(E2) = test_table(E2) + 1;
            end
            code_bit = 1;
            for bit = L:-1:L - bits_to_set + 1
                msb_hiding(xyz, i) = bitset(msb_hiding(xyz, i), bit, dict{E2, 2}(code_bit));
                code_bit = code_bit + 1;
                assigned_bit = bitget(msb_hiding(xyz, i), bit);
                rearrange_model(pointer_xyz, pointer_num) = bitset(rearrange_model(pointer_xyz, pointer_num), pointer_bit, assigned_bit);
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
            for bit = L - bits_to_set :-1: L - E2 + 1
                msb_hiding(xyz, i) = bitset(msb_hiding(xyz, i), bit, 1);

            end
        end
    end
end

fprintf('label_position = %d\n', pointer_num)


% Hide res_bit
for i = 1:numVertices
    if Reference(i) ~= 1
        for xyz = 1:3
            E2 = sum_length_table(xyz, i);
            if E2 <= (L-2)
                for bit = L - E2 - 1 :-1 : 1
                    res_bit = bitget(msb_vnew(xyz, i), bit);
                    rearrange_model(pointer_xyz, pointer_num) = bitset(rearrange_model(pointer_xyz, pointer_num), pointer_bit, res_bit);
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
            end
        end
    else
        for xyz = 1:3
            for bit = L :-1 : 1
                res_bit = bitget(msb_vnew(xyz, i), bit);
                rearrange_model(pointer_xyz, pointer_num) = bitset(rearrange_model(pointer_xyz, pointer_num), pointer_bit, res_bit);
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
        end
    end
end

fprintf('res_position = %d\n', pointer_num)

% Hide LM
LM_index = 1;
while LM_index <= length(LM)
    rearrange_model(pointer_xyz, pointer_num) = bitset(rearrange_model(pointer_xyz, pointer_num), pointer_bit, LM(LM_index));
    LM_index = LM_index + 1;
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

fprintf('LM_position = %d\n', pointer_num)

% Fill in redundant bits
EC = 0;
while pointer_num <= numVertices
    rearrange_model(pointer_xyz, pointer_num) = bitset(rearrange_model(pointer_xyz, pointer_num), pointer_bit, randi([0, 1]));
    EC = EC + 1;
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


reference_num = sum(Reference);
acc = EC/((numVertices - reference_num)*3*L)*100;
position_label_length = ceil(log2(3*L*numVertices));
key_length = 256;
final_EC_bpv = (EC - position_label_length - key_length) / numVertices;
total_EC = dec2bin((EC - position_label_length - key_length), position_label_length);

position_index = 1;

pointer_bit = 1;
pointer_xyz = 3;
pointer_num = numVertices;
while position_index <= position_label_length
    bit_value = double(total_EC(position_index)) - '0';
    rearrange_model(pointer_xyz, pointer_num) = bitset(rearrange_model(pointer_xyz, pointer_num), pointer_bit, bit_value);
    position_index = position_index + 1;
    pointer_bit = pointer_bit + 1;
    if pointer_bit > L
        pointer_xyz = pointer_xyz - 1;
        pointer_bit = 1;
        if pointer_xyz < 1
            pointer_num = pointer_num - 1;
            pointer_bit = 1;
            pointer_xyz = 3;
        end
    end

end

end