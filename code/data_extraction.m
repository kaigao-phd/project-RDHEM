function [extracted_secret_data, extracted_key] = data_extraction(numVertices, L, marked_model, data_hiding_key)

pointer_bit = 1;
pointer_xyz = 3;
pointer_num = numVertices;
position_index = 1;
extracted_position_label = [];

while position_index <= ceil(log2(3*L*numVertices))
    valid_bit = bitget(marked_model(pointer_xyz, pointer_num), pointer_bit);
    extracted_position_label = [extracted_position_label, valid_bit];
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

extracted_key = [];
key_index = 1;
while key_index <= 256
    valid_bit = bitget(marked_model(pointer_xyz, pointer_num), pointer_bit);
    extracted_key = [extracted_key, valid_bit];
    key_index = key_index + 1;
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

extracted_position_label = sprintf('%d', extracted_position_label);
EC = bin2dec(extracted_position_label);

extracted_data = [];
secret_index = 1;
while secret_index <= EC
    valid_bit = bitget(marked_model(pointer_xyz, pointer_num), pointer_bit);
    extracted_data = [extracted_data, valid_bit];
    secret_index = secret_index + 1;
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

rng(data_hiding_key)
random_data = randi([0,1], 1, EC);
extracted_secret_data = xor(extracted_data, random_data);
