function [encrypted_model, original_key, random_matrix] = model_encryption_v2(L, numVertices, rearrange_model, total_EC, derived_bits_alice)

% Generate 256bits model key
original_key = randi([0 1], 1, 256);
total_bits = numVertices * 3 * L; 

random_stream = zeros(1, total_bits);

key_bytes = zeros(1, 32, 'uint8');
for i = 1:32
    byte_val = bin2dec(char(original_key((i-1)*8+1:i*8) + '0'));
    key_bytes(i) = uint8(byte_val);
end

% 导入Java类
import javax.crypto.Cipher;
import javax.crypto.spec.SecretKeySpec;
import javax.crypto.spec.IvParameterSpec;

key_spec = SecretKeySpec(key_bytes, 'AES');

iv = zeros(1, 16, 'uint8');
iv_spec = IvParameterSpec(iv);

cipher = Cipher.getInstance('AES/CTR/NoPadding');
cipher.init(Cipher.ENCRYPT_MODE, key_spec, iv_spec);

num_blocks = ceil(total_bits / 128);
pos = 1;

for i = 1:num_blocks
    zero_block = zeros(1, 16, 'uint8');
    encrypted_block = cipher.update(zero_block);
    
    encrypted_bits = zeros(1, 128);
    for j = 1:16
        byte_val = typecast(encrypted_block(j), 'uint8'); 
        for k = 1:8
            encrypted_bits((j-1)*8+k) = bitget(byte_val, 9-k); 
        end
    end
    
    remaining = total_bits - pos + 1;
    fill_length = min(128, remaining);
    
    random_stream(pos:pos+fill_length-1) = encrypted_bits(1:fill_length);
    
    pos = pos + fill_length;
end

random_matrix = reshape(random_stream, [L, 3, numVertices]);
random_matrix = permute(random_matrix, [2 3 1]);

position_label = ceil(log2(3*L*numVertices));
reserved_bits = 256 + position_label;

encrypted_model = rearrange_model; 

pointer_bit = 1;
pointer_xyz = 3;
pointer_num = numVertices;
reserved_positions = zeros(reserved_bits, 3);

for i = 1:reserved_bits
    reserved_positions(i,:) = [pointer_xyz pointer_num pointer_bit];
    pointer_bit = pointer_bit + 1;
    if pointer_bit > L
        pointer_xyz = pointer_xyz - 1;
        pointer_bit = 1;
        if pointer_xyz < 1
            pointer_num = pointer_num - 1;
            pointer_xyz = 3;
        end
    end
end

for i = 1:numVertices
    for xyz = 1:3
        for bit = 1:L
            current_pos = [xyz i bit];
            is_reserved = any(all(reserved_positions == current_pos, 2));

            if ~is_reserved
                orig_bit = bitget(rearrange_model(xyz, i), bit);
                rand_bit = bitget(random_matrix(xyz, i), bit);
                encrypted_model(xyz, i) = bitset(encrypted_model(xyz, i), bit, bitxor(orig_bit, rand_bit));
            end
        end
    end
end

% Hiding position label and encrypted key

position_index = 1;
position_label_length = ceil(log2(3*L*numVertices));

pointer_bit = 1;
pointer_xyz = 3;
pointer_num = numVertices;
while position_index <= position_label_length
    bit_value = double(total_EC(position_index)) - '0'; 
    encrypted_model(pointer_xyz, pointer_num) = bitset(encrypted_model(pointer_xyz, pointer_num), pointer_bit, bit_value);
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

encrypted_key = xor(original_key, derived_bits_alice);

key_index = 1;
while key_index <= length(encrypted_key)
    encrypted_model(pointer_xyz, pointer_num) = bitset(encrypted_model(pointer_xyz, pointer_num), pointer_bit, encrypted_key(key_index));
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

end
