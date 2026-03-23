function [decrypted_model] = model_decryption_v2(numVertices, L, marked_model, extracted_key, derived_bits_alice)

% recover model encryption key
model_encryption_key = xor(extracted_key, derived_bits_alice);

total_bits = numVertices * 3 * L;

random_stream = zeros(1, total_bits);

key_bytes = zeros(1, 32, 'uint8');
for i = 1:32
    byte_val = bin2dec(char(model_encryption_key((i-1)*8+1:i*8) + '0'));
    key_bytes(i) = uint8(byte_val);
end

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

% validate
% random_matrix_2 = randi([0, 1], [3, numVertices, L]);
% random_matrix = random_matrix_2;

decrypted_model = marked_model;  

for i = 1:numVertices
    for xyz = 1:3
        for bit = 1:L
            enc_bit = bitget(marked_model(xyz, i), bit);
            rand_bit = bitget(random_matrix(xyz, i), bit);
            decrypted_model(xyz, i) = bitset(decrypted_model(xyz, i), bit, bitxor(enc_bit, rand_bit));
        end
    end
end

end
