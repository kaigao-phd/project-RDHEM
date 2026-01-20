function [decrypted_model] = model_decryption_v2(numVertices, L, marked_model, extracted_key, derived_bits_alice)

% recover model encryption key
model_encryption_key = xor(extracted_key, derived_bits_alice);

total_bits = numVertices * 3 * L; % 计算模型大小
% 初始化随机流
random_stream = zeros(1, total_bits);

% 将256位密钥转换为字节数组格式
key_bytes = zeros(1, 32, 'uint8');
for i = 1:32
    byte_val = bin2dec(char(model_encryption_key((i-1)*8+1:i*8) + '0'));
    key_bytes(i) = uint8(byte_val);
end

% 导入Java类
import javax.crypto.Cipher;
import javax.crypto.spec.SecretKeySpec;
import javax.crypto.spec.IvParameterSpec;

% 创建AES密钥规格对象
key_spec = SecretKeySpec(key_bytes, 'AES');

% 创建初始向量（全0）
iv = zeros(1, 16, 'uint8');
iv_spec = IvParameterSpec(iv);

% 获取AES/CTR/NoPadding加密器实例
cipher = Cipher.getInstance('AES/CTR/NoPadding');
cipher.init(Cipher.ENCRYPT_MODE, key_spec, iv_spec);

% 计算需要生成的块数（每块16字节）
num_blocks = ceil(total_bits / 128);
pos = 1;

% 生成密钥流
for i = 1:num_blocks
    % 创建全0块进行加密
    zero_block = zeros(1, 16, 'uint8');
    encrypted_block = cipher.update(zero_block);
    
    % 将加密块转换为比特流
    encrypted_bits = zeros(1, 128);
    for j = 1:16
        byte_val = typecast(encrypted_block(j), 'uint8'); % 确保是uint8类型
        for k = 1:8
            encrypted_bits((j-1)*8+k) = bitget(byte_val, 9-k); % 从最高位开始
        end
    end
    
    % 计算当前轮次可以填充的位数
    remaining = total_bits - pos + 1;
    fill_length = min(128, remaining);
    
    % 填充随机流
    random_stream(pos:pos+fill_length-1) = encrypted_bits(1:fill_length);
    
    % 更新位置
    pos = pos + fill_length;
end

% 重塑随机流为模型大小
random_matrix = reshape(random_stream, [L, 3, numVertices]);
random_matrix = permute(random_matrix, [2 3 1]);

% validate
% random_matrix_2 = randi([0, 1], [3, numVertices, L]);
% random_matrix = random_matrix_2;

% 初始化输出
decrypted_model = marked_model;  % 先复制所有值

% 执行加密，但跳过保留位
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