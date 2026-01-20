function [derived_bits_alice, derived_bits_bob] = get_alice_bob_key

import java.security.KeyPairGenerator
import java.security.spec.ECGenParameterSpec
import java.security.MessageDigest
import javax.crypto.KeyAgreement

% 为Alice生成EC密钥对
kpg = KeyPairGenerator.getInstance('EC');
spec = ECGenParameterSpec('secp256r1');
kpg.initialize(spec);
aliceKeyPair = kpg.generateKeyPair();
alicePrivateKey = aliceKeyPair.getPrivate();
alicePublicKey = aliceKeyPair.getPublic();

% 为Bob生成EC密钥对
bobKeyPair = kpg.generateKeyPair();
bobPrivateKey = bobKeyPair.getPrivate();
bobPublicKey = bobKeyPair.getPublic();

% Alice的密钥协商
aliceKeyAgreement = KeyAgreement.getInstance('ECDH');
aliceKeyAgreement.init(alicePrivateKey);
aliceKeyAgreement.doPhase(bobPublicKey, true);  % 使用Bob的公钥
aliceSharedSecret = aliceKeyAgreement.generateSecret();

% Bob的密钥协商
bobKeyAgreement = KeyAgreement.getInstance('ECDH');
bobKeyAgreement.init(bobPrivateKey);
bobKeyAgreement.doPhase(alicePublicKey, true);  % 使用Alice的公钥
bobSharedSecret = bobKeyAgreement.generateSecret();

% 使用SHA-256导出256位密钥（Alice和Bob会得到相同的结果）
md = MessageDigest.getInstance('SHA-256');
derived_key_alice = md.digest(aliceSharedSecret);  % 或使用bobSharedSecret，结果相同
derived_bits_alice = zeros(1, 256);
for i = 1:32
    bits = dec2bin(derived_key_alice(i), 8);
    derived_bits_alice((i-1)*8+1:i*8) = bits - '0';
end

% 使用SHA-256导出256位密钥（Alice和Bob会得到相同的结果）
md = MessageDigest.getInstance('SHA-256');
derived_key_bob = md.digest(bobSharedSecret);  % 或使用bobSharedSecret，结果相同
derived_bits_bob = zeros(1, 256);
for i = 1:32
    bits = dec2bin(derived_key_bob(i), 8);
    derived_bits_bob((i-1)*8+1:i*8) = bits - '0';
end
% length(derived_bits_alice)
% isequal(derived_bits_alice, derived_bits_bob)
end