function [derived_bits_alice, derived_bits_bob] = get_alice_bob_key

import java.security.KeyPairGenerator
import java.security.spec.ECGenParameterSpec
import java.security.MessageDigest
import javax.crypto.KeyAgreement

kpg = KeyPairGenerator.getInstance('EC');
spec = ECGenParameterSpec('secp256r1');
kpg.initialize(spec);
aliceKeyPair = kpg.generateKeyPair();
alicePrivateKey = aliceKeyPair.getPrivate();
alicePublicKey = aliceKeyPair.getPublic();

bobKeyPair = kpg.generateKeyPair();
bobPrivateKey = bobKeyPair.getPrivate();
bobPublicKey = bobKeyPair.getPublic();

aliceKeyAgreement = KeyAgreement.getInstance('ECDH');
aliceKeyAgreement.init(alicePrivateKey);
aliceKeyAgreement.doPhase(bobPublicKey, true); 
aliceSharedSecret = aliceKeyAgreement.generateSecret();

bobKeyAgreement = KeyAgreement.getInstance('ECDH');
bobKeyAgreement.init(bobPrivateKey);
bobKeyAgreement.doPhase(alicePublicKey, true);
bobSharedSecret = bobKeyAgreement.generateSecret();

md = MessageDigest.getInstance('SHA-256');
derived_key_alice = md.digest(aliceSharedSecret);  % or use bobSharedSecret
derived_bits_alice = zeros(1, 256);
for i = 1:32
    bits = dec2bin(derived_key_alice(i), 8);
    derived_bits_alice((i-1)*8+1:i*8) = bits - '0';
end


md = MessageDigest.getInstance('SHA-256');
derived_key_bob = md.digest(bobSharedSecret);  % or use bobSharedSecret
derived_bits_bob = zeros(1, 256);
for i = 1:32
    bits = dec2bin(derived_key_bob(i), 8);
    derived_bits_bob((i-1)*8+1:i*8) = bits - '0';
end
% length(derived_bits_alice)
% isequal(derived_bits_alice, derived_bits_bob)
end
