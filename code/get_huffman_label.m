function [dict, avglen, length_table, evaluated_EC] = get_huffman_label(m, numVertices, vlist, pred_vertex, Reference)

if m <= 2
    L = 8;
elseif m <= 4
    L = 16;
elseif m<= 9
    L = 32;
end


msb_vnew = round(vlist);
msb_vpred = round(pred_vertex);


EC = 0;
length_table = zeros(1, L+1);
for i = 1 : numVertices
    if Reference(i) ~= 1
        for xyz = 1:3
            bit_plane = L;
            E2 = 0;
            while (bit_plane > 0) && (bitget(msb_vnew(xyz, i), bit_plane) ==  bitget(msb_vpred(xyz, i), bit_plane))
                EC = EC + 1;
                E2 = E2 + 1;
                bit_plane = bit_plane - 1;
            end

            if E2 == 0
                length_table(L+1) = length_table(L+1) + 1;
            else
                length_table(E2) = length_table(E2) + 1;
            end

        end
    end
end


% huffman coding
symbols = (1: (L+1));
prob = zeros(1, (L+1));
for i = 1: (L+1)
    prob(i) = (length_table(i)/ sum(length_table));
end
[dict,avglen] = huffmandict(symbols,prob);

max_length = 0;
for i = 1:length(dict)
    if length(dict{i, 2}) > max_length
        max_length = length(dict{i, 2});
    end
end

reference_num = sum(Reference);
evaluated_EC = (EC - (numVertices - reference_num)*3*avglen + (numVertices - reference_num)*3-length_table(L)-length_table(L+1)*length(dict{L+1, 2}) - length(dict)*7 - 256 - ceil(log2(3*L*numVertices)))/ numVertices;

end