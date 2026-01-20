
%% Model owner
clear all
clc
format long g

% Read model
[numVertices,numFaces,flist,vlist] = model_read('bun_zipper.ply');

m = 5;

if m <= 2
    L = 8;
elseif m <= 4
    L = 16;
elseif m <= 9
    L = 32;
end
vlist=roundn(vlist,-m) *10^m;
vlist = round(vlist);

% Find_1_neighbor
[neighbor, num_neighbor_element] = get_neighbors(numFaces, flist, numVertices);
[sorted_num_neighbor_element,I] = sort(num_neighbor_element, 'descend');

% Initialize...
% Record reference vertex
Reference = zeros(1, numVertices);
Reference(I(1)) = 1;
% Record used vertex
used_vertex = zeros(1, numVertices);
used_vertex(I(1)) = 1;
used_vertex(neighbor(1, I(1))) = 1;
% Record predicted result of each vertex
pred_vertex = zeros(3, numVertices);
pred_vertex(:, neighbor(1, I(1))) = vlist(:, I(1));

% Process isolated vertices
[pred_vertex, used_vertex, isolated_vertices] = process_isolated_vertex(sorted_num_neighbor_element, I, num_neighbor_element, pred_vertex, used_vertex, numVertices, vlist);

% Process all vertices
[pred_vertex, used_vertex, Reference, LM, T_number, proportion, distance_mean, distance_vp, distance_hyb, LM_index] = process_vertex_v4(isolated_vertices, pred_vertex, used_vertex, I, neighbor, Reference, numVertices, vlist, m);

% Valid
% sum(used_vertex)

% Get Huffman information
[dict, avglen, length_table, evaluated_EC] = get_huffman_label(m, numVertices, vlist, pred_vertex, Reference);


[rearrange_model, final_EC_bpv, total_EC, acc, sum_length_table, test_table, length_huffman_table] = model_compress_rearrange(L, numVertices, vlist, pred_vertex, dict, Reference, LM);
length_huffman_table


[derived_bits_alice, derived_bits_bob] = get_alice_bob_key;

[encrypted_model, original_key, random_matrix] = model_encryption_v2(L, numVertices, rearrange_model, total_EC, derived_bits_alice);


% Data hider
[data_hiding_key, secret_data, marked_model] = get_marked_model(numVertices, L, encrypted_model);
% trimesh(flist',marked_model(1,:)',marked_model(2,:)',marked_model(3,:)');

% Receiver

% Extract secret data
[extracted_secret_data, extracted_key] = data_extraction(numVertices, L, marked_model, data_hiding_key);

% Valid extracted secret data
isequal(double(extracted_secret_data), secret_data)

% Recover model

[decrypted_model] = model_decryption_v2(numVertices, L, marked_model, extracted_key, derived_bits_alice);

%
[reconstruct_model, reconstruct_length] = model_recovery(flist, numFaces, numVertices, vlist, L, decrypted_model, m, dict, sum_length_table);


isequal(reconstruct_model, vlist)
trimesh(flist',reconstruct_model(1,:)',reconstruct_model(2,:)',reconstruct_model(3,:)');


