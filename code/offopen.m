function [numVertices,numFaces,flist,vlist]=offopen(fid)
if fid<0
    error('open faild!');
end
string = fgetl(fid);
num = fgetl(fid);
count = sscanf(num,'%d %d %d');
numVertices = count(1);
numFaces = count(2);
vlist = zeros(3,numVertices);
for vnum = 1:numVertices
    vert = fgetl(fid);
    vdata = sscanf(vert,'%f %f %f');
    vlist(:,vnum) = vdata;
end
flist = zeros(3,numFaces);
for fnum = 1:numFaces
    face = fgetl(fid);
    fdata = sscanf(face,'%d %d %d %d');
    flist(:,fnum) = fdata(2:end)+1;
end
flist=sort(flist);

% for i = 1:numFaces-1
%     temp = flist(:, 1);
%     for j = 1:numFaces-1-i
%         if flist(1, j) > flist(1, j+1)
%             temp(:, 1) = flist(:, j+1);
%             flist(:, j+1) = flist(:, j);
%             flist(:, j) = temp(:, 1);
%         elseif flist(1, j) == flist(1, j+1)
%             if flist(2, j) > flist(2, j+1)
%                 temp(:, 1) = flist(:, j+1);
%                 flist(:, j+1) = flist(:, j);
%                 flist(:, j) = temp(:, 1);
%             elseif flist(2, j) == flist(2, j+1)
%                 if flist(3, j) > flist(3, j+1)
%                     temp(:, 1) = flist(:, j+1);
%                     flist(:, j+1) = flist(:, j);
%                     flist(:, j) = temp(:, 1);
%                 end
%             end
%         end
%     end
% end

% 
% trimesh(flist',vlist(1,:)',vlist(2,:)',vlist(3,:)'); 
% axis off
% box off