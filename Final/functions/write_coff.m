function write_coff(filename, vertex, face);

% write_off - write a mesh to a OFF file
%
%   write_off(filename, vertex, face);
%
%   vertex must be of size [n,3]
%   face must be of size [p,3]
%
%   Copyright (c) 2003 Gabriel Peyr


if size(vertex,2)~=7
    vertex=vertex';
end
if size(vertex,2)~=7
    error('vertex does not have the correct format.');
end


if size(face,2)~=3
    face=face';
end
if size(face,2)~=3
    error('face does not have the correct format.');
end

fid = fopen(filename,'wt');
if( fid==-1 )
    error('Can''t open the file.');
    return;
end

% header
fprintf(fid, 'COFF\n');
fprintf(fid, '%d %d 0\n', size(vertex,1), size(face,1));

% write the points & faces
fprintf(fid, '%f %f %f %d %d %d %d\n', vertex');
fprintf(fid, '3 %d %d %d\n', face'-1);

fclose(fid);