function [ A ] = importA(filename,bands,N)

% fid = fopen(filename);
% tmp = textscan(fid, '%f', 'Delimiter', ',', 'BufSize', bufsize);
% fclose(fid);
% tmp = tmp{1};

fid = fopen(filename);
tmp = fread(fid, [N*length(bands),1], 'float32');
fclose(fid);

tmp = reshape(tmp, [N,length(bands)]);

% NOTE: cannot do this with spdiags() since the gpumatting representation
% and matlab are different, and Levin's Laplacian is not really symmetric

i = repmat(int32((1:N)'),[1,length(bands)]);
j = i + repmat(int32(bands),[N,1]);

% Just clip off all the bad indices.
i(i<1) = 1;
i(i>N) = N;
j(j<1) = 1;
j(j>N) = N;

%A = sparse(i(:),j(:),tmp(:),N,N);
% Have to do this to use integer (instead of double) indices.
A = accumarray([i(:), j(:)], tmp(:), [N,N], [], [], true);

end