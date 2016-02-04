N1d = 256;

tmpB = zeros(N1d,3);

% Simple Laplacian
% tmpB(:,1) = -ones(N1d,1);
% tmpB(1:end,2) = 2*ones(N1d,1);
% tmpB(:,3) = -ones(N1d,1);

% Without boundary effects
tmpB(:,1) = -ones(N1d,1);
tmpB(2:end-1,2) = 2*ones(N1d-2,1); tmpB(1,2) = 1; tmpB(end,2) = 1;
tmpB(:,3) = -ones(N1d,1);

% Create sparse matrix
A1d = spdiags(tmpB,[-1,0,1],N1d,N1d);
% Create rhs
b1d = zeros(N1d,1);
% Create initial solution x0
t1d = linspace(0,1,N1d)';
x01d = sin( 2*pi*t1d );

tmpn = N1d;
% NOTE: always doing the same number of small grids does NOT result in
% convergence that is independent of N1d! This only happens when the
% smallest grid is the same size, so the number of levels is O(log2(N1d)).
%up1d = cell(4,1);
%down1d = cell(4,1);
%for i=4:-1:1
% Always have the smallest grid be size 4.
up1d = cell(log2(N1d)-2,1);
down1d = cell(log2(N1d)-2,1);
for i=log2(N1d)-2:-1:1
    tmpijk = zeros(3*tmpn/2,3);
    % Beginning index.
    tmpijk(1:2,:) = [[1;1],[1;2],[1;0.5]];
    ijkndx = 3;
    bigndx = 3;
    for smallndx=2:tmpn/2
        tmpijk(ijkndx:ijkndx+2,:) = [repmat(smallndx,3,1),(bigndx-1:bigndx+1)',[0.5,1,0.5]'];
        ijkndx = ijkndx + 3;
        bigndx = bigndx + 2;
    end
    down1d{i} = sparse(tmpijk(1:ijkndx-1,1), tmpijk(1:ijkndx-1,2), tmpijk(1:ijkndx-1,3), tmpn/2, tmpn);
    up1d{i} = down1d{i}';
    tmpn = tmpn/2;
end

clear tmpB t1d tmpn tmpijk ijkndx bigndx;