%% Test the convergence rate of various solvers at different resolutions.

iter = 10000;
numRes = 3;

Lres = cell(numRes,1);
Ares = cell(numRes,1);
bres = cell(numRes,1);
x0res = cell(numRes,1);
alphagtres = cell(numRes,1);
scribsres = cell(numRes,1);
upres = cell(numRes,1);
downres = cell(numRes,1);

Itmp = I;
alphagttmp = alphagt(:,:,1);
x0tmp = reshape(x0,[532,800]);
scribstmp = scribs;

for i=numRes:-1:1
    Itmp = imresize(Itmp,0.5);
    alphagttmp = imresize(alphagttmp,0.5);
    alphagtres{i} = alphagttmp;
    scribstmp = imresize(scribstmp,0.5);
    scribsres{i} = scribstmp;
    x0tmp = imresize(x0tmp,0.5);
    x0res{i} = x0tmp(:);
    
    [ctmp,vtmp] = scribData(scribsres{i});
    ctmp = ctmp(:);
    vtmp = vtmp(:);

    Lres{i} = getLevinLap(Itmp, eps, 1);
    Ares{i} = Lres{i}+gamma*spdiags(ctmp,[0],size(Lres{i},1),size(Lres{i},2));
    bres{i} = gamma*vtmp;

    %N = size(I,1)*size(I,2);

    [upres{i},downres{i}] = dataTransferOperators(Itmp,Ares{i},4);
end

clear Itmp alphatmp ctmp vtmp x0tmp

resCG = cell(numRes,1);
mseCG = cell(numRes,1);
resVcycle = cell(numRes,1);
mseVcycle = cell(numRes,1);

for i=1:numRes
    [ ~, resCG{i}, mseCG{i} ] = multiGridRoughCG( x0res{i}, Ares{i}, bres{i}, {}, {}, 10000, alphagtres{i} );
    [ ~, resVcycle{i}, mseVcycle{i} ] = vcycle( x0res{i}, Ares{i}, bres{i}, upres{i}, downres{i}, 10000, alphagtres{i} );
end
