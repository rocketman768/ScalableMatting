names = {'doll','donkey','elephant','net','pineapple','plant','plasticbag','troll'};
eps = 1e-6;
gamma = 10;

for ii=1:length(names)
    I = double(imread(sprintf('Evaluation/%s.png',names{ii})))./255;
    load(sprintf('Evaluation/%s_labels.mat',names{ii}));
    
    for jj=1:3
        scribs = double(imread(sprintf('Evaluation/Trimap%d/%s.png',jj,names{ii})))./255;
        scribs = scribs(:,:,1);
        [c,v] = scribData(scribs);
        c = c(:);
        v = v(:);

        L = getLevinLap(I, eps, 1);
        A = L+gamma*spdiags(c,[0],size(L,1),size(L,2));
        b = gamma*v;

        N = size(I,1)*size(I,2);

        [up8,down8] = naiveTransferOperators(I,8);
        [upmerge,downmerge] = mergeTransferOperators(I,labels,4);
        
        % Just solve the problem directly on the small scale segmentation space.
        x0_segment = upmerge{4}*((downmerge{4}*(A*upmerge{4}))\(downmerge{4}*b));
        
        % Vcycle for a while.
        x = x0_segment;
        for iter=1:1000
            fprintf(1,'im: %s jj: %d iter: %.4d\n', names{ii}, jj, iter);
            x = vcycle( x, A, b, up8, down8, 8 );
        end
        
        myImwrite( reshape(x,[size(I,1),size(I,2)]), sprintf('Evaluation/Output/Trimap%d/%s.png',jj,names{ii}), [0,1] );
    end
end
