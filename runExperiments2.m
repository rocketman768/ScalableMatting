gamma = 10;
numLevels = 4;
iterations = 100;
t = 1:iterations+1;

graymap = zeros(256,3);
for i=1:256
    graymap(i,:) = (i-1)*ones(1,3)/255;
end

for imnum=2:27
    fprintf(1,'==Experiment %.2d==\n',imnum);
    
    imname = sprintf('experiments/%.2d.png',imnum);
    scribname = sprintf('experiments/dbscribs/%.2d_scribs.png',imnum);
    gtname = sprintf('experiments/%.2d_gt.png',imnum);
    labelname = sprintf('experiments/%.2d_labels.mat',imnum);
    
    % Load shit.
    I = double(imread(imname))./255;
    % Size shit.
    rows = size(I,1);
    cols = size(I,2);
    N = rows*cols;
    gt = double(imread(gtname))./255;
    gt = gt(:,:,1);
    scribs = double(imread(scribname))./255;
    scribs = scribs(:,:,1);
    % DANGEROUS!
    scribs = imresize(scribs, [rows,cols]);
    %labels = load(labelname);
    load(labelname);
    
    % Scribbles and shit.
    [c,b] = scribData(scribs);
    c = reshape(c,[N,1]);
    b = gamma*reshape(b,[N,1]);
    
    [up,down] = mergeTransferOperators(I, labels, numLevels);
    [ upnaive, downnaive ] = naiveTransferOperators( I, numLevels );
    
    % Huge sparse matrix shit.
    A = getLevinLap(I, 1e-4, 1) + gamma*spdiags(c,[0], N, N);
    
    % Initial point.
    %x0 = 0.5*ones([N,1]);
    Asmall = down{numLevels}*A*up{numLevels};
    bsmall = down{numLevels}*b;
    x0 = up{numLevels}*(Asmall\bsmall);
    
    % Gradient descent for various levels.
    rgrad = zeros(iterations+1,numLevels+1);
    rcg = zeros(iterations+1,numLevels+1);
    for ell = numLevels+1:-1:1
        upsmall = {};
        downsmall = {};
        for i = ell:numLevels
            upsmall{end+1} = upnaive{i};
            downsmall{end+1} = downnaive{i};
        end
        
        % =============Gradient Descent================
        [xgrad, rgrad(:,numLevels+1-ell+1), xk] = multiGridRoughGradient(x0,A,b,upsmall,downsmall,iterations,gt);
        
        myImwrite(reshape(xgrad,[rows,cols]), sprintf('results2/%.2d_grad%d.png',imnum,numLevels+1-ell+1), [0,1]);
        
        % Make a movie.
%         aviobj = avifile(sprintf('results/%.2d_grad%d.avi',imnum,numLevels+1-ell+1),'COLORMAP',graymap);
%         for framenum=1:iterations+1
%             aviobj = addframe(aviobj,myIm2Frame(reshape(xk(:,framenum),[rows,cols]),graymap));
%         end
%         aviobj = close(aviobj);
        
        % ================CG Descent===================
        [xcg, rcg(:,numLevels+1-ell+1), xk] = multiGridRoughCG(x0,A,b,upsmall,downsmall,iterations,gt);
        
        myImwrite(reshape(xcg,[rows,cols]), sprintf('results2/%.2d_cg%d.png',imnum,numLevels+1-ell+1), [0,1]);
        
        % Make a movie.
%         aviobj = avifile(sprintf('results/%.2d_cg%d.avi',imnum,numLevels+1-ell+1),'COLORMAP',graymap);
%         for framenum=1:iterations+1
%             aviobj = addframe(aviobj,myIm2Frame(reshape(xk(:,framenum),[rows,cols]),graymap));
%         end
%         aviobj = close(aviobj);
        
        % ==================Vcycle=====================
        if( ell < numLevels )
            xv = x0;
            rv = zeros([iterations,1]);
            for i=1:100
                xv = vcycle( xv, A, b, upsmall, downsmall, numLevels+1-ell );
                rv(i) = sqrt(norm(xv-gt(:),2)/N);
            end
            myImwrite(reshape(xv,[rows,cols]), sprintf('results2/%.2d_vcycle%d.png',imnum,numLevels+1-ell+1), [0,1]);
        end
    end
    
    % ================Combined Descent===================
    [xcomb, rcombined, xk] = multiGridRoughCG(x0,A,b,upsmall,downsmall,iterations,gt,19);
    
    myImwrite(reshape(xcomb,[rows,cols]), sprintf('results2/%.2d_comb.png',imnum), [0,1]);
    
    % Make a movie.
%     aviobj = avifile(sprintf('results/%.2d_comb.avi',imnum),'COLORMAP',graymap);
%     for framenum=1:iterations+1
%         aviobj = addframe(aviobj,myIm2Frame(reshape(xk(:,framenum),[rows,cols]),graymap));
%     end
%     aviobj = close(aviobj);
    
    % Save raw data.
    save(sprintf('results2/%.2d_data.mat', imnum), 'rgrad', 'rcg', 'rv', 'rcombined');
    
    % Make a figure.
    fig = figure(1);
    plot(t,rgrad(:,1),'b-', t,rgrad(:,end),'c-', t,rcg(:,1),'r-', t,rcg(:,end),'m-', t,rcombined,'g-', 1:length(rv),rv,'k-');
    ylim([0,0.1]);
    xlabel('Iterations');
    ylabel('Ground Truth Error');
    legend('Gradient', 'MG Gradient', 'CG', 'MG CG', 'Combined');
    %print(fig,sprintf('results2/%.2d_plot.png',imnum),'-dpng');
    saveas(fig, sprintf('results2/%.2d_plot.fig',imnum));
    %close all;
    
    % Clean up shit.
    clear A c b M xk;
end