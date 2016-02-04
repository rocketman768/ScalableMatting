function [ A, B ] = mardia(I,epsilon,win_size)
  
  if (~exist('epsilon','var'))
    epsilon=1e-7;
  end  
  if (isempty(epsilon))
    epsilon=1e-7;
  end
  if (~exist('win_size','var'))
    win_size=1;
  end     
  if (isempty(win_size))
    win_size=1;
  end     

  neb_size=(win_size*2+1)^2;
  [h,w,c]=size(I);

  A = zeros(h,w);
  B = zeros(h,w);
  len=0;
  for j=1+win_size:w-win_size
    for i=win_size+1:h-win_size
      %if (consts(i,j))
      %  continue
      %end
      winI=I(i-win_size:i+win_size,j-win_size:j+win_size,:);
      winI=reshape(winI,neb_size,c);
      win_mu=mean(winI,1)';
      %win_var=inv(winI'*winI/neb_size-win_mu*win_mu' +epsilon/neb_size*eye(c));
      sigma = winI'*winI/neb_size-win_mu*win_mu' +epsilon/neb_size*eye(c);
      
      winI=winI-repmat(win_mu',neb_size,1);
      tvals=(winI*(sigma\winI'));
      
      % A should be chi-squared with 1/6 * c*(c+1)*(c+2) degrees of freedom
      A(i,j) = sum(sum(tvals.^3))/(6*neb_size);
      % B should be N(0,1)
      B(i,j) = sqrt(neb_size/(8*c*(c+2)))*(sum(diag(tvals).^2)/neb_size - c*(c+2));
      
      len=len+neb_size^2;
    end
  end


end

