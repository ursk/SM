%classical MRF with im2col on images
% profile on; cmrf(1,'test',2); profile off; profile viewer

function mrf(LOADFILE, fname, step)    
% call e.g.  cmrf(1,'test',2) % for 8x8, 25x25, and 100/1000 images


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %define nonlinearities and gradients
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    sum2 = inline('sum(sum(x))');
    
    G   = inline(' -log(cosh(x))');
    g   = inline(' -tanh(x)');
    gp  = inline(' tanh(x).^2-1');
    gpp = inline(' -2* tanh(x) .* (tanh(x).^2-1)');
    
    % Black and Roth nonlinearities
    %implement cauchy here to be consistent
    G   = inline(' -log(1+1/2*x.^2)');
    g   = inline(' -x./(1+1/2*x.^2) ');
    gp  = inline(' -  1./(1+1/2*x.^2)  +    x.^2 ./ (1+1/2*x.^2).^2 ');
    gpp = inline(' 3* x./(1+1/2*x.^2).^2  - 2* x.^3 ./ (1+1/2*x.^2).^3');
    
    
    rev = inline('x( size(x,1) : -1 : 1   ,  :)');
    rev2 = inline('x(:  ,  size(x,2) : -1 : 1 )');
    
    obj=[]; 
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % define model parameters 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    p.n=8; % subfilter size
    p.m=24; %2*p.n-1; %minimum size
    p.K=25; %complete model
    p.T=10000; %number of images
    p.batch=50; %save interval and batch size 
    p.nrm=(p.m-p.n+1)^2 * p.batch * p.K;
    p.method='non'; % 'sto'chastic  or 'non'-stochastic updates
    p.seed=0;

    rand('seed', p.seed); randn('seed', p.seed);
    fname=[fname, '.mat'];

    V = randn(p.n.^2,p.K);
    V = V ./ repmat(sqrt(sum(V.^2)),p.n^2,1);

    %%%%%%%%%%%%% LOAD FILE ??? %%%%%%%%%%%% 
    if LOADFILE==1; eval(['load ' fname]); end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % data now depends on p.parameters
    [X,wM]=data(p.T, p.m, p.m^2, 0); %zero phase and no CGC
    %X=dwM*X; %zero phase whitened data
    X=X(:,find(std(X)>.5));
    p.T=size(X,2);
    X=X(:,randperm(p.T));


     %create a mask to extact elements from the big matrix
    fo=zeros((p.m-p.n+1),(p.m-p.n+1)); fo(1:p.n, 1:p.n)=1; fo=fo(:);
    M=zeros((p.m-p.n+1)^2,(p.m-2*p.n+2)^2);
    for row=1:(p.m-2*p.n+2)^2
      M(:,row)=fo;
      fo=fo([end 1:end-1 ]); %cyclic permute
      if mod(row, p.m-2*p.n+2)==0 % jump a few pixels to next row
        fo=fo([end-p.n+2:end 1:end-p.n+1 ]);
      end
    end
     

    %%%%%%%%%%%%%%%%%%
    % start iteration
    %%%%%%%%%%%%%%%%%%
    iter=1e10;  J=0; dJ=0; 
    p
    for i=1:iter
        if mod(i,1)==0; fprintf('.'); end
        if mod(i,100)==0; fprintf(':'); end
        
        %stocastic gradient: Image and nonliearities
        if p.method=='sto'
          I=X(:,ceil(p.T*rand())); %one image, whitened
        else  
          I=X(:,mod(i,p.batch)+1);  % DEBUG: cycle images 1:10
        end
        Xall=im2col(reshape(I,p.m,p.m), [p.n,p.n] ); % 64x(p.m-p.n+1)^2
      	XV=Xall'*V;     % matrix with all filter resoponses 
    	gALL=g(XV);     % gALL_=rev(gALL); % do all this outside the loop...
    	gpALL=gp(XV);   % gpALL_=rev(gpALL);
    	gppALL=gpp(XV); % size 13x13 (from 20-8+1)
      	
      	%%% First loop for Score function %%%
      	Psi=0; Psip=0; %
      	%tic
      	for k=1:p.K % filters
     	  gMAT{k}=im2col( reshape(gALL(:,k), (p.m-p.n+1),(p.m-p.n+1))  ,[p.n,p.n] );
      	  gpMAT{k}=im2col( reshape(gpALL(:,k), (p.m-p.n+1),(p.m-p.n+1))  ,[p.n,p.n] );
      	  gppMAT{k}=im2col( reshape(gppALL(:,k), (p.m-p.n+1),(p.m-p.n+1))  ,[p.n,p.n] );

      	  Psi   = Psi  + rev(V(:,k))' *  gMAT{k};
      	  Psip  = Psip + rev(V(:,k).^2)' *  gpMAT{k};
    	end
    	%toc
    	PsiK  = repmat(Psi, p.n^2, 1); %matrix
      	PsipK = repmat(Psip, p.n^2, 1); 

      	%%% Second loop for Gradients %%%
      	%tic
      	for k=1:p.K 
       	  dPsi1 = rev(gMAT{k}); 
       	  A0    = gpMAT{k}    .*  repmat( rev(V(:,k)),1,(p.m-2*p.n+2)^2);  %has zeros??
       	  A2=(Xall' * A0).*M; %mask of valid elements
       	  A3=reshape(  A2(find(A2))  ,p.n^2,(p.m-2*p.n+2)^2);% remove most elements
       	  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       	  PsidPsiK(k,:)=sum(  PsiK .* (A3+dPsi1),  2);
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          dPsip1 = 2*repmat(V(:,k),1,(p.m-2*p.n+2)^2) .* rev(gpMAT{k});
          B0    = gppMAT{k}  .*  repmat( rev(V(:,k).^2),1,(p.m-2*p.n+2)^2); 
          B2=(Xall' * B0).*M; %mask of valid elements
          B3=reshape(  B2(find(A2))  ,p.n^2,(p.m-2*p.n+2)^2);% remove most elements
       	  dPsipK(k,:)=sum(B3+dPsip1, 2); %25x144
        end
        %toc
        %% increment the update step with this image
        J   = J   +  sum(.5* Psi.^2 + Psip); %sum out the pixels
       	dJ  = dJ  +  (PsidPsiK) + (dPsipK); %Psi is for all k
     	
     	
    	
    	
    	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    	% write cycle: compute energy and plot stuff
    	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    	if mod(i,p.batch)==0
    	    %%% update step
    	    
    	    V=V - step/p.nrm * dJ';	%AVERAGE update v 
    	    V = V ./ repmat(sqrt(sum(V.^2)),p.n^2,1);
            obj=[obj J/p.nrm]; %AVERAGE error over 10 trials
            J=0; dJ=0; %reset gradients for next batch
    	
        	%new objective: It's J, not E
            fprintf( '(%1.3f)', obj(end) );	
            % create the plots
            subplot(1,2,1); visual(V,6,sqrt(p.K)); colorbar;
        	subplot(1,2,2); plot(obj,'.-'); %ylim([-5 0]);  
        	 drawnow;
        	%save to file
        	eval(['save ' fname ' V obj I wM p']);
    	end
    end
return
