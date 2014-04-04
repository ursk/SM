% Demo to do overcomplete ICA with score matching. 

function ica(LOADFILE, fname, step)    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %define nonlinearities and gradients
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    sum2 = inline('sum(sum(x))');
    
    % old nonlinearityes
    G   = inline(' -log(cosh(x))');
    g   = inline(' -tanh(x)');
    gp  = inline(' -1+tanh(x).^2');
    gpp = inline(' -2* tanh(x) .* (tanh(x).^2-1)');
    
    %cauchy
    G   = inline(' -log(1+1/2*x.^2)');
    g   = inline(' -x./(1+1/2*x.^2) ');
    gp  = inline(' -  1./(1+1/2*x.^2)  +    x.^2 ./ (1+1/2*x.^2).^2 ');
    gpp = inline(' 3* x./(1+1/2*x.^2).^2  - 2* x.^3 ./ (1+1/2*x.^2).^3');
    
    rev = inline('x( size(x,1) : -1 : 1   ,  :)');
    rev2 = inline('x(:  ,  size(x,2) : -1 : 1 )');
    
    % nonlinearity used in 2-layer model:
    G   = inline(' -(1+x.^2).^.5'); %use diff() to check
    g   = inline(' -x .* (1+x.^2).^-.5');
    gp  = inline(' x.^2 .* (1+x.^2).^-(3/2)   -   (1+x.^2).^-(1/2)      ');
    gpp = inline(' 3* x .* (1+x.^2).^-(3/2) -  3* x.^3 .* (1+x.^2).^-(5/2)  ');
    

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % define model parameters 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    p.n=16; % window size
    p.K=100; %number of filter 
    p.cols=6;
    p.T=10e3;
    
    [X, wM, dwM]=data(p.T,p.n, p.K, 8); % returns non-witenend data
    
    V=randn(p.K); %one filters: the true filter bank
    V = V ./ repmat(sqrt(sum(V.^2)),p.K,1);
  
    
    %%%%%%%%%%%%%%%%%%
    % start main loop
    %%%%%%%%%%%%%%%%%%
    obj=[];
    
    if LOADFILE==1
    	eval(['load ' fname]);
    end
    tic

    for i=1:1e5
        fprintf('.');

        
        %compute gradients
        y=V'*X; % 10k filter outputs - for 2 filters!
        Psi=( V*g(y) ); % score function: implicit sum over filters 
        Psip= V.^2 * gp(y); % gradient of the score function 
        J=  sum2(  1/2*  Psi.^2    +   Psip   )/p.T; % Score matching objective J sums i and t
        
        % parts of the derivative of the objective wrt. parameter matrix (4 terms)
        dPsi1= X *   (  (V'*Psi).*gp(y)  )'   ; % 
        dPsi2= Psi*g(y)' ; % 
        dPsip1= X * (gpp(y)'  .* repmat(sum(V.^2),p.T,1) ); %
        dPsip2= 2* V.* repmat(sum(gp(y),2), 1, p.K)';

        % update
        V = V - step * (dPsi1+dPsi2+dPsip1+dPsip2)/p.T; % gradient step to reduce the function
        obj= [obj J];
        V = V ./ repmat(sqrt(sum(V.^2)),p.K,1); % constrain the filters to unit norm

        
        % write and plot and stuff
        if mod(i,10)==1
          fprintf('(%2.2f)', obj(end));
          %subplot(1,2,1); plot(obj,'.-'); ylim([-1.5 0])
          %subplot(1,2,2); visual(dwM*V,2,p.cols); colorbar %inv(wM)*
          %drawnow;
          eval(['save ' fname ' wM dwM V obj']);
        end
    
    end
    toc    	
    	
return


