function [samp,numdraw] = HSeta_rs_smallM(M,varargin)

    % enhanced rejection sampler to sample from the density
    % h(x) \propto exp(-M x)/(1+x), x > 0
    % this code is for the case when M = eps < 1

    % Input:
    % eps < 1 determines the density h
    % a and b decide the approximant
    % a < 1 determines the first and second segment of the approximation(default: 1/5)
    % b > 1 determines the third and final segment of the approximation(default: 10)

    % Output
    % samp is a sample from h
    % numdraw is the number of samples (from h_L) to get one sample
    
    % ordinary rejection when M > 1
    % some somersaults to handle the indexing
    n = length(M);
    nid0 = (1:n)';
    M1 = M(M>1);
    nid1 = nid0(M>1);
    nid2 = nid0(M<=1);
    
    x1 = HSeta1(M1);
    xoutall = zeros(n,1);
    xoutall(nid1) = x1;
    
    % if we weren't using the ordinary rejection, the code would just start
    % here
    M = M(M<=1);
    n = length(M);
    

    if ~isempty(varargin)
        if length(varargin)==1
            a = varargin{1};
        elseif length(varargin)==2
            b = varargin{2};
        else
            warning('you passed the wrong number of arguments');
        end
    else
       a = 1/5.*ones(n,1);
       b = 10.*ones(n,1);
    end

    [nu, numat] = approxconstants(M,a,b);
    nuprob = numat./nu;                                          % determine mixture weights
    cumprob = cumsum(nuprob,2);
    cumprob = [zeros(n,1) cumprob];

    A = neglogpost(M,a./M); I = neglogpost(M,1./M); B = neglogpost(M,b./M);
    lambda2 = M.*(I-A)./(1-a);  H2 = 1 - exp(-(I-A));
    lambda3 = M.*(B-I)./(b-1);  H3 = 1 - exp(-(B-I));
    nid = (1:n)';
    
    xout = zeros(n,1);
    
    breakloop = false; numdraw = ones(n,1);

    while ~breakloop

        % draw random uniform
        u = unifrnd(0,1,[n 1]);

        % sample from h_L
        %ind = find(mnrnd(1,nuprob)==1);                          % indicator for which region to sample
        % mnrnd is really slow. not sure why
        ind = sum(1.*(cumprob<u),2);

        uaux = rand(n,1);

        x = zeros(n,1); rat = zeros(n,1);
        sel = (ind==1); Mi = M(sel); 
        x(sel) = (1+a(sel)./Mi).^uaux(sel) - 1;
        rat(sel) = exp(-Mi.*x(sel));

        sel = (ind==2); Mi = M(sel); ai = a(sel); l2i = lambda2(sel);
        x(sel) = ai./Mi -rlog(-uaux(sel).*H2(sel))./l2i;
        xi = x(sel);
        rat(sel) = exp(A(sel) + l2i.*(xi-ai./Mi) - Mi.*xi-rlog(xi));

        sel = (ind==3); Mi = M(sel); l3i = lambda3(sel);
        x(sel) = 1./Mi - rlog(-uaux(sel).*H3(sel))./l3i;
        xi = x(sel);
        rat(sel) = exp(I(sel) + l3i.*(xi-1./Mi) - Mi.*xi-rlog(xi));

        sel = (ind==4); Mi = M(sel); 
        x(sel) = b(sel)./Mi-rlog(-uaux(sel))./Mi;
        xi = x(sel);
        rat(sel) = exp(B(sel) + Mi.*(xi-b(sel)./Mi) - Mi.*xi - rlog(xi));

        done = u<rat;
%         ndone0 = ndone;
%         ndone = ndone + sum(done);
        %xout(ndone0+1:ndone,:) = [nid(done) x(done)];
        xout(nid(done)) = x(done);
        breakloop = all(done);

        numdraw(nid(~done)) = numdraw(nid(~done))+1;

        a = a(~done); b = b(~done); M = M(~done); lambda2 = lambda2(~done);
        lambda3 = lambda3(~done); H2 = H2(~done); H3 = H3(~done); I = I(~done);
        A = A(~done); B = B(~done); cumprob = cumprob(~done,:); n = sum(~done);
        nid = nid(~done);

    end
    %xout = sortrows(xout);
    %samp = xout(:,2);
    xoutall(nid2) = xout;
    samp = xoutall;
end








