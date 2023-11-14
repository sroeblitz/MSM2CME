function P=computeP(alphaMat)

    nrows=size(alphaMat,1);
    P=zeros(nrows,nrows);
    varmat=zeros(nrows,nrows);
    for r=1:nrows
        %r
        alpha=alphaMat(r,:);
        sig=sum(alpha);
        mu=alpha/sig;
        varmat(r,:)=mu.*(1-mu)/(sig+1);
        P(r,:)=mu;
    end

end
