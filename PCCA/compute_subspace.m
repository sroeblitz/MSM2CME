function [EVS,lambda,kmax]=compute_subspace(M,kmax,target,opts)
%compute invariant subspace of dimension kmax for matrix M corresponding to
%eigenvalues closest to target

% Compute eigenvectors and sort them according to distance from 1
[EVS,la]=eigs(M,kmax,target,opts);
la=diag(la);
[~,index]=sort(abs(target-la),'ascend');
EVS=EVS(:,index);
lambda=la(index);

% extract real and imaginary parts of eigenvectors and eigenvalues
k=size(EVS,2);
i=1;
while (i<=k)
    if (abs(imag(lambda(i)))~=0)
        EVS(:,i)=real(EVS(:,i));
        %lambda(i)=real(lamba(i));
        if i<k
            EVS(:,i+1)=imag(EVS(:,i+1));
            %lambda(i+1)=imag(lambda(i+1));
            i=i+2;
        else
            disp('reduce kmax by 1; otherwise split of complex eigenpair');
            kmax=kmax-1;
            lambda=lambda(1:end-1);
            EVS(:,end)=[];
            return;
        end
    else
        i=i+1;
    end
end
