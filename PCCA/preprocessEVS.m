function EVS=preprocessEVS(EVS,lambda)
% separates real and imaginary parts of complex eigenvectors

k=size(EVS,2);  %number of columns in EVS
i=1;
while (i<size(EVS,2))
    if (abs(imag(lambda(i)))>eps)
        EVS(:,i)=real(EVS(:,i));
        EVS(:,i+1)=imag(EVS(:,i+1));
        i=i+2;
    else
        i=i+1;
    end
end
if (i<=k)
    if abs(imag(EVS(:,k)))>eps
        disp('reduce kmax by 1; otherwise split of complex eigenpair');
        kmax=kmax-1;
        lambda=lambda(1:end-1);
        EVS(:,end)=[];
    end
end


end
