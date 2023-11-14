function [Pc,chi,A,ovec]=computeMSM(P,nc,target,centers,Nx,Ny)

    

    %sum(varmat,2)
    % perform clustering by PCCA+
    
    kmin=nc;                 % minimum number of clusters
    kmax=nc;                 % maximum number of clusters
    flag=0;                 % flag=1: full optimization, flag=0: unscaled initial guess
    solver='nelder-mead';   % solver for unconstrained optimization problem
    %no display of eigenvalue solver
    opts.disp=0;
    %save('P.mat','P')
    [EVS,lambda,kmax]=compute_subspace(P,kmax,target,opts);
%     if max(max(abs(imag(EVS))))>0
%         warning('complex eigenvectiors')
%     end
    
    lambda
    %EVS=preprocessEVS(EVS,diag(la));


    %density for normalization of eigenvectors
    %pi=ones(n,1)/n;
    [ovec,~]=eigs(P',1,1+10*eps);
    %idx=find(abs(ovec)<eps);
    %ovec(idx)=0;
    [~,idx]=find(abs(ovec)==max(abs(ovec)));
    ovec=ovec*sign(ovec(idx));
    %ovec=max(ovec,eps);
    ovec=max(ovec,0);
    ovec=ovec/sum(ovec);

%     n=size(P,1);
%     ovec=ones(n,1);
%     ovec=ovec/sum(ovec);

    %compute A such that chi=EVS*A by PCCA+
    [chi,A,~,~]=pcca(EVS,ovec,kmin,kmax,flag,solver);
    idx=[1:nc];
    
    cluster_weights=chi'*ovec

    %number of clusters:
    %nc=size(chi,2);

    disp (' ')
    disp ('Coarse-grainded transition matrix:') 


    %Pc=inv(chi'*diag(ovec)*chi)*chi'*diag(ovec)*P*chi %matrix close to singular
    Pc=inv(A)*diag(lambda(1:nc))*A     %ill-conditioned for la close to zero
    %Pc=inv(A'*A)*A'*diag(lambda(1:nc))*A   %matrix close to singular
        
    stop=1;

    %color centroids according to cluster membership
    figure(10)
    for k=1:size(P,1)
        [~,idx]=max(chi(k,:));
        %for cluster=1:nc
        if idx==1
            plot(centers(1,k),centers(2,k),'r*')
            hold on
        end
        if idx==2
            plot(centers(1,k),centers(2,k),'bx')
            hold on
        end
        if idx==3
            plot(centers(1,k),centers(2,k),'ko')
            hold on    
        end
        if idx==4
            plot(centers(1,k),centers(2,k),'mo')
            hold on    
       end
        
    end

    %generate surface plots of chi
    [Xq,Yq] = meshgrid(0:1:Nx,0:1:Ny);
    for k=1:nc
        figure(2+k)
        %plot3(centers(1,:),centers(2,:),chi(:,k),'ro')
        %hold on
        F=scatteredInterpolant(centers(1,:)',centers(2,:)',chi(:,k));
        F.Method = 'natural';
        vq=F(Xq,Yq);
        %mesh(Xq,Yq,vq)
        contour(Xq,Yq,vq)
    end

    %plot partial densities
    n=size(chi,1);
    chic=zeros(n,nc);
    [~,idx]=max(chi,[],2);
    for i=1:n
        chic(i,idx(i))=1;
    end
    for k=1:nc
        cluster=find(chic(:,k)==1);
        Pk=P(cluster,cluster);
        [pik,~]=eigs(Pk',1,1.0001);
        %pik=max(0,pik);
        pik=pik/sum(pik);
        pik_full=zeros(n,1);
        pik_full(cluster)=pik;
        figure(20+k)
        %plot3(centers(1,:),centers(2,:),chi(:,k),'ro')
        %hold on
        F=scatteredInterpolant(centers(1,:)',centers(2,:)',pik_full);
        F.Method = 'natural';
        vq=F(Xq,Yq);
        %mesh(Xq,Yq,vq)
        contourf(Xq,Yq,vq)
    end

end




