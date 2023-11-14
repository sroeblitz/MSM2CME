function [Pmats,chimats,eigvals,weights]=computeMSM_with_error_dirichlet(alphaMat,nc,target,numMSM,centers,Nx,Ny,chiref)

    nrows=size(alphaMat,1);
    
    k=0;
    while k<numMSM
    
        k=k+1
        P=zeros(nrows,nrows);
        for r=1:nrows
            P(r,:)=dirichlet_sample(alphaMat(r,:),1);
        end
        
        % perform clustering by PCCA+
        kmin=nc;                 % minimum number of clusters
        kmax=nc;                 % maximum number of clusters
        flag=0;                 % flag=1: full optimization, flag=0: unscaled initial guess
        solver='nelder-mead';   % solver for unconstrained optimization problem
        opts.disp=0;            %no display of eigenvalue solver
        [EVS,lambda,kmax]=compute_subspace(P,kmax,target,opts);
    %     if max(max(abs(imag(EVS))))>0
    %         warning('complex eigenvectiors')
    %     end
        
        lambda
           
        %density for normalization of eigenvectors
        [ovec,~]=eigs(P',1,1+10*eps);
        [~,idx]=find(abs(ovec)==max(abs(ovec)));
        ovec=ovec*sign(ovec(idx));
        %ovec=max(ovec,eps);
        ovec=max(ovec,0);
        ovec=ovec/sum(ovec);
    
        %compute A such that chi=EVS*A by PCCA+
        [chi,A,~,~]=pcca(EVS,ovec,kmin,kmax,flag,solver);
        
        
        %make chi coarse
        chic=zeros(nrows,nc);
        [~,idx]=max(chi,[],2);
        for i=1:nrows
            chic(i,idx(i))=1;
        end    
        %ensure that the order of the clusters is the same as in the
        %reference solution
        idx=reorder_chic(chic,chiref);
        chi=chi(:,idx);
        %lambda=lambda(idx);
        A=A(idx,:);
 
    
%         %generate surface plots of chi
%         [Xq,Yq] = meshgrid(0:1:Nx,0:1:Ny);
%         for i=1:nc
%             figure(2+i)
%             %plot3(centers(1,:),centers(2,:),chi(:,i),'ro')
%             %hold on
%             F=scatteredInterpolant(centers(1,:)',centers(2,:)',chi(:,i));
%             F.Method = 'natural';
%             vq=F(Xq,Yq);
%             %mesh(Xq,Yq,vq)
%             contour(Xq,Yq,vq)
%             title('chi%d',i)
%         end
    
    
        cluster_weights=chi'*ovec
    
           
        disp (' ')
        disp ('Coarse-grainded transition matrix:') 
    
    
        %Pc=inv(chi'*diag(ovec)*chi)*chi'*diag(ovec)*P*chi %matrix close to singular
        Pc=inv(A)*diag(lambda(1:nc))*A    
        %Pc=inv(A'*A)*A'*diag(lambda(1:nc))*A   %matrix close to singular
            
        Pmats(k).pc=Pc;
        chimats(k).chi=chi;
        eigvals(k).eig=lambda(1:nc);
        weights(k).w=cluster_weights;
       
    
      
        if(max(max(abs(Pc))>1))
            k=max(0,k-1);
        end
        if length(unique(idx))<nc
            warning('reordering of membership vectors went wrong')
            k=max(0,k-1);
        end
    
    end


end

