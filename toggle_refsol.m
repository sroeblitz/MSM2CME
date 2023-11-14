function chi=toggle_refsol()

    tau=5000;
    nc=2;
    target=10*eps;
    N=300;

    %rate matrix with column sum zero
    %Q=toggle_assembleQ(N);
    load('Q-toggle-exact.mat')
    assert(max(abs(sum(Q,1)))<1e-12,'Columns sums of Q are not zero!')
    
    eigs(Q,3,'smallestabs')

    n=size(Q,1);
    opts.disp=0;
    [EVS,lambda,kmax]=compute_subspace(Q',nc,target,opts);
    [ovec,~]=eigs(Q,1,target);
    idx=find(abs(ovec)==max(abs(ovec)));
    ovec=ovec*sign(ovec(idx));
    %ovec=max(ovec,eps);
    ovec=abs(ovec);
    ovec=ovec/sum(ovec);
    logvec=-log(ovec);

    %plot the stationary density
    figure(10)
    %h=surf(reshape(ovec,N+1,N+1));
    %set(h,'LineStyle','none');
    h=contour(reshape(ovec,N+1,N+1));
    colorbar()
    %set(gca,'ColorScale','log')
    %title(sprintf('Global stationary density'));
    xlabel('x_1')
    ylabel('x_2')
%     %fname= sprintf('Fig_%d_Case_%d.png', 10+nc+1, Case);
%     %saveas(gcf, fullfile(filepath,fname));

    [chi,A,~,~]=pcca(EVS,ovec,nc,nc,0,'nelder-mead');
    cluster_weights=chi'*ovec

    %Qc=inv(chi'*diag(ovec)*chi)*chi'*diag(ovec)*Q*chi
    F1=zeros(2,2);
    F2=zeros(2,2);
    for i=1:2
        for j=1:2
            F1(i,j)=sum(chi(:,i).*ovec.*chi(:,j));
            F2(i,j)=(chi(:,i).*ovec)'*Q*chi(:,j);
        end
    end
    Qc=inv(F1)*F2
    

    %expm(tau*Qc)
    Pc=inv(A)*diag(exp(tau*lambda(1:nc)))*A  

    Qc2=inv(A)*diag(lambda(1:nc))*A

    %figure(1)
    %plot(eigs(Q',5,target),'*')

    % Plot the membership functions in 2D
    for k=1:nc
        chi_2D=reshape(chi(:,k),N+1,N+1);
        figure(k+1)
        %contour(-log(chi1_2D))
        h=contour(chi_2D);
        %h=contourf(chi_2D);
        colorbar()
        %set(h,'LineStyle','none')
        %title(sprintf('Membership functions, chi_%d',k))
        xlabel('x_1')
        ylabel('x_2')
        %fname= sprintf('Fig_%d_Case_%d.png', k+1, Case);
        %saveas(h, fullfile(filepath,fname));
    end

    

    % Make the fuzzy clustering crisp
    chi_crisp=zeros(n,nc);
    for i=1:n
        [val,idx]=max(chi(i,:));
        chi_crisp(i,idx)=1;
    end

    pik_2D_sum=zeros(N+1,N+1);
    for k=1:nc
            [pik_full]=pikFullfnc(k,chi_crisp,Q,target,n);
            %findSet_2(k,pik_full,filepath,Case,N)
            %------------------------------------------------------------------
            %pik_full=chi(:,k).*pi;
            pik_2D=reshape(pik_full,N+1,N+1);
            pik_2D_sum=pik_2D_sum+pik_2D;
            figure(10+k)
            h=surf(pik_2D);
            set(h,'LineStyle','none');
            colorbar
            %h=contour(pik_2D);
            %colorbar
            title(sprintf('Partial stationary densities of cluster %d', k))
            %fname= sprintf('Fig_%d_Case_%d.png', 10+k, Case);
            %saveas(gcf, fullfile(filepath,fname));
            %saveas(gcf,sprintf('Fig_%d_Case_%d.png', 10+k, Case))
    end
    figure(10+nc+1)
    h=surf(pik_2D_sum);
    set(h,'LineStyle','none');
    colorbar
    set(gca,'ColorScale','log')
    title(sprintf('Sum of local stationary densities'));
    xlabel('x_1')
    ylabel('x_2') 
    
end