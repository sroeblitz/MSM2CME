function  visualize_errors2(Pmats,chimats,eigvals,weights,numMSM,centers,chiref,Nx,Ny)
%visualization of errors for model A (toggle switch)


    %reference solution on full grid
    %chifull=toggle_refsol();
    load('toggle_chiref.mat')
    %reference solution at center points
    chiR=Voronoi_chi(centers,Nx,Ny,chifull);
    %reorder chiR according to membership from average model
    idx=reorder_chic(chiR,chiref)
    if length(unique(idx))<2
        warning('reordering of full chi went wrong')
        return;
    end
    chiR=chiR(:,idx);


 
    w1=zeros(1,numMSM);
    w2=zeros(1,numMSM);
    lambda1=zeros(1,numMSM);
    lambda2=zeros(1,numMSM);
    P11=zeros(1,numMSM);
    P22=zeros(1,numMSM);
    chi1=zeros(1,numMSM);
    chi2=zeros(1,numMSM);
    for k=1:numMSM
        P=Pmats(k).pc;
        w=weights(k).w;
        w1(k)=w(1);
        w2(k)=w(2);
        P11(k)=P(1,1);
        P22(k)=P(2,2);
        la=eigvals(k).eig;
        lambda1(k)=la(1);
        lambda2(k)=la(2);
        chi=chimats(k).chi;
% %         chi1(k)=norm(diag(pi)*(chi(:,1)-chiR(:,1)));
% %         chi2(k)=norm(diag(pi)*(chi(:,2)-chiR(:,2)));
% %         chi3(k)=norm(diag(pi)*(chi(:,3)-chiR(:,3)));
%         chi1(k)=norm((chi(:,1)-chiR(:,1)));
%         chi2(k)=norm((chi(:,2)-chiR(:,2)));
    end
    
    figure(11)
    boxplot([w1',w2'])
    ylim([0 1]);
    xlabel('cluster')
    ylabel('weight')
    figure(12)
    %boxplot([-tau./log(lambda2)',-tau./log(lambda3)'])
    boxplot([lambda1',lambda2'])
    xlabel('cluster')
    %ylim([0.95 1]);
    ylabel('eigenvalue')
    figure(13)
    boxplot([P11',P22'])
    ylabel('holding probability')
    xlabel('cluster')
%     figure(14)
%     boxplot([chi1',chi2']);
%     ylabel('error in membership')
%     xlabel('cluster')



end
