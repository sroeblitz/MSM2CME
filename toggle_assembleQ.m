function Q=toggle_assembleQ(N)
%% Code assembles the CME matrix for the toggle switch
%% Author: Susanna Roeblitz (susanna.roblitz@uib.no)
 
%%parameters 
    c    = [ 3e3  1.1e4 1e-3  3e3  1.1e4  1e-3  ];
    beta = [ 2 2 ];

%% Assembly CME matrix
    Q=sparse((N+1)^2,(N+1)^2);
    for i=1:(N-1)
        for j=1:(N-1)
            stateIndex=i*(N+1)+j+1;
            state_im1=(i-1)*(N+1)+j+1;
            state_ip1=(i+1)*(N+1)+j+1;
            state_jp1=i*(N+1)+j+2;
            state_jm1=i*(N+1)+j;
            alpha1p=c(1)/(c(2)+j^beta(1));   %(a1*(i-1)^n1/(k1^n1+(i-1)^n1)+S1)*(1./(1+(j/p2)^l2))+b1;
            alpha2p=c(4)/(c(5)+i^beta(2));    %a2*((j-1)^n2/(k2^n2+(j-1)^n2))+S2*(1/(1+(i/p1)^l1))+b2;
            alpha1m=c(1)/(c(2)+j^beta(1));   %(a1*i^n1/(k1^n1+i^n1)+S1)*(1./(1+(j/p2)^l2))+b1;
            alpha2m=c(4)/(c(5)+i^beta(2));    %a2*(j^n2/(k2^n2+j^n2))+S2*(1/(1+(i/p1)^l1))+b2;
            Q(stateIndex,state_im1)=alpha1p;
            Q(stateIndex,state_ip1)=c(3)*(i+1);
            Q(stateIndex,state_jm1)=alpha2p;
            Q(stateIndex,state_jp1)=c(6)*(j+1);
            Q(stateIndex,stateIndex)=-alpha1m-c(3)*i-alpha2m-c(6)*j;
        end
    end
    %deal with 1st and last row and column seperately
    i=0;
    for j=1:(N-1)
        stateIndex=i*(N+1)+j+1;
        state_ip1=(i+1)*(N+1)+j+1;
        state_jm1=i*(N+1)+j;
        state_jp1=i*(N+1)+j+2;
        alpha2p=c(4)/(c(5)+i^beta(2));    %a2*((j-1)^n2/(k2^n2+(j-1)^n2))+S2*(1/(1+(i/p1)^l1))+b2;
        alpha1m=c(1)/(c(2)+j^beta(1));   %(a1*i^n1/(k1^n1+i^n1)+S1)*(1./(1+(j/p2)^l2))+b1;
        alpha2m=c(4)/(c(5)+i^beta(2));    %a2*(j^n2/(k2^n2+j^n2))+S2*(1/(1+(i/p1)^l1))+b2;
        Q(stateIndex,state_ip1)=c(3)*(i+1);
        Q(stateIndex,state_jm1)=alpha2p;
        Q(stateIndex,state_jp1)=c(6)*(j+1);
        Q(stateIndex,stateIndex)=-alpha1m-alpha2m-c(6)*j;
    end
    j=0;
    for i=1:(N-1)
        stateIndex=i*(N+1)+j+1;
        state_im1=(i-1)*(N+1)+j+1;
        state_ip1=(i+1)*(N+1)+j+1;
        state_jp1=i*(N+1)+j+2;
        alpha1p=c(1)/(c(2)+j^beta(1));   %(a1*(i-1)^n1/(k1^n1+(i-1)^n1)+S1)*(1./(1+(j/p2)^l2))+b1;
        alpha1m=c(1)/(c(2)+j^beta(1));   %(a1*i^n1/(k1^n1+i^n1)+S1)*(1./(1+(j/p2)^l2))+b1;
        alpha2m=c(4)/(c(5)+i^beta(2));    %a2*(j^n2/(k2^n2+j^n2))+S2*(1/(1+(i/p1)^l1))+b2;
        Q(stateIndex,state_im1)=alpha1p;
        Q(stateIndex,state_ip1)=c(3)*(i+1);
        Q(stateIndex,state_jp1)=c(6)*(j+1);
        Q(stateIndex,stateIndex)=-alpha1m-c(3)*i-alpha2m;
    end
    i=N;
    for j=1:(N-1)
        stateIndex=i*(N+1)+j+1;
        state_im1=(i-1)*(N+1)+j+1;
        state_jm1=i*(N+1)+j;
        state_jp1=i*(N+1)+j+2;
        alpha1p=c(1)/(c(2)+j^beta(1));   %(a1*(i-1)^n1/(k1^n1+(i-1)^n1)+S1)*(1./(1+(j/p2)^l2))+b1;
        alpha2p=c(4)/(c(5)+i^beta(2));    %a2*((j-1)^n2/(k2^n2+(j-1)^n2))+S2*(1/(1+(i/p1)^l1))+b2;
        alpha1m=c(1)/(c(2)+j^beta(1));   %(a1*i^n1/(k1^n1+i^n1)+S1)*(1./(1+(j/p2)^l2))+b1;
        alpha2m=c(4)/(c(5)+i^beta(2));    %a2*(j^n2/(k2^n2+j^n2))+S2*(1/(1+(i/p1)^l1))+b2;
        Q(stateIndex,state_im1)=alpha1p;
        Q(stateIndex,state_jm1)=alpha2p;
        Q(stateIndex,state_jp1)=c(6)*(j+1);
        Q(stateIndex,stateIndex)=-c(3)*i-alpha2m-c(6)*j;    %-alpha1m
    end
    j=N;
    for i=1:(N-1)
        stateIndex=i*(N+1)+j+1;
        state_im1=(i-1)*(N+1)+j+1;
        state_ip1=(i+1)*(N+1)+j+1;
        state_jm1=i*(N+1)+j;
        alpha1p=c(1)/(c(2)+j^beta(1));   %(a1*(i-1)^n1/(k1^n1+(i-1)^n1)+S1)*(1./(1+(j/p2)^l2))+b1;
        alpha2p=c(4)/(c(5)+i^beta(2));    %a2*((j-1)^n2/(k2^n2+(j-1)^n2))+S2*(1/(1+(i/p1)^l1))+b2;
        alpha1m=c(1)/(c(2)+j^beta(1));   %(a1*i^n1/(k1^n1+i^n1)+S1)*(1./(1+(j/p2)^l2))+b1;
        alpha2m=c(4)/(c(5)+i^beta(2));    %a2*(j^n2/(k2^n2+j^n2))+S2*(1/(1+(i/p1)^l1))+b2;
        Q(stateIndex,state_im1)=alpha1p;
        Q(stateIndex,state_ip1)=c(3)*(i+1);
        Q(stateIndex,state_jm1)=alpha2p;
        Q(stateIndex,stateIndex)=-alpha1m-c(3)*i-c(6)*j;    %-alpha2m
    end
    i=0;
    j=0;
    stateIndex=i*(N+1)+j+1;
    state_ip1=(i+1)*(N+1)+j+1;
    state_jp1=i*(N+1)+j+2;
    alpha1m=c(1)/(c(2)+j^beta(1));   %(a1*i^n1/(k1^n1+i^n1)+S1)*(1./(1+(j/p2)^l2))+b1;
    alpha2m=c(4)/(c(5)+i^beta(2));    %a2*(j^n2/(k2^n2+j^n2))+S2*(1/(1+(i/p1)^l1))+b2;
    Q(stateIndex,state_ip1)=c(3)*(i+1);
    Q(stateIndex,state_jp1)=c(6)*(j+1);
    Q(stateIndex,stateIndex)=-alpha1m-alpha2m;
    %%%%%%%%%%%%%%%%%%%%%%%
    i=0;
    j=N;
    stateIndex=i*(N+1)+j+1;
    state_ip1=(i+1)*(N+1)+j+1;
    state_jm1=i*(N+1)+j;
    alpha2p=c(4)/(c(5)+i^beta(2));    %a2*((j-1)^n2/(k2^n2+(j-1)^n2))+S2*(1/(1+(i/p1)^l1))+b2;
    alpha1m=c(1)/(c(2)+j^beta(1));   %(a1*i^n1/(k1^n1+i^n1)+S1)*(1./(1+(j/p2)^l2))+b1;
    alpha2m=c(4)/(c(5)+i^beta(2));    %a2*(j^n2/(k2^n2+j^n2))+S2*(1/(1+(i/p1)^l1))+b2;
    Q(stateIndex,state_ip1)=c(3)*(i+1);
    Q(stateIndex,state_jm1)=alpha2p;
    Q(stateIndex,stateIndex)=-alpha1m-c(6)*j; %-alpha2m
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    i=N;
    j=0;
    stateIndex=i*(N+1)+j+1;
    state_im1=(i-1)*(N+1)+j+1;
    state_jp1=i*(N+1)+j+2;
    alpha1p=c(1)/(c(2)+j^beta(1));   %(a1*(i-1)^n1/(k1^n1+(i-1)^n1)+S1)*(1./(1+(j/p2)^l2))+b1;
    alpha2m=c(4)/(c(5)+i^beta(2));    %a2*(j^n2/(k2^n2+j^n2))+S2*(1/(1+(i/p1)^l1))+b2;
    Q(stateIndex,state_im1)=alpha1p;
    Q(stateIndex,state_jp1)=c(6)*(j+1);
    Q(stateIndex,stateIndex)=-alpha2m-c(3)*i;     %-alpha1m
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    i=N;
    j=N;
    stateIndex=i*(N+1)+j+1;
    state_im1=(i-1)*(N+1)+j+1;
    state_jm1=i*(N+1)+j;
    alpha1p=c(1)/(c(2)+j^beta(1));   %(a1*(i-1)^n1/(k1^n1+(i-1)^n1)+S1)*(1./(1+(j/p2)^l2))+b1;
    alpha2p=c(4)/(c(5)+i^beta(2));    %a2*((j-1)^n2/(k2^n2+(j-1)^n2))+S2*(1/(1+(i/p1)^l1))+b2;
    Q(stateIndex,state_im1)=alpha1p;
    Q(stateIndex,state_jm1)=alpha2p;
    Q(stateIndex,stateIndex)=-c(3)*i-c(6)*j;    %-alpha1m-alpha2m
end