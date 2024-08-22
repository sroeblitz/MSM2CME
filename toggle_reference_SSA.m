%SSA sampling for model A (toggle switch)
%generate N_rows SSA trajectories of lengths #cells * N_vert * tau
%and compute some statistics like the weights of the clusters


N_ssa=0;
tau=5000;
N_rows=10;
N_vert=500;
num_cells=60;
%num_trajectories = num_cells * N_vert; 
%tf=num_trajectories*tau;



UL=zeros(1,N_rows);
LR=zeros(1,N_rows);
UL_rel=zeros(1,N_rows);
LR_rel=zeros(1,N_rows);
diagonal=zeros(1,N_rows);
crossings=zeros(1,N_rows);
out_count=zeros(1,N_rows);
num_hor=zeros(1,N_rows);
Pii=zeros(1,N_rows);


cell=0;
voronoi_table=[];
N0=[];
for n=1:N_rows
    n
    x0=[1;1];
    X_total=[];
    %break sampling down in smaller parts to save some intermediate results
    for k=1:15
        k
        tf=10^7; 
        save_name = sprintf('SSA_results/toggle_reference_SSA_tf%d_%d.mat', k*tf,n);
        %tic
        [X,tvec,t_exit]=toggle_SSA_new(x0,N_ssa,tf,cell,voronoi_table,N0);
        %toc
        save(save_name,'X')
        %load(save_name);
        x0=X(:,end);
        X_total=[X_total,X];
    end


    count_matrix=sparse(301,301);
    
    num_hor(n)=size(X_total,2);
    %tic
    for k=1:num_hor(n)
        i=X_total(1,k);
        j=X_total(2,k);
        if i<301 && j<301
            count_matrix(i+1,j+1)=count_matrix(i+1,j+1)+1;
            if k>1 && k<num_hor(n) && (i==j)
                if ((X_total(1,k-1)>X_total(2,k-1)) && X_total(1,k+1)<X_total(2,k+1)) ||...
                        ((X_total(1,k-1)<X_total(2,k-1)) && X_total(1,k+1)>X_total(2,k+1))
                    crossings(n)=crossings(n)+1;
                end
            end
        else
            out_count(n)=out_count(n)+1;
        end
    end
    %toc
    
    %compute statistical weights
    for i=0:300
        for j=i+1:300
            UL(n)=UL(n)+count_matrix(i+1,j+1);
        end
        for j=0:(i-1)
            LR(n)=LR(n)+count_matrix(i+1,j+1);
        end
        diagonal(n)=diagonal(n)+count_matrix(i+1,i+1);
    end
    
    UL_rel(n)=(UL(n)+0.5*diagonal(n))/(UL(n)+LR(n)+diagonal(n))
    LR_rel(n)=1-UL_rel(n)

    %compute holding probabilities
    qij=crossings(n)/2/(15*tf);
    Q=[-qij,qij;qij,-qij];
    P=expm(5000*Q);
    Pii(n)=P(1,1);    
     
    %plot density
    contour(count_matrix')

%     save_name = sprintf('SSA_results/toggle_reference_SSA_results');
%     save(save_name,'UL','num_hor','out_count','crossings')   

end

%save_name = sprintf('SSA_results/toggle_reference_SSA_results');
%save(save_name,'UL','num_hor','out_count','crossings')

figure(1)
boxplot([UL_rel',LR_rel'])
figure(2)
boxplot([Pii',Pii'])


