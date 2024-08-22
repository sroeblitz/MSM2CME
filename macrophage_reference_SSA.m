%SSA sampling for model B (macrophage polarization)
%generate N_rows SSA trajectories of lengths #cells * N_vert * tau
%and compute some statistics like the weights of the clusters

Nx=100;
Ny=100;

N_ssa=0;
tau=10;
N_rows=10;
N_vert=200;
num_cells=117;
num_trajectories = num_cells * N_vert; 
tf=num_trajectories*tau;


c1=zeros(1,N_rows);
c2=zeros(1,N_rows);
c3=zeros(1,N_rows);
out_count=zeros(1,N_rows);
num_hor=zeros(1,N_rows);

cell=0;
voronoi_table=[];
N0=[];
for n=1:N_rows
    n
    x0=round(rand(2,1).*[Nx;Ny])
    X_total=[];
    %break sampling down in smaller parts to save some intermediate results
    for k=1:1
        k
        %tf=10^7; 
        save_name = sprintf('SSA_results/macrophage_reference_SSA_tf%d_%d.mat', k*tf,n);
        tic
        [X,tvec,t_exit]=macrophage_SSA_new(x0,N_ssa,tf,cell,voronoi_table,N0);
        toc
        save(save_name,'X')
        %load(save_name);
        x0=X(:,end);
        X_total=[X_total,X];
    end


    count_matrix=sparse(Nx+1,Ny+1);
    
    num_hor(n)=size(X_total,2);
    %tic
    for k=1:num_hor(n)
        i=X_total(1,k);
        j=X_total(2,k);
        if i<(Nx+1) && j<(Ny+1)
            count_matrix(i+1,j+1)=count_matrix(i+1,j+1)+1;
        else
            out_count(n)=out_count(n)+1;
        end
    end
    %toc
    
    %compute statistical weights
    %cluster 1:
    for i=13:100
        for j=0:10
            c1(n)=c1(n)+count_matrix(i+1,j+1);
        end
    end
    %cluster 2:
    for i=0:100
        for j=11:100
            c2(n)=c2(n)+count_matrix(i+1,j+1);
        end
    end
    %cluster 3:
    for i=0:12
        for j=0:10
            c3(n)=c3(n)+count_matrix(i+1,j+1);
        end
    end

    %plot density
    contour(count_matrix')

    %save_name = sprintf('SSA_results/macrophage_reference_SSA_results');
    %save(save_name,'c1','c2','c3','num_hor','out_count')   

end

c1_rel=c1./(c1+c2+c3)
c2_rel=c2./(c1+c2+c3)
c3_rel=c3./(c1+c2+c3)

figure(1)
boxplot([c1_rel',c2_rel',c3_rel'])


