%this file produces the figures for the macrophage model (model B)

addpath("PCCA/")
addpath("Dirichlet/")

seed=2;
rng(seed);

Nx=100; Ny=100;
nc=3;
%macrophage_refsol();

%N=100, Nvert=200, tf=10
    load('Results/model1_seed2_alphaMat__N100__N_ssa500__N_vert0200_tf10.mat');
    load('Results/model1_seed2_horSampling__N100__N_ssa500__N_vert0200_tf10.mat');
    tf=10;
    N=100;
% %N=100, Nvert=500, tf=5    
%     load('Results/model1_seed2_alphaMat__N100__N_ssa500__N_vert0500_tf5.mat');
%     load('Results/model1_seed2_horSampling__N100__N_ssa500__N_vert0500_tf5.mat');
%     tf=5;
% %N=20, Nvert=500; tf=10
%     load('Results/model1_seed2_alphaMat__N20__N_ssa500__N_vert0500_tf10.mat');
%     load('Results/model1_seed2_horSampling__N20__N_ssa500__N_vert0500_tf10.mat');
%     nc=2;



    
    target=1+eps;               %target eigenvalue

    cells=find(voronoi_table.log_ans==0);
    cell_centers=voronoi_table.centers(cells,:);
    
    %compute average MSM model (maximum likelihood)
    Pmle=diag(1./sum(Pmle,2))*Pmle;
    [Pc,chi,A,ovec,pi_partial]=computeMSM(Pmle,nc,target,cell_centers',Nx,Ny);

    %visualize partial densities
    figure(101)
    visualize_density(pi_partial(:,1),voronoi_table,N,s,Nx,Ny)
    figure(102)
    visualize_density(pi_partial(:,2),voronoi_table,N,s,Nx,Ny)
    figure(103)
    visualize_density(pi_partial(:,3),voronoi_table,N,s,Nx,Ny)

    %compute MSM with the mean of the Dirichlet distrubution
    %P=computeP(alphaMat);
    %[Pc,chi,A,ovec]=computeMSM(P,nc,target,cell_centers',Nx,Ny);

    %make chi coarse
    chic=zeros(size(chi,1),size(chi,2));
    [~,idx]=max(chi,[],2);
    for i=1:size(chic,1)
        chic(i,idx(i))=1;
    end 

    numMSM=20;
    %use chic as reference for sorting chi correctly
    [Pmats,chimats,eigvals,weights]=computeMSM_with_error_dirichlet(alphaMat,nc,target,numMSM,cell_centers',Nx,Ny,chic);
    %generate boxplots
    visualize_errors3(Pmats,chimats,eigvals,weights,numMSM,cell_centers',chic,Nx,Ny);  
