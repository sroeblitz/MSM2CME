    Nx=300;
    Ny=300;
    
    nc=2;

%experiment 1: toggle switch model parameters
%problem dimension
    D=2;
    %number of basis functions (Voronoi cells)
    N=20;
    %number of sampling points for horizontal sampling
    N_ssa=500;
    %number of trajectories
    num_traj=5;
    %maximum number of repetitions before R needs to be small
    max_rep=5;
    %maximum deviation of R from 1
    Rthr=0.03;
    %maximum number of cells to be refined
    max_ref=10;
    %maximum refinement level
    max_level=5;    
    %number of vertical sampling points
    N_vert=500; 
    %number of candidate rows for P
    N_rows=10;
    global N_total
    N_total=N;  