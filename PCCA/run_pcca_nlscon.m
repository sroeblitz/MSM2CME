function [x_estim, f]=run_pcca_nlscon(alpha,EVS,k)

    problem='problem_pcca_nlscon';
    
    
    %additional function parameter
    par.evs=EVS;
    par.k=k;

    
    %initial guess 
    x_guess=alpha;
    
    
    %model function with initial parameter guess
    fstart = feval(problem, x_guess, '', par) ;
    norm(fstart);
    

    % run nlscon
    x_estim = main_nlscon(x_guess, par, problem);
    
    
    %model function with optimized parameters
    f = feval(problem, x_estim, '', par) ;   
    norm(f);
    
    

    
  

    