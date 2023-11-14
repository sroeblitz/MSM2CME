function [X,tvec,t_exit]=toggle_SSA(x0,N_ssa,tf,cell,voronoi_table,N)
    
    t_exit=Inf;
    
    c1=3e+3;    %c1
    c2=1.1e+4;  %c2
    c3=1e-3;    %c3
    c4=3e+3;    %c4
    c5=1.1e+4;  %c5
    c6=1e-3;    %c6
    beta=2;       %beta
    gamma=2;       %gamma

    

    %stoichiometric vectors
    Nr=[1,0;-1,0;0,1;0,-1]';

    X=[]; tvec=[];
    
    z=x0(1); y=x0(2);
    x=x0;
    X(:,1)=x;
    t=0;
    tvec(1)=t;
    i=1;
    while i<N_ssa || t<tf
       i=i+1;
       alpha=[c1/(c2+y^beta), c3*z, c4/(c5+z^gamma),c6*y];	%propensities
       %check feasibility off all reactions			
       for k=1:4
            if min(x+Nr(:,k))<0 %reaction outside the domain
                alpha(k)=0;
            end
       end
       W=sum(alpha);
        if W==0
           warning('negative copy number')
           break 
        end
        if W>0
            tau=-1/W*log(rand(1));  %when does the next reaction take place?
            r=find(cumsum(alpha)>=W*rand(1),1); %which reaction does take place?
            while alpha(r)==0   %jump over non-feasible ractions
                r=r+1;
            end
            %ensure that trajectory stays within cell
            if cell>0 
                cell_idx=membership(x+Nr(:,r),voronoi_table,N);
                if cell_idx==cell
                    x=x+Nr(:,r);    %update state vector
                else %trajectory has left the cell
                    if isinf(t_exit)    %only report very first exit
                          t_exit=t;
                    end
                end      
                t=t+tau;    %update time
                z=x(1); y=x(2);
                X(:,i)=x;
                tvec(i)=t;
            else
                x=x+Nr(:,r);    %update state vector
                t=t+tau;    %update time
                z=x(1); y=x(2);
                X(:,i)=x;
                tvec(i)=t;
            end
        end
    end 
end
