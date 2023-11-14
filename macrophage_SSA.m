function [X,tvec,t_exit]=macrophage_SSA(x0,N_ssa,tf,cell,voronoi_table,N)

    parameters_tristab();

    t_exit=Inf;

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
        %propensities
        alpha=[(a1.*(z.^n1./(z.^n1+k1^n1))+S1).*(1./(1+(y./p2).^l2))+b1,q1*z,a2.*(y.^n2./(y.^n2+k2.^n2))+S2.*(1./(1+(z./p1).^l1))+b2,q2*y];
	%check feasability of each reaction        
	for k=1:4
            if min(x+Nr(:,k))<0 %reaction leads to negative copy numbers
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
            r=find(cumsum(alpha)>=W*rand(1),1);  %which reaction does take place?
            while alpha(r)==0   %jump over non-feasible ractions
                r=r+1;
            end
	    %ensure that trajectory remains within the cell		
            if cell>0 
                cell_idx=membership(x+Nr(:,r),voronoi_table,N);
                if cell_idx==cell
                    x=x+Nr(:,r);    %update state vector
                 else %trajectory has left the cell
                    if isinf(t_exit)
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

