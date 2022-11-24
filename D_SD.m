function z0 = D_SD(blta, W, Nnode, Ws)

     
    M = size(W,2);% num of all nodes
    u = ones(M,1);
    I = eye(M);
    Wc = W - I;
              
    L_eq = 0;  L_ineq = zeros(M,1); sa = zeros(M,1);  
    for itr = 1:100
        %%%%%%%%%%%%%%%%
%             sa = blta + L_eq*u  - (W-I)'*L_ineq;%   sa = blta_i + L_eq - (W-I)'_i*L_ineq_i
%             s = sign(sa);
        for n = 1:M
        sa(n,1) = blta(n,1) + L_eq - W(n,:)*L_ineq + L_ineq(n,1);
        end
        s = sign(sa);
        %%%%%%%%%%%%%%%%%            
         % consensus iteration
        us = s;
        for k = 1:100
            us = Ws*us;
        end
        us = us * M;

        dL_eq = 2*Nnode - M - round(mean(us));%since us is an integer, round is necessary
        for n = 1:M
            dL_ineq(n,1) = Wc(n,:)*u + W(n,:)*s - s(n,1);
        end

        if norm(dL_eq)~=0  
            dL_eq = dL_eq/norm(dL_eq); 
        end
        if norm(dL_ineq)~=0  
            dL_ineq = dL_ineq/norm(dL_ineq); 
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%
        for n = 1:M
            W_dl(n,1) = W(n,:)*dL_ineq - dL_ineq(n,1);
        end
        % consensus iteration
        W_DL = W_dl;
        for k = 1:100
            W_DL = Ws*W_DL;
        end
        W_DL = mean(W_DL * M);

        for n=1:M
            q(n,1) = sa(n,1);
            d(n,1) = dL_eq - W_dl(n,1);                
        end
        b = dL_eq*(2*Nnode-M) + W_DL;
%             b = ( dL_eq*(2*Nnode-M) + dL_ineq'*(W-I)*u );
        %%%%%%%%%%%%%%%%%%%%
%             f(itr) = norm(sa,1) + b;
        % consensus for q and d are omitted, at all nodes:
        % eliminate zero elements
        d_nz = d;
        q_nz = q;


        ind = d~=0;           
        if sum(ind)~=0
            d_nz(d==0) = [];
            q_nz(d==0) = [];

           [val_qd, idx_d] = sort( q_nz./d_nz);
           d_nz_srt = d_nz(idx_d);

            val_obj_i = zeros(1:val_qd,1);
            for i = 1:length(val_qd)
                val_obj_i(i) = abs( b+sum(abs(d_nz_srt))-2*sum(abs(d_nz_srt(i:length(d_nz_srt)))) );               
            end
    %         val_obj_tt = abs( b-sum(abs(d_nz_srt))+2*sum((2*cumsum(abs(d_nz_srt)))) ); 
            [ ~ ,idx_i] = min(val_obj_i);
            galma =   eps*sign(val_qd(idx_i))-val_qd(idx_i);  
        else
            galma = 0;
        end

%             g(itr) = galma;
        galma = max(galma,0);


        L_eq = L_eq + galma*dL_eq;
        L_ineq = L_ineq + min( galma*dL_ineq, 0);

    end

     sa = blta + L_eq*u  - (W-I)'*L_ineq;

    [~,index]=sort(sa);
    z0 = zeros(M,1);
    for i = M-Nnode+1:M
        z0(index(i)) = 1;
    end
         
end