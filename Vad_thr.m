function threshold = Vad_thr(x)

    
    Lf = 1024;%Length of frame
    Lo = Lf/2;%Length of overlap
    
    Nf = length(x)/Lf; %Num of frames

    for i = 1: Nf
        x_tmp = x((i-1)*Lo+(1:Lf));
        
        P(i) = 1/Lf * sum(x_tmp.^2);
    end

    threshold = min(P) * 1.24;   
% threshold = min(P) * 1.32;   
end