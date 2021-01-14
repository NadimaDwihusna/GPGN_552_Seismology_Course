% Function for Far Field Source Wavelet
function[x_o]=x_o(R,beta,w,T,t)

x_o = zeros(1,length(t));
for i = 1:length(t)
    if t(i) < (R/beta)
        x_o(i) = 0;
    elseif (R/beta) <= t(i) && t(i) <= ((R/beta)+T)
        x_o(i) = sin(w*(t(i)-(R/beta)));
    else
        x_o(i) = 0;
    end
end
      
        


