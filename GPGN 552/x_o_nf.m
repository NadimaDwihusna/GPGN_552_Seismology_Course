% Function for Near Field Source Wavelet
function [x_o_nf] = x_o_nf(R,alpha,beta,w,T,t)

x_o_nf = zeros(1,length(t));
for i = 1:length(t)
    if t(i) < (R/alpha)
        x_o_nf(i) = 0;
    elseif (R/alpha) <= t(i) && t(i) <= (R/alpha + T)
        upper_bound(i) = ((R/alpha*cos(w*(t(i)-R/alpha)))/w) + ((sin(w*(t(i))-R/alpha))/(w^2));
        lower_bound(i) = ((R/alpha*cos(w*(R/alpha-R/alpha)))/w) + ((sin(w*(R/alpha)-R/alpha))/(w^2));
        x_o_nf(i) = -(upper_bound(i) - lower_bound(i));
    elseif (R/beta) <= t(i)&& t(i) <= ((R/beta)+T)
        upper_bound(i) = ((R/beta*cos(w*(t(i)-R/beta)))/w) + ((sin(w*(t(i))-R/beta))/(w^2));
        lower_bound(i) = ((R/beta*cos(w*(R/beta-R/beta)))/w) + ((sin(w*(R/beta)-R/beta))/(w^2));
        x_o_nf(i) = (upper_bound(i) - lower_bound(i));
    else
        x_o_nf(i) = 0;
    end
end