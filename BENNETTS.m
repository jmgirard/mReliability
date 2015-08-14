function [S, p_a, p_e] = BENNETTS(M, type)
% Calculate Bennett's Agreement Coefficient (S)
%   [S, p_a, p_e] = BENNETTS(M, type)
%
%   (c) Jeffrey M Girard, 2015
%   
%   Reference: Bennett, E. M., Alpert, R., & Goldstein, A. C. (1954).
%   Communication through limited response questioning.
%   The Public Opinion Quarterly, 18(3), 303–308.

%% Calculate variables
M(all(~isfinite(M),2),:) = [];
[n,~] = size(M);
x = unique(M(:));
x(~isfinite(x)) = [];
q = length(x);

%% Calculate weights
w = nan(q,q);
for k = 1:q
    for l = 1:q
        switch type
            case 'nominal'
                w = eye(q);
            case 'ordinal'
                if k==l
                    w(k,l) = 1;
                else
                    M_kl = nchoosek((max(k,l) - min(k,l) + 1),2);
                    M_1q = nchoosek((max(1,q) - min(1,q) + 1),2);
                    w(k,l) = 1 - (M_kl / M_1q);
                end
            case 'linear'
                if k==l
                    w(k,l) = 1;
                else
                    dist = abs(x(k) - x(l));
                    maxdist = max(x) - min(x);
                    w(k,l) = 1 - (dist / maxdist);
                end
            case 'quadratic'
                if k==l
                    w(k,l) = 1;
                else
                    w(k,l) = 1 - (x(k) - x(l))^2 / (max(x) - min(x))^2;
                end
            case 'radical'
                if k==l
                    w(k,l) = 1;
                else
                    w(k,l) = 1 - sqrt(abs(x(k) - x(l))) / sqrt(max(x) - min(x));
                end
            case 'ratio'
                w(k,l) = 1 - (((x(k) - x(l)) / (x(k) + x(l)))^2) / (((max(x) - min(x)) / (max(x) + min(x)))^2);
            otherwise
                error('Type must be nominal, ordinal, linear, quadratic, radical, or ratio');
        end
    end
end

%% Calculate percent agreement for each item and overall
p_a_i = zeros(n,1);
for i = 1:n
    r_i = sum(isfinite(M(i,:)));
    if r_i >= 2
        for k = 1:q
            r_ik = sum(M(i,:)==x(k));
            rstar_ik = 0;
            for l = 1:q
                w_kl = w(k,l);
                r_il = sum(M(i,:)==x(l));
                rstar_ik = rstar_ik + (w_kl * r_il);
            end
            p_a_i(i) = p_a_i(i) + (r_ik * (rstar_ik - 1)) / (r_i * (r_i - 1));
        end
    end
end
p_a = sum(p_a_i) / sum(sum(isfinite(M),2)>=2);

%% Calculate percent chance agreement for each item and overall
p_e_i = zeros(n,1);
for i = 1:n
    T_w = sum(sum(w));
    p_e_i(i) = T_w / (q ^ 2);
end
p_e = mean(p_e_i);

%% Calculate S point estimate
S = (p_a - p_e) / (1 - p_e);

end