function [matrix_FDR Nsig_FDR pval tval matrix_uncorr001 Nsig_uncorr001] = network_onesample_ttest(matrix,alpha_FDR);

% ========================================================================
% Ruslan Masharipov, October, 2023
% email: ruslan.s.masharipov@gmail.com
% ========================================================================

N = length(matrix);

for roii = 1:N
    for roij = 1:N
        if roii == roij
           pval(roii,roij) = 0;
           tval(roii,roij) = 0;
        else
           [h,pval(roii,roij),ci,stat] = ttest(shiftdim(matrix(roii,roij,:)));
           tval(roii,roij) = stat.tstat;
        end
    end
end

matrix_mask = ones(N,N);
matrix_mask = tril(matrix_mask,-1);
a_mask = matrix_mask(:);

for i = 1:size(matrix,3)
    tmp = matrix(:,:,i);
    a_matrix = tmp(:);
    a_matix_masked(i,:) = a_matrix(find(a_mask)); 
end

[h,p,ci,stats] = ttest(a_matix_masked);
[pFDR,pN] = FDR(p,alpha_FDR);

matrix_FDR = pval<pFDR;
matrix_uncorr001 = pval<0.001;

Nsig_FDR = sum(lower_triangle(double(matrix_FDR)));
Nsig_uncorr001 = sum(lower_triangle(double(matrix_uncorr001)));

end