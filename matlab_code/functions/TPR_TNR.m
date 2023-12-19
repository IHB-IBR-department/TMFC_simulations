function [TPR,TNR,PPV,NPV,FDR] = TPR_TNR(test_matrix,ground_truth)

% ========================================================================
% Ruslan Masharipov, October, 2023
% email: ruslan.s.masharipov@gmail.com
% ========================================================================

T = ground_truth;
P = test_matrix;
P(1:1+size(P,1):end) = 0;
FP = sum(sum((1-T).*P));
FN = sum(sum(T.*(1-P)));
TP = sum(sum(T.*P));
TN = sum(sum((1-T).*(1-P)));
TPR = TP./(TP+FN);
TNR = TN./(TN+FP);
PPV = TP./(TP+FP);
NPV = TN./(TN+FN);
FDR = 1-PPV;

end