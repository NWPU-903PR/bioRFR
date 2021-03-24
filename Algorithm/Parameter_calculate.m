function [x_effN,beta_effN,x_effT,beta_effT] = Parameter_calculate(expN,expT,th)

% Input:
%     expN:The N*P expression matrix of reference samples, which contains 
%          N genes and P samples
%     expT:The N*Q expression matrix of other time points samples, which contains 
%          N genes and Q samples
%     th:The logarithm threshold of remainding edge number when constructing
%        personalized state transition networks
% Output:
%     x_effN (x_effT):A P*1 (Q*1) vector which contains the calculated
%                     state parameters of each reference sample (other 
%                     time points sample)
%     beta_effN (beta_effT):A P*1 (Q*1) vector which contains the calculated
%                           resilience parameters of each reference sample 
%                           (other time points sample)

    th1 = ceil(exp(th));
    pccref = corr(expN');
    pccref(isnan(pccref)) = 0;
    pccref = (pccref + pccref') ./ 2;
    x_effN = zeros(size(expN,2),1);
    beta_effN = zeros(size(expN,2),1);
    parfor i=1:size(expN,2)
        expN1 = expN;
        expN1(:,i) = [];
        meanN = mean(expN1,2);
        stdN = std(expN1,0,2);
        stdN(stdN == 0) = 1;
        dist = (expN(:,i) - meanN) ./ stdN;
        pccref1 = corr(expN1');
        pccref1(isnan(pccref1)) = 0;
        pccref1 = (pccref1 + pccref1') ./ 2;
        spcc = abs(dist * dist' - pccref1);
        spcc = spcc - diag(diag(spcc));
        [~,~,v_spcc] = find(triu(spcc,1));
        if th > 14
            percentiles = prctile(v_spcc,100-th1/(size(spcc,1)^2-size(spcc,1))*200);
        else
            percentiles = min(maxk(v_spcc, th1));
        end
        v_spcc = 0;
        spcc = spcc .* (spcc >= percentiles);
        [row,col,val] = find(triu(spcc));
        spcc = (spcc ~= 0);
        row_nonzero = sum(spcc,2);
        col_nonzero = sum(spcc,1)';
        edge_degree = row_nonzero(row) + col_nonzero(col) - 2;
        x_effN(i) = edge_degree' * val / sum(edge_degree);
        beta_effN(i) = edge_degree' * edge_degree / sum(edge_degree);
    end

    meanN = mean(expN,2);
    stdN = std(expN,0,2);
    stdN(stdN == 0) = 1;
    x_effT = zeros(size(expT,2),1);
    beta_effT = zeros(size(expT,2),1);
    parfor i=1:size(expT,2)
        dist = (expT(:,i) - meanN) ./ stdN;
        spcc = abs(dist * dist' - pccref);
        spcc = spcc - diag(diag(spcc));
        [~,~,v_spcc] = find(triu(spcc,1));
        if th > 14
            percentiles = prctile(v_spcc,100-th1/(size(spcc,1)^2-size(spcc,1))*200);
        else
            percentiles = min(maxk(v_spcc, th1));
        end
        v_spcc = 0;
        spcc = spcc .* (spcc >= percentiles);
        [row,col,val] = find(triu(spcc,1));
        spcc = (spcc ~= 0);
        row_nonzero = sum(spcc,2);
        col_nonzero = sum(spcc,1)';
        edge_degree = row_nonzero(row) + col_nonzero(col) - 2;
        x_effT(i) = edge_degree' * val / sum(edge_degree);
        beta_effT(i) = edge_degree' * edge_degree / sum(edge_degree);
    end
end

