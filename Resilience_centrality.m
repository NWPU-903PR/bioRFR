function [resN_beta,resT_beta] = Resilience_centrality(expN,expT,th)

% Input:
%     expN:The N*P expression matrix of reference samples, which contains 
%          N genes and P samples
%     expT:The N*Q expression matrix of other time point samples, which contains 
%          N genes and Q samples
%     th:The logarithm threshold of remainding edge number when constructing
%        personalized state transition networks
% Output:
%     resN_beta:A N*P result matrix, where the element (i,j)
%               represents the resilience centrality of gene i in reference
%               sample j
%     resT_beta:A N*Q result matrix, where the element (i,j) 
%               represents the resilience centrality of gene i in other time 
%               points sample j

    th1 = ceil(exp(th));
    pccref = corr(expN');
    pccref(isnan(pccref)) = 0;
    pccref = (pccref + pccref') ./ 2;
    resN_x = zeros(size(expN,1),size(expN,2));
    resN_beta = zeros(size(expN,1),size(expN,2));
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
        nonzero_el = [row,col,val];
        spcc = (spcc ~= 0);
        row_nonzero = sum(spcc,2);
        col_nonzero = sum(spcc,1)';
        edge_degree = row_nonzero(row) + col_nonzero(col) - 2;
        x_effN = edge_degree' * val / sum(edge_degree);
        beta_effN = edge_degree' * edge_degree / sum(edge_degree);
        x_effN_del = zeros(size(expN,1),1);
        beta_effN_del = zeros(size(expN,1),1);
        for j = unique([unique(row);unique(col)])'
            nonzero_el_temp = nonzero_el((row ~= j) & (col ~= j),:);
            edge_degree = row_nonzero(nonzero_el_temp(:,1)) + ...
                col_nonzero(nonzero_el_temp(:,2)) - 2;
            x_effN_del(j) = x_effN - (edge_degree' * nonzero_el_temp(:,3) ...
                / sum(edge_degree));
            beta_effN_del(j) = beta_effN - (edge_degree' * edge_degree / ...
                sum(edge_degree));
        end
        resN_x(:,i) = x_effN_del;
        resN_beta(:,i) = beta_effN_del;
    end

    meanN = mean(expN,2);
    stdN = std(expN,0,2);
    stdN(stdN == 0) = 1;
    resT_x = zeros(size(expT,1),size(expT,2));
    resT_beta = zeros(size(expT,1),size(expT,2));
    parfor i = 1:size(expT,2)
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
        [row,col,val] = find(triu(spcc));
        nonzero_el = [row,col,val];
        spcc = (spcc ~= 0);
        row_nonzero = sum(spcc,2);
        col_nonzero = sum(spcc,1)';
        edge_degree = row_nonzero(row) + col_nonzero(col) - 2;
        x_effT = edge_degree' * val / sum(edge_degree);
        beta_effT = edge_degree' * edge_degree / sum(edge_degree);
        x_effT_del = zeros(size(expT,1),1);
        beta_effT_del = zeros(size(expT,1),1);
        for j = unique([unique(row);unique(col)])'
            nonzero_el_temp = nonzero_el((row ~= j) & (col ~= j),:);
            edge_degree = row_nonzero(nonzero_el_temp(:,1)) + ...
                col_nonzero(nonzero_el_temp(:,2)) - 2;
            x_effT_del(j) = x_effT - (edge_degree' * ...
                nonzero_el_temp(:,3) / sum(edge_degree));
            beta_effT_del(j) = beta_effT - ...
                (edge_degree' * edge_degree / sum(edge_degree));
        end
        resT_x(:,i) = x_effT_del;
        resT_beta(:,i) = beta_effT_del;
    end
end

