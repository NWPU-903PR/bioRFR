function [minlocs,maxlocs] = Estimate_landscape(N_x,N_beta,T_x,T_beta,local_th)

%Input:
%   N_x (T_x):A P*1 (Q*1) vector which contains the status parameters of 
%             each reference sample (other time points sample)
%   N_beta (T_beta):A P*1 (Q*1) vector which contains the resilience parameters of 
%             each reference sample (other time points sample)
%Output:
%   minlocs:The minimal values on landscape
%   maxlocs:The maximal values on landscape

    beta = log([N_beta;T_beta]);
    x = log([N_x;T_x]);
    gridx1 = linspace(min(x)-0.5,max(x)+0.5,100);
    gridx2 = linspace(min(beta)-0.5,max(beta)+0.5,100);
    [x1,x2] = meshgrid(gridx1, gridx2);
    x1 = x1(:);
    x2 = x2(:);
    xi = [x1 x2];
    [f,xi] = ksdensity([x,beta],xi);
    [f_beta,~] = ksdensity(beta,gridx2);
    f_beta = repmat(f_beta',100,1);
    f_beta(f_beta < 0.3) = 0.3;
    f = f ./ f_beta;
    f = -log(f+0.01)./2;
    [xq,yq,z] = computeGrid(xi(:,2),xi(:,1),f,500);

    minlocs = cell(500,1);
    for i = 1:500
        [~,locs,~,p] = findpeaks(-z(:,i));
        minlocs{i} = locs(p >= local_th);
    end
    minlocs_temp = [];
    for i = 1:500
        if isempty(minlocs{i})
            continue;
        else
            minlocs_temp = [minlocs_temp;[minlocs{i},repmat(i,length(minlocs{i}),1)]];
        end
    end
    minlocs = xq(sub2ind(size(xq),minlocs_temp(:,1),minlocs_temp(:,2)));
    minlocs = [minlocs,yq(sub2ind(size(yq),minlocs_temp(:,1),minlocs_temp(:,2)))];

    maxlocs = cell(500,1);
    for i = 1:500
        [~,locs,~,p] = findpeaks(z(:,i));
        maxlocs{i} = locs(p >= local_th);
    end
    maxlocs_temp = [];
    for i = 1:500
        if isempty(maxlocs{i})
            continue;
        else
            maxlocs_temp = [maxlocs_temp;[maxlocs{i},repmat(i,length(maxlocs{i}),1)]];
        end
    end
    maxlocs = xq(sub2ind(size(xq),maxlocs_temp(:,1),maxlocs_temp(:,2)));
    maxlocs = [maxlocs,yq(sub2ind(size(yq),maxlocs_temp(:,1),maxlocs_temp(:,2)))];
end

