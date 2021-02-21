Cancer_type = 'BRCA';

% Input time series gene expression data of cell differentiation
expN = csvread(['Data/expN_',Cancer_type,'.csv'],1,1);
expT = csvread(['Data/expT_',Cancer_type,'.csv'],1,1);
th = 12.5;
% Start parallel pool
poolobj = gcp('nocreate');
delete(poolobj);
parpool(25);

% Calculating resilience and status parameters for each sample
[x_effN,beta_effN,x_effT,beta_effT] = Parameter_calculate(expN,expT,th);
% Saving results
fid = fopen(['Data/Intermediate/Parameters_Results/x_effN_',Cancer_tyhpe,'.csv'],'wt');
fprintf(fid,'%d\n',x_effN);
fclose(fid);
fid = fopen(['Data/Intermediate/Parameters_Results/beta_effN_',Cancer_tyhpe,'.csv'],'wt');
fprintf(fid,'%d\n',beta_effN);
fclose(fid);
fid = fopen(['Data/Intermediate/Parameters_Results/x_effT_',Cancer_tyhpe,'.csv'],'wt');
fprintf(fid,'%d\n',x_effT);
fclose(fid);
fid = fopen(['Data/Intermediate/Parameters_Results/beta_effT_',Cancer_tyhpe,'.csv'],'wt');
fprintf(fid,'%d\n',beta_effT);
fclose(fid);

% Estimating system's potential landscape
[minlocs,maxlocs] = Estimate_landscape(x_effN,beta_effN,x_effT,beta_effT,0.2);
% Saving result
fid = fopen(['Data/Intermediate/Landscape_Results/minlocs_',Cancer_tyhpe,'.csv'],'wt');
fprintf(fid,'%d,%d\n',minlocs');
fclose(fid);
fid = fopen(['Data/Intermediate/Landscape_Results/maxlocs_',Cancer_tyhpe,'.csv'],'wt');
fprintf(fid,'%d,%d\n',maxlocs');
fclose(fid);

% Calculating resilience centrality of each gene for every samples
[resN_beta,resT_beta] = Resilience_centrality(expN,expT,th);
%Saving result
csvwrite(['Data/Intermediate/Resilience_Centrality_Results/resN_beta_',Cancer_tyhpe,'.csv'],resN_beta)
csvwrite(['Data/Intermediate/Resilience_Centrality_Results/resT_beta_',Cancer_tyhpe,'.csv'],resT_beta)
