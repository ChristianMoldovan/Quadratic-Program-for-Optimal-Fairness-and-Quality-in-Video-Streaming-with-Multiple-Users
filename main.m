gurobi_home=getenv('GUROBI_HOME');
run([gurobi_home '\matlab\gurobi_setup']); % set Gurobi path
runID = 1586;
objective = 2;
numberOfUsers = 3;
userdelay = [1;50;100];
traffictrace = ones(1,1050).*100000;
content = 0;
load segmentsizes.mat;
segmentsize(:,:,1) = segmentsizes{7};
for u=2:numberOfUsers
    segmentsize(:,:,u) = segmentsize(:,:,1);
end
run_length = 581.5500;
segmentlength = 5;
alpha = 1;
beta = 0.5;
gamma = 0.5;


stats = runFairDash(runID, objective, numberOfUsers, userdelay, traffictrace+1,[],350,1,content,segmentsize,run_length,0,segmentlength,alpha,beta,gamma);
