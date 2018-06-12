% stat=runDashOptimizerRealDataMultipleUsers3checked(runID,objective,numberOfUsers,bw_faktor,params,n,bundle)
% switches as factor of optimization term
function [stat,result,result2,model,model2]=runFairDashm(runID,objective,numberOfUsers,userdelay,traffictrace,params,n,bundle,content,segmentsizes,run_length,sumofstalling,segmentlength,alpha,beta,gamma)
if nargin<12, segmentlength = 5;end
if nargin<11, sumofstalling = 0;end
if nargin<10, run_length = 50;end
if nargin<9,segmentsizes=readvideo();end
if nargin<8,content=false;end
if nargin<7,bundle=1;end
if nargin<6, n=50;end
if nargin<5 || isempty(params),
    params.outputflag = 0;
    params.TimeLimit = 10;
end

VOLUME=1;LAYER=2;SSIM=3;CODEC_BW=4;TIMEHIGHESTLAYER=5;
if nargin<4,bw_faktor=2;end
if nargin<3,numberOfUsers=3;end
if nargin<2,objective=LAYER;end
if nargin<1,
    runID = 0; % which traffic pattern simulation (id=0:29)
end
R=4;
% numberOfUsers = 1;
% net_dl_total = bw_faktor;
%% Parameters
% V(t) : total amount data V(t) received by client during time [0,t]
% Sij  : size of segment i from representation j
% R=3  : number of presentations per segment
% n    : number of segments in the video
% Ts   : time from the start of the download (t=0) until playback starts
% tau  : duration of a segment
% Di   : playback deadlines for individual segments
% xij  : Optimization variable if client downloads segment i from presentation j
% alpha: determines importance of mean resolution and #switches
%
% Di = Ts + (i-1)*tau;
%% set and check parameters
Sij = segmentsizes;
Suij = Sij;
% Suij(:,:,1) = segmentsizes;
% Suij(:,:,2) = 2*segmentsizes;
% traffictrace = 0.1*traffictrace;

n = size(Sij,1); % number of segments

tau=segmentlength;%*bundle; % duration of segment (sec)

V = cumsum(traffictrace);

timeV = (0:length(V)-1)';
%% Define Volume function
Vt = @(xi) interp1(timeV,V,xi); % total amount data V(t) received by client during time [0,t]
Vtrev = @(yi) interp1(V,timeV,yi); %reverse function of V(t): Vtrev(volume)=time

%% check: Ts>=Te
% Te = Vtrev(Sij(1,4)*numberOfUsers); % we account for smallest initial delay to download highest quality
% Te = Vtrev(Sij(1,4)*numberOfUsers)+sumofstalling; % we account for smallest initial delay to download highest quality
% Ts = Te; % no initial delay
Ts = userdelay.*tau;
%% calculate playback deadlines Di
i=1:n;
Di = (Ts + (i-1)*tau)'; % or read it from file
%% MIP model
names = {};
for i=1:n % segments
    for j=1:R+1 % representations
        for u=1:numberOfUsers
            names{j,i,u}=sprintf('x_S%d_R%d_U%d',i,j,u); % Xij(i,j)
        end
    end
end
%% first contraint: download one presenation for each segment
first.A = sparse(n*numberOfUsers,n*R*numberOfUsers+n*numberOfUsers);
for i=1:n
    for u=1:numberOfUsers
        first.A((u-1)*n + i,(u-1)*n*R + (i-1)*R+1+(0:R-1))=1;
    end
end
first.rhs = ones(n*numberOfUsers,1); % >=
first.sense = repmat('=',1,n*numberOfUsers);
%% second constraint: consider the downloadable volume and segment sizes
% for all users l from 1 to U it is
%   for all users u from l to U it is
second.A = sparse(n*numberOfUsers,n*R*numberOfUsers+n*numberOfUsers);
second.rhs = zeros(n*numberOfUsers,1);
for l = 1:numberOfUsers
    %     second.A((l-1)*n + 1:l*n,:)=zeros(n,numberOfUsers*n*R);
    for u = l:numberOfUsers
        timedelay = floor((Ts(u) - Ts(l))/tau);
        if (l-1)*n + 1 + timedelay > l*n % if segment comes after l
            break;
        end
        second.A((l-1)*n + 1 + timedelay,(u-1)*n*R + (1:R))=Suij(1,1:R,u);
        for i=2:n
            if (l-1)*n + i + timedelay > l*n % if segment comes after l
                break;
            end
            second.A((l-1)*n + i + timedelay,(u-1)*n*R + (1:R*(i-1)))=second.A((l-1)*n + i-1 + timedelay,(u-1)*n*R + (1:R*(i-1)));
            second.A((l-1)*n + i + timedelay,(u-1)*n*R + (i-1)*R+1+(0:R-1))=Suij(i,1:R,u);
            %             spy(second.A) % graphical help for debugging
        end
    end
    %     second.rhs((l-1)*n + (1:n))=Vt(Di(:,l))-Vt(Di(1,l));
end
vertVt = Vt(Di)-Vt(Di(1,:)-tau);
second.rhs = vertVt(:);
% second.rhs = Vt(Di)'; % respect playout deadlines
second.sense = repmat('<',1,n*numberOfUsers);
%% third constraint - Fairness is a positive value
% similar quality at every point in time for all players
third.A = sparse(2*n*numberOfUsers,n*R*numberOfUsers+n*numberOfUsers);
for i = 1:n
    for u=1:numberOfUsers
        for j=1:R
            for u_schlange = 1:numberOfUsers
                % mean quality over all users at segment i
                third.A((u-1)*n + i,(u_schlange-1)*n*R + (i-1)*R + j) = j/numberOfUsers;
            end
            % quality of user u at segment i
            third.A((u-1)*n + i,(u-1)*n*R + (i-1)*R + j) = -j;
        end
        % fairness
        third.A((u-1)*n + i,n*R*numberOfUsers+(u-1)*n + i) = 1;
    end
end
for i = 1:n
    for u=1:numberOfUsers
        for j=1:R
            for u_schlange = 1:numberOfUsers
                % mean quality over all users at segment i
                third.A(n*numberOfUsers+(u-1)*n + i,(u_schlange-1)*n*R + (i-1)*R + j) = -j/numberOfUsers;
            end
            % quality of user u at segment i
            third.A(n*numberOfUsers+(u-1)*n + i,(u-1)*n*R + (i-1)*R + j) = j;
        end
        % fairness
        third.A(n*numberOfUsers+(u-1)*n + i,n*R*numberOfUsers+(u-1)*n + i) = 1;
    end
end
third.rhs = zeros(2*n*numberOfUsers,1);
third.sense = repmat('>',1,2*n*numberOfUsers);
%% Model
model.varnames = names(:)';
model.vtype='I'; % binary variables
% constraints
model.A =  vertcat(first.A, second.A, third.A);% sparse([1 2 3; 1 1 0]);
model.rhs=[first.rhs; second.rhs; third.rhs];
model.sense=[first.sense second.sense third.sense];

%% testModel
% model.varnames = names(:)';
% model.vtype='B'; % binary variables
% % constraints
% model.A =  vertcat(first.A);%,second.A);% sparse([1 2 3; 1 1 0]);
% model.rhs=[first.rhs];%; second.rhs];
% model.sense=[first.sense];% second.sense];
%% objective
% SSIMvalues=round(SSIMvalues(1:n,:)*20);
switch objective
    case VOLUME
        tmp=Sij';
        model.obj = repmat(tmp(:),numberOfUsers,1);
    case LAYER
%         beta = 0.8;
        
        LayerVal=repmat(((1:R)-1)./R./n./numberOfUsers,n*numberOfUsers,1)'; % mean layer
        fairness = -1.* ones(n*numberOfUsers,1)./(R-1)./n./numberOfUsers;
        model.obj = [alpha.*LayerVal(:);gamma.*fairness];
    case CODEC_BW
        LayerVal=repmat((1:R).^2,n*numberOfUsers,1)'/n; % mean layer
        %LayerVal=repmat([mean(Sij)/1e6],n*numberOfUsers,1)'/n; % mean layer
        model.obj = LayerVal(:);
    case SSIM
        SSIMvalues=textread('data/SSIM.txt'); % mean layer
        SSIMvalues=SSIMvalues(1:n,:);
        LayerVal=SSIMvalues(1:n,:)';
        model.obj = repmat(LayerVal(:),numberOfUsers,1);
    case TIMEHIGHESTLAYER
        LayerVal=repmat(1000.^((1:R)-1),n*numberOfUsers,1)'/n; % mean layer
        model.obj = LayerVal(:);
    otherwise
        %         warning(sprintf('objective %d not existing',objective));
        return;
end

% minimize switches
% model.Q = sparse(R*n*numberOfUsers+n*numberOfUsers,R*n*numberOfUsers+n*numberOfUsers);
% for u = 1:numberOfUsers
%     for i=1:n-1
%         for j=1:R
%             xij = (u-1)*n*R +(i-1)*R+j;
%             xi1j = (u-1)*n*R +i*R+j;
%             model.Q(xij,xij)=model.Q(xij,xij)+1/2;
%             model.Q(xi1j,xi1j)=model.Q(xi1j,xi1j)+1/2;
%             model.Q(xij,xi1j)=model.Q(xij,xi1j)-2/2;
%             %        model2.Q(xi1j,xij)=-0;
%         end
%     end
% end
% model.Q = -model.Q;

% minimize switches considering the amplitude; positive half amplitude for upswitches
model.Q = sparse(R*n*numberOfUsers+n*numberOfUsers,R*n*numberOfUsers+n*numberOfUsers);
for u = 1:numberOfUsers
    for i=1:n-1
        for j=2:R
                xij = (u-1)*n*R +(i-1)*R+j;
            for k=1:j
                xi1k = (u-1)*n*R +i*R+k;
                model.Q(xij,xi1k)=0.5*(j-k);
                model.Q(xi1k,xij)=0.5*(j-k);
            end
        end
        for j=1:n-1
                xij = (u-1)*n*R +(i-1)*R+j;
            for k=j:R
                xi1k = (u-1)*n*R +i*R+k;
                model.Q(xij,xi1k)=0.5*(k-j);
                model.Q(xi1k,xij)=0.5*(k-j);
            end
        end
    end
end
model.Q = -model.Q;

% maximize Fairness, i.e. difference of quality over users per time slot.
% 1/U * sum_i(sum_u(ABS(sum_j(j*((U-1)*x_u + x_1 + x_2 + ... + x_u-1 + x_u+1 + ... +x_U)))))
% F = sparse(R*n*numberOfUsers,R*n*numberOfUsers);
% for u = 1:numberOfUsers
%     for i=1:n
%         for j=1:R
%                 xij = (u-1)*n*R +(i-1)*R+j;
%             for k=1:j
%                 xi1k = (u-1)*n*R +i*R+k;
%                 model.Q(xij,xi1k)=0.5*(j-k);
%                 model.Q(xi1k,xij)=0.5*(j-k);
%             end
%         end
%     end
% end
% model.Q = -model.Q;

% model.obj = beta .* model.obj;
model.Q = beta .* model.Q./n./numberOfUsers./(R-1);

% minimize qualityvariance of users
% quality = sqrt(sum((model.obj - mean(model.obj)).^2)./n);

model.modelsense = 'max';

%TODO: Implement QoS
% session QoE Fairness
% Idea: just change the order of iteration

%% solve it
result = gurobi(model, params)
%% show it
if ~isfield(result,'x')
    warning('stalling occurs!');
    stat=0;
    return;
else
    %     fprintf('Obj: %d\n', result.objval);
    playedSegment=zeros(n,numberOfUsers);
    downloadedSegments = zeros(n,numberOfUsers);
    for u=1:numberOfUsers
        ind=(u-1)*n*R+(1:n*R);
        solution=reshape(result.x(ind),R,n)';
        [tmp,j]=ind2sub([R,n],find(result.x(ind)));
        downloadedSegments(:,u)=sum(solution,2); % should be
        playedSegment(j,u)=tmp;
    end
end
%%
downloadedSegments1 = downloadedSegments;
playedSegment1=playedSegment;
x=sum(playedSegment);
% printMe('mean quality',x/n)
cx=std(x)/mean(x);
% fprintf('Fairness Index Quality: \t%f\n',1/(1+cx.^2));

stat.meanQualityPerUserStep1=x/n;
stat.meanQualityStep1=mean(x/n);
% stat.fairnessQualityStep1=1/(1+cx.^2);
stat.fairnessQualityStep1 = 1 - 2*std(stat.meanQualityPerUserStep1)/3;

switches=sum(diff(playedSegment)~=0);
% printMe('#switches',switches,'\t%d');
cxswitches=std(switches)/mean(switches);
% fprintf('Fairness Index Switches: \t%f\n',1/(1+cxswitches.^2));

stat.meanSwitchesPerUserStep1=switches;
stat.meanSwitchesStep1=mean(switches);
stat.fairnessSwitchesStep1=1/(1+cxswitches.^2);
%% minimize switches: first contraint
%% minimize switches: second constraint
%% minimize switches: third constraint
%%
% fprintf('=== Optimum: Run %04d Step2 ==================================================\n',runID);
x=sum(playedSegment);
% printMe('mean quality',x/n)
cx=std(x)/mean(x);
% fprintf('Fairness Index Quality: \t%f\n',1/(1+cx.^2));

switches=sum(diff(playedSegment)~=0);
% printMe('#switches',switches,'\t%d');
cxswitches=std(switches)/mean(switches);
% fprintf('Fairness Index Switches: \t%f\n',1/(1+cxswitches.^2));

cx=std(x)/mean(x);
% fprintf('Fairness Index Quality: \t%f\n',1/(1+cx.^2));

stat.meanQualityPerUserStep2=x/n;
stat.meanQualityStep2=mean(x/n);
stat.fairnessQualityStep2=1 - 2*std(stat.meanQualityPerUserStep1)/3;;

x=sum(playedSegment.^2);
cx=std(x)/mean(x);
stat.meanQualityPerUserStep2_with149values=x/n;
stat.meanQualityStep2_with149values=mean(x/n);
stat.fairnessQualityStep2_with149values=1/(1+cx.^2);

stat.meanSwitchesPerUserStep2=switches;
stat.meanSwitchesStep2=mean(switches);
stat.fairnessSwitchesStep2=1/(1+cxswitches.^2);

stat.segmentsPerLayer2=hist(playedSegment,1:R);
stat.playedSegment=playedSegment;


stat.bundle=bundle;
% stat.bw_faktor=bw_faktor;
stat.numberOfUsers=numberOfUsers;
stat.runID=runID;
stat.objective=objective;
% stat.Te=Te;
stat.result = result;
stat.model = model;

% save(int2str(bw_faktor*100));
% plot(stat.playedSegment)