% Step 5. Inter-Subject Correlation computations

%% Setting parameters and loading data

% some ISC processing parameters
gamma = 0.1; % shrinkage parameter; smaller gamma for less regularization
Ncomp = 3;  % number of components to dispaly (all D are computed)
fs=500; % in Hz

% Load the data
load('h_spb_pr.mat','PR');
X=PR;
% X=X(:,setdiff(1:D,[5 25]),:)
clear PR


%% Computing cross-covariance between subjects, correlated components, ISCs and forward model

% Take only subjects with good signal.
% e.g. suppose that only subjects 4, 5, 7 and 12 have good signal.
% X=X(:,:,[4 5 7 12]);

% T samples, D channels, N subjects
[T,D,N] = size(X); 

% now start the ISC code proper

% compute cross-covariance between all subjects i and j
Rij = permute(reshape(cov(X(:,:)),[D N  D N]),[1 3 2 4]); 

% compute within- and between-subject covariances
Rw =       1/N* sum(Rij(:,:,1:N+1:N*N),3);  % pooled over all subjects
Rb = 1/(N-1)/N*(sum(Rij(:,:,:),3) - N*Rw);  % pooled over all pairs of subjects

% shrinkage regularization of Rw
Rw_reg = (1-gamma)*Rw + gamma*mean(eig(Rw))*eye(size(Rw));

% compute correlated components W using regularized Rw, sort components by ISC
[W,ISC]=eig(Rb,Rw_reg); [ISC,indx]=sort(diag(ISC),'descend'); W=W(:,indx);
% ISC

% compute forward model ("scalp projections") A
 A=Rw*W/(W'*Rw*W);


%% Computing ISC resolved by subject

% +++ If multiple stimuli are available, then Rij as computed for each stimulus
% should be used in the following to compute ISC_persubject, and
% ISC_persecond +++

% Compute ISC resolved by subject, see Cohen et al.
% Every row corresponds to one component, every column
% corresponds to each subject.
for i=1:N
    Rw=0; for j=1:N, if i~=j, Rw = Rw+1/(N-1)*(Rij(:,:,i,i)+Rij(:,:,j,j)); end; end
    Rb=0; for j=1:N, if i~=j, Rb = Rb+1/(N-1)*(Rij(:,:,i,j)+Rij(:,:,j,i)); end; end
    ISC_persubject(:,i) = diag(W'*Rb*W)./diag(W'*Rw*W);
end

% Sum the 3 strongest components and save it as a vector
ISC_persubject_sum=sum(ISC_persubject(1:3,:),1);
save(['h_spb_ISC_persubject.mat'],'ISC_persubject_sum');


%% Computing ISC across time windows

Nsec  = 1;  % time-window (in seconds) over which to compute time-reposeved ISC
overlap = 0.8*Nsec;
t = 1:floor((Nsec-overlap)*fs+1):(T-fs*Nsec);
for tau=1:length(t)
    Xt = X(t(tau):(t(tau)+Nsec*fs-1),:,:);
    Rij = permute(reshape(cov(Xt(:,:)),[D N  D N]),[1 3 2 4]);
    Rw =  1/N* sum(Rij(:,:,1:N+1:N*N),3);  % pooled over all subjects
    Rb = 1/(N-1)/N*(sum(Rij(:,:,:),3) - N*Rw);  % pooled over all pairs of subjects
    ISC_persecond(:,tau) = diag(W'*Rb*W)./diag(W'*Rw*W);
end


%% Plotting ISC time series

% Plot the ISC time series (sum of the first 3 components)
avg_vis=0.5;
for tau=2:length(t)
    avg_vis(tau)=avg_vis(tau-1)+0.2;
end
ISC_persecond_sum=sum(ISC_persecond(1:3,:),1);

figure(1)
subplot(3,1,1)
plot(avg_vis,ISC_persecond_sum(:,:))
ylim([-0.05 0.25]);
xlabel('Seconds');
title('SPb Highway video');
ylabel('ISC');


% Plot the ISC time series of the first 3 components
avg_vis=0.5;
for tau=2:length(t)
    avg_vis(tau)=avg_vis(tau-1)+0.2;
end

figure(2)
for i=1:3
    subplot(3,1,i)
    plot(avg_vis,ISC_persecond(i,:))
    ylim([-0.05 0.2]);
    xlabel('Seconds')
    title(['Component ' num2str(i)])
    ylabel('ISC')
end


%% Checking for statistical significance

num_permutations = 100;
permutation_step = 10;

fprintf('Running %d permutations...', num_permutations)
chance_val = [];
for iter = 1:num_permutations
    Xr = phase_randomized(X);
    
    if mod(iter, permutation_step) == 0
        fprintf('\n%d permutations done...', iter)
    end

    t = 1:floor((Nsec-overlap)*fs+1):(T-fs*Nsec);
    for tau=1:length(t)
        Xt = Xr(t(tau):(t(tau)+Nsec*fs-1),:,:);
        Rij = permute(reshape(cov(Xt(:,:)),[D N  D N]),[1 3 2 4]);
        Rw = 1/N* sum(Rij(:,:,1:N+1:N*N),3);  % pooled over all subjects
        Rb = 1/(N-1)/N*(sum(Rij(:,:,:),3) - N*Rw);  % pooled over all pairs of subjects
        rISC_perwindow(:,tau) = diag(W'*Rb*W)./diag(W'*Rw*W);
    end
    chance_val(:,:,iter) = rISC_perwindow';
end

%% Getting p-values

% get the average ISC value of the null distribution for further plotting

%avg_null_isc = mean(mean(chance_val, 3));
%avg_null_isc2 = mean(chance_val, 3);

fprintf('Getting the p-values...')
pvals=stat_surrogate_pvals(chance_val,ISC_persecond');

% Adjust for multiple comparisons

save('pvals_h_spb.mat','pvals');
%p1=pval_adjust(pvals(:,1), 'BH');
%p2=pval_adjust(pvals(:,2), 'BH');
%p3=pval_adjust(pvals(:,3), 'BH');

p1=pvals(:,1);
p2=pvals(:,2);
p3=pvals(:,3);

%Check what is the percentage of significant windows in the 1st component
length(find(p1<0.05))/length(p1) % 0.56 => ISC in 56% of windows is significant
length(find(p2<0.05))/length(p2)
length(find(p3<0.05))/length(p3)

% Since we do not adjust pvalues, we can use more strict threshold
length(find(p1<0.01))/length(p1)
length(find(p2<0.01))/length(p2)
length(find(p3<0.01))/length(p3)

%% Plotting ISC time series with the threshold — null distribution average

% Plot the ISC time series of the first 3 components + avg ISC of the null
% distribution

%sorted_null_isc = sort(mean(mean(chance_val, 3)), 'descend');
avg_null_isc2 = mean(chance_val, 3)'; % we calculate the avg ISCs of the null
% distribution for each time-window

avg_vis=0.5;
for tau=2:length(t)
    avg_vis(tau)=avg_vis(tau-1)+0.2;
end
figure(3)
for i=1:3
    subplot(3,1,i)
    plot(avg_vis,ISC_persecond(i,:))
    hold on
    % Plot the average ISC of the null distribution with a thinner grey line
    plot(avg_vis, avg_null_isc2(i,:), 'color', [0.5 0.5 0.5], 'LineWidth', 0.25)

    % Fill the area under the average ISC of the null distribution with a
    % translucent grey color
    fill([avg_vis, fliplr(avg_vis)], [avg_null_isc2(i,:), fliplr(-0.05*ones(1,length(avg_vis)))], [0.5 0.5 0.5], 'FaceAlpha', 0.3);
    ylim([-0.05 0.2]);
    xlabel('Seconds')
    title(['Component ' num2str(i)])
    ylabel('ISC')
end

% Plot the ISC time series (sum of the first 3 components) + threshold
avg_vis=0.5;
for tau=2:length(t)
    avg_vis(tau)=avg_vis(tau-1)+0.2;
end
avg_null_isc2_sum=sum(avg_null_isc2(1:3,:),1);

figure(4)
subplot(3,1,1)
plot(avg_vis,ISC_persecond_sum(:,:))
hold on
% Plot the average ISC of the null distribution with a thinner grey line
plot(avg_vis, avg_null_isc2_sum(:,:), 'color', [0.5 0.5 0.5], 'LineWidth', 0.25)

% Fill the area under the average ISC of the null distribution with a
% translucent grey color
fill([avg_vis, fliplr(avg_vis)], [avg_null_isc2_sum(:,:), fliplr(-0.05*ones(1,length(avg_vis)))], [0.5 0.5 0.5], 'FaceAlpha', 0.3);

ylim([-0.05 0.25]);
xlabel('Seconds');
title('SPb Highway video');
ylabel('ISC');

%% Plotting ISC time series with the threshold — null distribution 0.95 percentile

% Plot the ISC time series of the first 3 components + 0.05 threshold of the null
% distribution

%sorted_null_isc = sort(mean(mean(chance_val, 3)), 'descend');
avg_null_isc2 = mean(chance_val, 3)'; % we calculate the avg ISCs of the null
% distribution for each time-window

% Calculate 0.05 threshold for each time-window
null_threshold_05 = prctile(chance_val, 95, 3)';

% Calculate the max value of the null distribution for each time-window
max_null_isc = max(chance_val, [], 3)';

avg_vis=0.5;
for tau=2:length(t)
    avg_vis(tau)=avg_vis(tau-1)+0.2;
end

figure(5)
for i=1:3
    subplot(3,1,i)
    plot(avg_vis,ISC_persecond(i,:))
    hold on
    
    % Plot the 95th percentile of the null distribution with a thinner grey line
    plot(avg_vis, null_threshold_05(i,:), 'color', [0.5 0.5 0.5], 'LineWidth', 0.25)

    % Fill the area under the average ISC of the null distribution with a
    % translucent grey color
    fill([avg_vis, fliplr(avg_vis)], [null_threshold_05(i,:), fliplr(-0.05*ones(1,length(avg_vis)))], [0.5 0.5 0.5], 'FaceAlpha', 0.3);
    ylim([-0.05 0.2]);
    xlabel('Seconds')
    title(['Component ' num2str(i)])
    ylabel('ISC')
end

% Plot the ISC time series (sum of the first 3 components) + threshold
avg_vis=0.5;
for tau=2:length(t)
    avg_vis(tau)=avg_vis(tau-1)+0.2;
end
null_threshold_05_sum=sum(null_threshold_05(1:3,:),1);

figure(6)
subplot(3,1,1)
plot(avg_vis,ISC_persecond_sum(:,:))
hold on
% Plot the 95th percentile of the null distribution with a thinner grey line
plot(avg_vis, null_threshold_05_sum(:,:), 'color', [0.5 0.5 0.5], 'LineWidth', 0.25)
% Fill the area under the average ISC of the null distribution with a
% translucent grey color
fill([avg_vis, fliplr(avg_vis)], [null_threshold_05_sum(:,:), fliplr(-0.05*ones(1,length(avg_vis)))], [0.5 0.5 0.5], 'FaceAlpha', 0.3);
ylim([-0.05 0.25]);
xlabel('Seconds');
title('Sea video');
ylabel('ISC');

%% Plotting ISC time series with the threshold — null distribution 0.99 percentile

% Plot the ISC time series of the first 3 components + 0.01 threshold of the null
% distribution

% Calculate the upper 0.01 threshold for each time-window
null_threshold_01 = prctile(chance_val, 99, 3)';

% Calculate the max value of the null distribution for each time-window
% max_null_isc = max(chance_val, [], 3)';

avg_vis=0.5;
for tau=2:length(t)
    avg_vis(tau)=avg_vis(tau-1)+0.2;
end

figure(7)
for i=1:3
    subplot(3,1,i)
    plot(avg_vis,ISC_persecond(i,:))
    hold on
    
    % Plot the 95th percentile of the null distribution with a thinner grey line
    plot(avg_vis, null_threshold_01(i,:), 'color', [0.5 0.5 0.5], 'LineWidth', 0.25)

    % Fill the area under the average ISC of the null distribution with a
    % translucent grey color
    fill([avg_vis, fliplr(avg_vis)], [null_threshold_01(i,:), fliplr(-0.05*ones(1,length(avg_vis)))], [0.5 0.5 0.5], 'FaceAlpha', 0.3);
    ylim([-0.05 0.2]);
    xlabel('Seconds')
    title(['Component ' num2str(i)])
    ylabel('ISC')
end

% Plot the ISC time series (sum of the first 3 components) + threshold
avg_vis=0.5;
for tau=2:length(t)
    avg_vis(tau)=avg_vis(tau-1)+0.2;
end
null_threshold_01_sum=sum(null_threshold_01(1:3,:),1);

figure(8)
subplot(3,1,1)
plot(avg_vis,ISC_persecond_sum(:,:))
hold on
% Plot the 99th percentile of the null distribution with a thinner grey line
plot(avg_vis, null_threshold_01_sum(:,:), 'color', [0.5 0.5 0.5], 'LineWidth', 0.25)
% Fill the area under the average ISC of the null distribution with a
% translucent grey color
fill([avg_vis, fliplr(avg_vis)], [null_threshold_01_sum(:,:), fliplr(-0.05*ones(1,length(avg_vis)))], [0.5 0.5 0.5], 'FaceAlpha', 0.3);
ylim([-0.05 0.25]);
xlabel('Seconds');
title('Sea video');
ylabel('ISC');

%% Plotting the forward model on a topoggraphic map

% show some results
global_min = -0.2;
global_max = 0.2;
figure(9)
for i=1:3
    subplot(1,3,i);
    topoplot(A(:,i),'BioSemi64_edit.loc','electrodes','on'); title(['Component ' num2str(i)])
    clim([global_min, global_max]);
end
colorbar

% plot key components on the top of each other
%notBoxPlot(ISC_persubject(1:Ncomp,:)'); xlabel('Component'); ylabel('ISC'); title('Per subjects');
%plot(ISC_persecond(1:Ncomp,:)'); xlabel('Time (s)'); ylabel('ISC'); title('Per second');

