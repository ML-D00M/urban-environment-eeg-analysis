% some ISC processing parameters
alpha = 0.1; % shrinkage parameter; smaller alpha for less regularization
Ncomp = 3;  % number of components to dispaly (all D are computed)
fs=500; % in Hz

%% sea
% Load the data
load('sea_pr_delta.mat','F');
X=F;
clear F

% T samples, D channels, N subjects
[T,D,N] = size(X); 

% now start the ISC code properly

% compute cross-covariance between all subjects i and j
Rij = permute(reshape(cov(X(:,:)),[D N  D N]),[1 3 2 4]); 

% compute within- and between-subject covariances
Rw =       1/N* sum(Rij(:,:,1:N+1:N*N),3);  % pooled over all subjects
Rb = 1/(N-1)/N*(sum(Rij(:,:,:),3) - N*Rw);  % pooled over all pairs of subjects

% shrinkage regularization of Rw
Rw_reg = (1-alpha)*Rw + alpha*mean(eig(Rw))*eye(size(Rw));

% compute correlated components W using regularized Rw, sort components by ISC
[W,ISC]=eig(Rb,Rw_reg); [ISC,indx]=sort(diag(ISC),'descend'); W=W(:,indx);
% ISC

% compute forward model ("scalp projections") A
A_s=Rw*W/(W'*Rw*W);

% ensure that we have only the positive values in A
for j=1:60
    for i=1:60
        if A_s(j,i) < 0
            A_s(j,i) = - A_s(j,i)
        end
    end
end
%% h_msk
% Load the data
load('h_msk_pr_delta.mat','F');
X=F;
clear F

% T samples, D channels, N subjects
[T,D,N] = size(X); 

% now start the ISC code properly

% compute cross-covariance between all subjects i and j
Rij = permute(reshape(cov(X(:,:)),[D N  D N]),[1 3 2 4]); 

% compute within- and between-subject covariances
Rw =       1/N* sum(Rij(:,:,1:N+1:N*N),3);  % pooled over all subjects
Rb = 1/(N-1)/N*(sum(Rij(:,:,:),3) - N*Rw);  % pooled over all pairs of subjects

% shrinkage regularization of Rw
Rw_reg = (1-alpha)*Rw + alpha*mean(eig(Rw))*eye(size(Rw));

% compute correlated components W using regularized Rw, sort components by ISC
[W,ISC]=eig(Rb,Rw_reg); [ISC,indx]=sort(diag(ISC),'descend'); W=W(:,indx);
% ISC

% compute forward model ("scalp projections") A
 A_hm=Rw*W/(W'*Rw*W);

% ensure that we have only the positive values in A
for j=1:60
    for i=1:60
        if A_hm(j,i) < 0
            A_hm(j,i) = - A_hm(j,i)
        end
    end
end
%% h_spb
% Load the data
load('h_spb_pr_delta.mat','F');
X=F;
clear F

% T samples, D channels, N subjects
[T,D,N] = size(X); 

% now start the ISC code properly

% compute cross-covariance between all subjects i and j
Rij = permute(reshape(cov(X(:,:)),[D N  D N]),[1 3 2 4]); 

% compute within- and between-subject covariances
Rw =       1/N* sum(Rij(:,:,1:N+1:N*N),3);  % pooled over all subjects
Rb = 1/(N-1)/N*(sum(Rij(:,:,:),3) - N*Rw);  % pooled over all pairs of subjects

% shrinkage regularization of Rw
Rw_reg = (1-alpha)*Rw + alpha*mean(eig(Rw))*eye(size(Rw));

% compute correlated components W using regularized Rw, sort components by ISC
[W,ISC]=eig(Rb,Rw_reg); [ISC,indx]=sort(diag(ISC),'descend'); W=W(:,indx);
% ISC

% compute forward model ("scalp projections") A
 A_hs=Rw*W/(W'*Rw*W);

% ensure that we have only the positive values in A
for j=1:60
    for i=1:60
        if A_hs(j,i) < 0
            A_hs(j,i) = - A_hs(j,i)
        end
    end
end
%% b_msk
% Load the data
load('b_msk_pr_delta.mat','F');
X=F;
clear F

% T samples, D channels, N subjects
[T,D,N] = size(X); 

% now start the ISC code properly

% compute cross-covariance between all subjects i and j
Rij = permute(reshape(cov(X(:,:)),[D N  D N]),[1 3 2 4]); 

% compute within- and between-subject covariances
Rw =       1/N* sum(Rij(:,:,1:N+1:N*N),3);  % pooled over all subjects
Rb = 1/(N-1)/N*(sum(Rij(:,:,:),3) - N*Rw);  % pooled over all pairs of subjects

% shrinkage regularization of Rw
Rw_reg = (1-alpha)*Rw + alpha*mean(eig(Rw))*eye(size(Rw));

% compute correlated components W using regularized Rw, sort components by ISC
[W,ISC]=eig(Rb,Rw_reg); [ISC,indx]=sort(diag(ISC),'descend'); W=W(:,indx);
% ISC

% compute forward model ("scalp projections") A
 A_bm=Rw*W/(W'*Rw*W);

% ensure that we have only the positive values in A
for j=1:60
    for i=1:60
        if A_bm(j,i) < 0
            A_bm(j,i) = - A_bm(j,i)
        end
    end
end
%% b_spb
% Load the data
load('b_spb_pr_delta.mat','F');
X=F;
clear F

% T samples, D channels, N subjects
[T,D,N] = size(X); 

% now start the ISC code properly

% compute cross-covariance between all subjects i and j
Rij = permute(reshape(cov(X(:,:)),[D N  D N]),[1 3 2 4]); 

% compute within- and between-subject covariances
Rw =       1/N* sum(Rij(:,:,1:N+1:N*N),3);  % pooled over all subjects
Rb = 1/(N-1)/N*(sum(Rij(:,:,:),3) - N*Rw);  % pooled over all pairs of subjects

% shrinkage regularization of Rw
Rw_reg = (1-alpha)*Rw + alpha*mean(eig(Rw))*eye(size(Rw));

% compute correlated components W using regularized Rw, sort components by ISC
[W,ISC]=eig(Rb,Rw_reg); [ISC,indx]=sort(diag(ISC),'descend'); W=W(:,indx);
% ISC

% compute forward model ("scalp projections") A
 A_bs=Rw*W/(W'*Rw*W);

% ensure that we have only the positive values in A
for j=1:60
    for i=1:60
        if A_bs(j,i) < 0
            A_bs(j,i) = - A_bs(j,i)
        end
    end
end
%% p_msk
% Load the data
load('p_msk_pr_delta.mat','F');
X=F;
clear F

% T samples, D channels, N subjects
[T,D,N] = size(X); 

% now start the ISC code properly

% compute cross-covariance between all subjects i and j
Rij = permute(reshape(cov(X(:,:)),[D N  D N]),[1 3 2 4]); 

% compute within- and between-subject covariances
Rw =       1/N* sum(Rij(:,:,1:N+1:N*N),3);  % pooled over all subjects
Rb = 1/(N-1)/N*(sum(Rij(:,:,:),3) - N*Rw);  % pooled over all pairs of subjects

% shrinkage regularization of Rw
Rw_reg = (1-alpha)*Rw + alpha*mean(eig(Rw))*eye(size(Rw));

% compute correlated components W using regularized Rw, sort components by ISC
[W,ISC]=eig(Rb,Rw_reg); [ISC,indx]=sort(diag(ISC),'descend'); W=W(:,indx);
% ISC

% compute forward model ("scalp projections") A
 A_pm=Rw*W/(W'*Rw*W);

% ensure that we have only the positive values in A
for j=1:60
    for i=1:60
        if A_pm(j,i) < 0
            A_pm(j,i) = - A_pm(j,i)
        end
    end
end
%% p_spb
% Load the data
load('p_spb_pr_delta.mat','F');
X=F;
clear F

% T samples, D channels, N subjects
[T,D,N] = size(X); 

% now start the ISC code properly

% compute cross-covariance between all subjects i and j
Rij = permute(reshape(cov(X(:,:)),[D N  D N]),[1 3 2 4]); 

% compute within- and between-subject covariances
Rw =       1/N* sum(Rij(:,:,1:N+1:N*N),3);  % pooled over all subjects
Rb = 1/(N-1)/N*(sum(Rij(:,:,:),3) - N*Rw);  % pooled over all pairs of subjects

% shrinkage regularization of Rw
Rw_reg = (1-alpha)*Rw + alpha*mean(eig(Rw))*eye(size(Rw));

% compute correlated components W using regularized Rw, sort components by ISC
[W,ISC]=eig(Rb,Rw_reg); [ISC,indx]=sort(diag(ISC),'descend'); W=W(:,indx);
% ISC

% compute forward model ("scalp projections") A
 A_ps=Rw*W/(W'*Rw*W);

% ensure that we have only the positive values in A
for j=1:60
    for i=1:60
        if A_ps(j,i) < 0
            A_ps(j,i) = - A_ps(j,i)
        end
    end
end
%% mov
% Load the data
load('mov_pr_delta.mat','F');
X=F;
clear F

% T samples, D channels, N subjects
[T,D,N] = size(X); 

% now start the ISC code properly

% compute cross-covariance between all subjects i and j
Rij = permute(reshape(cov(X(:,:)),[D N  D N]),[1 3 2 4]); 

% compute within- and between-subject covariances
Rw =       1/N* sum(Rij(:,:,1:N+1:N*N),3);  % pooled over all subjects
Rb = 1/(N-1)/N*(sum(Rij(:,:,:),3) - N*Rw);  % pooled over all pairs of subjects

% shrinkage regularization of Rw
Rw_reg = (1-alpha)*Rw + alpha*mean(eig(Rw))*eye(size(Rw));

% compute correlated components W using regularized Rw, sort components by ISC
[W,ISC]=eig(Rb,Rw_reg); [ISC,indx]=sort(diag(ISC),'descend'); W=W(:,indx);
% ISC

% compute forward model ("scalp projections") A
 A_m=Rw*W/(W'*Rw*W);

% ensure that we have only the positive values in A
for j=1:60
    for i=1:60
        if A_m(j,i) < 0
            A_m(j,i) = - A_m(j,i)
        end
    end
end
%% show some results

% sea and mov
for i=1:5
subplot(2,5,i);
    topoplot(A_s(:,i),'BioSemi64_edit.loc','electrodes','on');
    title(['Component ' num2str(i)])
    subplot(2,5,i+5);
    topoplot(A_m(:,i),'BioSemi64_edit.loc','electrodes','on');
    title(['Component ' num2str(i)])
end
colorbar

% 6 urban
for i=1:5
    subplot(6,5,i);
    topoplot(A_hm(:,i),'BioSemi64_edit.loc','electrodes','on');
    title(['Component ' num2str(i)])
    subplot(6,5,i+5);
    topoplot(A_hs(:,i),'BioSemi64_edit.loc','electrodes','on');
    title(['Component ' num2str(i)])
    subplot(6,5,i+10);
    topoplot(A_bm(:,i),'BioSemi64_edit.loc','electrodes','on');
    title(['Component ' num2str(i)])
    subplot(6,5,i+15);
    topoplot(A_bs(:,i),'BioSemi64_edit.loc','electrodes','on');
    title(['Component ' num2str(i)])
    subplot(6,5,i+20);
    topoplot(A_pm(:,i),'BioSemi64_edit.loc','electrodes','on');
    title(['Component ' num2str(i)])
    subplot(6,5,i+25);
    topoplot(A_ps(:,i),'BioSemi64_edit.loc','electrodes','on');
    title(['Component ' num2str(i)])
end
colorbar

% parks
for i=1:5
    subplot(2,5,i);
    topoplot(A_pm(:,i),'BioSemi64_edit.loc','electrodes','on');
    title(['Component ' num2str(i)])
    subplot(2,5,i+5);
    topoplot(A_ps(:,i),'BioSemi64_edit.loc','electrodes','on');
    title(['Component ' num2str(i)])
end
colorbar

% highways
for i=1:5
    subplot(2,5,i);
    topoplot(A_hm(:,i),'BioSemi64_edit.loc','electrodes','on');
    title(['Component ' num2str(i)])
    subplot(2,5,i+5);
    topoplot(A_hs(:,i),'BioSemi64_edit.loc','electrodes','on');
    title(['Component ' num2str(i)])
end
colorbar

% boulevards
for i=1:5
    subplot(2,5,i);
    topoplot(A_bm(:,i),'BioSemi64_edit.loc','electrodes','on');
    title(['Component ' num2str(i)])
    subplot(2,5,i+5);
    topoplot(A_bs(:,i),'BioSemi64_edit.loc','electrodes','on');
    title(['Component ' num2str(i)])
end
colorbar

% msk
for i=1:5
    subplot(3,5,i);
    topoplot(A_hm(:,i),'BioSemi64_edit.loc','electrodes','on');
    title(['Component ' num2str(i)])
    subplot(3,5,i+5);
    topoplot(A_bm(:,i),'BioSemi64_edit.loc','electrodes','on');
    title(['Component ' num2str(i)])
    subplot(3,5,i+10);
    topoplot(A_pm(:,i),'BioSemi64_edit.loc','electrodes','on');
    title(['Component ' num2str(i)])
end
colorbar

% spb
for i=1:5
    subplot(3,5,i);
    topoplot(A_hs(:,i),'BioSemi64_edit.loc','electrodes','on');
    title(['Component ' num2str(i)])
    subplot(3,5,i+5);
    topoplot(A_bs(:,i),'BioSemi64_edit.loc','electrodes','on');
    title(['Component ' num2str(i)])
    subplot(3,5,i+10);
    topoplot(A_ps(:,i),'BioSemi64_edit.loc','electrodes','on');
    title(['Component ' num2str(i)])
end
colorbar

%% The whole signal topoplots comparing to delta

%% sea
% Load the data
load('sea_pr.mat','PR');
X=PR;
clear PR

% T samples, D channels, N subjects
[T,D,N] = size(X); 

% now start the ISC code properly

% compute cross-covariance between all subjects i and j
Rij = permute(reshape(cov(X(:,:)),[D N  D N]),[1 3 2 4]); 

% compute within- and between-subject covariances
Rw =       1/N* sum(Rij(:,:,1:N+1:N*N),3);  % pooled over all subjects
Rb = 1/(N-1)/N*(sum(Rij(:,:,:),3) - N*Rw);  % pooled over all pairs of subjects

% shrinkage regularization of Rw
Rw_reg = (1-alpha)*Rw + alpha*mean(eig(Rw))*eye(size(Rw));

% compute correlated components W using regularized Rw, sort components by ISC
[W,ISC]=eig(Rb,Rw_reg); [ISC,indx]=sort(diag(ISC),'descend'); W=W(:,indx);
% ISC

% compute forward model ("scalp projections") A
 AA_s=Rw*W/(W'*Rw*W);

% ensure that we have only the positive values in A
for j=1:60
    for i=1:60
        if AA_s(j,i) < 0
            AA_s(j,i) = - AA_s(j,i)
        end
    end
end
%% h_msk
% Load the data
load('h_msk_pr.mat','PR');
X=PR;
clear PR

% T samples, D channels, N subjects
[T,D,N] = size(X); 

% now start the ISC code properly

% compute cross-covariance between all subjects i and j
Rij = permute(reshape(cov(X(:,:)),[D N  D N]),[1 3 2 4]); 

% compute within- and between-subject covariances
Rw =       1/N* sum(Rij(:,:,1:N+1:N*N),3);  % pooled over all subjects
Rb = 1/(N-1)/N*(sum(Rij(:,:,:),3) - N*Rw);  % pooled over all pairs of subjects

% shrinkage regularization of Rw
Rw_reg = (1-alpha)*Rw + alpha*mean(eig(Rw))*eye(size(Rw));

% compute correlated components W using regularized Rw, sort components by ISC
[W,ISC]=eig(Rb,Rw_reg); [ISC,indx]=sort(diag(ISC),'descend'); W=W(:,indx);
% ISC

% compute forward model ("scalp projections") A
 AA_hm=Rw*W/(W'*Rw*W);

% ensure that we have only the positive values in A
for j=1:60
    for i=1:60
        if AA_hm(j,i) < 0
            AA_hm(j,i) = - AA_hm(j,i)
        end
    end
end
%% h_spb
% Load the data
load('h_spb_pr.mat','PR');
X=PR;
clear PR

% T samples, D channels, N subjects
[T,D,N] = size(X); 

% now start the ISC code properly

% compute cross-covariance between all subjects i and j
Rij = permute(reshape(cov(X(:,:)),[D N  D N]),[1 3 2 4]); 

% compute within- and between-subject covariances
Rw =       1/N* sum(Rij(:,:,1:N+1:N*N),3);  % pooled over all subjects
Rb = 1/(N-1)/N*(sum(Rij(:,:,:),3) - N*Rw);  % pooled over all pairs of subjects

% shrinkage regularization of Rw
Rw_reg = (1-alpha)*Rw + alpha*mean(eig(Rw))*eye(size(Rw));

% compute correlated components W using regularized Rw, sort components by ISC
[W,ISC]=eig(Rb,Rw_reg); [ISC,indx]=sort(diag(ISC),'descend'); W=W(:,indx);
% ISC

% compute forward model ("scalp projections") A
 AA_hs=Rw*W/(W'*Rw*W);

% ensure that we have only the positive values in A
for j=1:60
    for i=1:60
        if AA_hs(j,i) < 0
            AA_hs(j,i) = - AA_hs(j,i)
        end
    end
end
%% b_msk
% Load the data
load('b_msk_pr.mat','PR');
X=PR;
clear PR

% T samples, D channels, N subjects
[T,D,N] = size(X); 

% now start the ISC code properly

% compute cross-covariance between all subjects i and j
Rij = permute(reshape(cov(X(:,:)),[D N  D N]),[1 3 2 4]); 

% compute within- and between-subject covariances
Rw =       1/N* sum(Rij(:,:,1:N+1:N*N),3);  % pooled over all subjects
Rb = 1/(N-1)/N*(sum(Rij(:,:,:),3) - N*Rw);  % pooled over all pairs of subjects

% shrinkage regularization of Rw
Rw_reg = (1-alpha)*Rw + alpha*mean(eig(Rw))*eye(size(Rw));

% compute correlated components W using regularized Rw, sort components by ISC
[W,ISC]=eig(Rb,Rw_reg); [ISC,indx]=sort(diag(ISC),'descend'); W=W(:,indx);
% ISC

% compute forward model ("scalp projections") A
 AA_bm=Rw*W/(W'*Rw*W);

% ensure that we have only the positive values in A
for j=1:60
    for i=1:60
        if AA_bm(j,i) < 0
            AA_bm(j,i) = - AA_bm(j,i)
        end
    end
end
%% b_spb
% Load the data
load('b_spb_pr.mat','PR');
X=PR;
clear PR

% T samples, D channels, N subjects
[T,D,N] = size(X); 

% now start the ISC code properly

% compute cross-covariance between all subjects i and j
Rij = permute(reshape(cov(X(:,:)),[D N  D N]),[1 3 2 4]); 

% compute within- and between-subject covariances
Rw =       1/N* sum(Rij(:,:,1:N+1:N*N),3);  % pooled over all subjects
Rb = 1/(N-1)/N*(sum(Rij(:,:,:),3) - N*Rw);  % pooled over all pairs of subjects

% shrinkage regularization of Rw
Rw_reg = (1-alpha)*Rw + alpha*mean(eig(Rw))*eye(size(Rw));

% compute correlated components W using regularized Rw, sort components by ISC
[W,ISC]=eig(Rb,Rw_reg); [ISC,indx]=sort(diag(ISC),'descend'); W=W(:,indx);
% ISC

% compute forward model ("scalp projections") A
 AA_bs=Rw*W/(W'*Rw*W);
 
% ensure that we have only the positive values in A
for j=1:60
    for i=1:60
        if AA_bs(j,i) < 0
            AA_bs(j,i) = - AA_bs(j,i)
        end
    end
end

%% p_msk
% Load the data
load('p_msk_pr.mat','PR');
X=PR;
clear PR

% T samples, D channels, N subjects
[T,D,N] = size(X); 

% now start the ISC code properly

% compute cross-covariance between all subjects i and j
Rij = permute(reshape(cov(X(:,:)),[D N  D N]),[1 3 2 4]); 

% compute within- and between-subject covariances
Rw =       1/N* sum(Rij(:,:,1:N+1:N*N),3);  % pooled over all subjects
Rb = 1/(N-1)/N*(sum(Rij(:,:,:),3) - N*Rw);  % pooled over all pairs of subjects

% shrinkage regularization of Rw
Rw_reg = (1-alpha)*Rw + alpha*mean(eig(Rw))*eye(size(Rw));

% compute correlated components W using regularized Rw, sort components by ISC
[W,ISC]=eig(Rb,Rw_reg); [ISC,indx]=sort(diag(ISC),'descend'); W=W(:,indx);
% ISC

% compute forward model ("scalp projections") A
 AA_pm=Rw*W/(W'*Rw*W);

% ensure that we have only the positive values in A
for j=1:60
    for i=1:60
        if AA_pm(j,i) < 0
            AA_pm(j,i) = - AA_pm(j,i)
        end
    end
end
%% p_spb
% Load the data
load('p_spb_pr.mat','PR');
X=PR;
clear PR

% T samples, D channels, N subjects
[T,D,N] = size(X); 

% now start the ISC code properly

% compute cross-covariance between all subjects i and j
Rij = permute(reshape(cov(X(:,:)),[D N  D N]),[1 3 2 4]); 

% compute within- and between-subject covariances
Rw =       1/N* sum(Rij(:,:,1:N+1:N*N),3);  % pooled over all subjects
Rb = 1/(N-1)/N*(sum(Rij(:,:,:),3) - N*Rw);  % pooled over all pairs of subjects

% shrinkage regularization of Rw
Rw_reg = (1-alpha)*Rw + alpha*mean(eig(Rw))*eye(size(Rw));

% compute correlated components W using regularized Rw, sort components by ISC
[W,ISC]=eig(Rb,Rw_reg); [ISC,indx]=sort(diag(ISC),'descend'); W=W(:,indx);
% ISC

% compute forward model ("scalp projections") A
 AA_ps=Rw*W/(W'*Rw*W);

% ensure that we have only the positive values in A
for j=1:60
    for i=1:60
        if AA_ps(j,i) < 0
            AA_ps(j,i) = - AA_ps(j,i)
        end
    end
end
%% mov
% Load the data
load('mov_pr.mat','PR');
X=PR;
clear PR

% T samples, D channels, N subjects
[T,D,N] = size(X); 

% now start the ISC code properly

% compute cross-covariance between all subjects i and j
Rij = permute(reshape(cov(X(:,:)),[D N  D N]),[1 3 2 4]); 

% compute within- and between-subject covariances
Rw =       1/N* sum(Rij(:,:,1:N+1:N*N),3);  % pooled over all subjects
Rb = 1/(N-1)/N*(sum(Rij(:,:,:),3) - N*Rw);  % pooled over all pairs of subjects

% shrinkage regularization of Rw
Rw_reg = (1-alpha)*Rw + alpha*mean(eig(Rw))*eye(size(Rw));

% compute correlated components W using regularized Rw, sort components by ISC
[W,ISC]=eig(Rb,Rw_reg); [ISC,indx]=sort(diag(ISC),'descend'); W=W(:,indx);
% ISC

% compute forward model ("scalp projections") A
 AA_m=Rw*W/(W'*Rw*W);

% ensure that we have only the positive values in A
for j=1:60
    for i=1:60
        if AA_m(j,i) < 0
            AA_m(j,i) = - AA_m(j,i)
        end
    end
end
%% show some results

% sea
for i=1:5
subplot(2,5,i);
    topoplot(A_s(:,i),'BioSemi64_edit.loc','electrodes','on');
    title(['Component ' num2str(i)])
    subplot(2,5,i+5);
    topoplot(AA_s(:,i),'BioSemi64_edit.loc','electrodes','on');
    title(['Component ' num2str(i)])
end
colorbar

% mov
for i=1:5
subplot(2,5,i);
    topoplot(A_m(:,i),'BioSemi64_edit.loc','electrodes','on');
    title(['Component ' num2str(i)])
    subplot(2,5,i+5);
    topoplot(AA_m(:,i),'BioSemi64_edit.loc','electrodes','on');
    title(['Component ' num2str(i)])
end
colorbar

% h_msk
for i=1:5
subplot(2,5,i);
    topoplot(A_hm(:,i),'BioSemi64_edit.loc','electrodes','on');
    title(['Component ' num2str(i)])
    subplot(2,5,i+5);
    topoplot(AA_hm(:,i),'BioSemi64_edit.loc','electrodes','on');
    title(['Component ' num2str(i)])
end
colorbar

% h_spb
for i=1:5
subplot(2,5,i);
    topoplot(A_hs(:,i),'BioSemi64_edit.loc','electrodes','on');
    title(['Component ' num2str(i)])
    subplot(2,5,i+5);
    topoplot(AA_hs(:,i),'BioSemi64_edit.loc','electrodes','on');
    title(['Component ' num2str(i)])
end
colorbar

% b_msk
for i=1:5
subplot(2,5,i);
    topoplot(A_bm(:,i),'BioSemi64_edit.loc','electrodes','on');
    title(['Component ' num2str(i)])
    subplot(2,5,i+5);
    topoplot(AA_bm(:,i),'BioSemi64_edit.loc','electrodes','on');
    title(['Component ' num2str(i)])
end
colorbar

% b_spb
for i=1:5
subplot(2,5,i);
    topoplot(A_bs(:,i),'BioSemi64_edit.loc','electrodes','on');
    title(['Component ' num2str(i)])
    subplot(2,5,i+5);
    topoplot(AA_bs(:,i),'BioSemi64_edit.loc','electrodes','on');
    title(['Component ' num2str(i)])
end
colorbar

% p_msk
for i=1:5
subplot(2,5,i);
    topoplot(A_pm(:,i),'BioSemi64_edit.loc','electrodes','on');
    title(['Component ' num2str(i)])
    subplot(2,5,i+5);
    topoplot(AA_pm(:,i),'BioSemi64_edit.loc','electrodes','on');
    title(['Component ' num2str(i)])
end
colorbar

% p_spb
for i=1:5
subplot(2,5,i);
    topoplot(A_ps(:,i),'BioSemi64_edit.loc','electrodes','on');
    title(['Component ' num2str(i)])
    subplot(2,5,i+5);
    topoplot(AA_ps(:,i),'BioSemi64_edit.loc','electrodes','on');
    title(['Component ' num2str(i)])
end
colorbar