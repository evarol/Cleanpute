    clear all
clc
close all




%% PARAMETERS
%%%%%%%%%%%%%
lambda=1;%% Bulk proportion sparsity
lambda2=1;% imputation sparsity
lambda3=1;%% bulk purity
lambda4=1;%% soft constraint that imputed gene counts should match the confident gene counts
hard_constraint=1; %%whether the "confident" gene counts should stay as is
covariance_correction=1; %%whether the covariance structure of single cell reads should be taken into account
iter=100; %%number of iterations of algorithm
imputation_cutoff=10; %% Any single cell reads below this will be counted as unreliable and will be imputed
%%%%%%%%%%%


%% helper functions
vec=@(x)(x(:));
addpath(genpath('/Users/erdem/Dropbox/matlab_toolboxes/export-fig/')); %% Replace with local copy of export_fig to print figs

%% read data
tissue_single=readtable('./data/Tissue_level_pseudobulk_11620.csv');
neuron_single=readtable('./data/Neuron_level_pseudobulk_11620.csv');
bulk = readtable('./data/sc_matched_bulk_11620.csv');

%% extract matrices from data
tissue_genes=table2array(tissue_single(:,1));
tissue_cells=tissue_single.Properties.VariableNames(2:end);
tissue_counts=table2array(tissue_single(:,2:end));

neuron_genes=table2array(neuron_single(:,1));
neuron_cells=neuron_single.Properties.VariableNames(2:end);
neuron_counts=table2array(neuron_single(:,2:end));

bulk_genes=table2array(bulk(:,1));
bulk_cells=bulk.Properties.VariableNames(2:end);
bulk_counts=table2array(bulk(:,2:end));

[a,b]=ismember(tissue_cells,neuron_cells);

neuron_counts(:,b(a==1))=[];
neuron_cells(b(a==1))=[];

b=find(strcmpi(tissue_cells,'neuron'));

tissue_counts(:,b)=[];
tissue_cells(b)=[];

[a,b]=ismember(tissue_genes,neuron_genes);

neuron_counts=neuron_counts(b(a==1),:);
neuron_genes=neuron_genes(b(a==1));

[a,b]=ismember(tissue_genes,bulk_genes);

bulk_counts=bulk_counts(b(a==1),:);
bulk_genes=bulk_genes(b(a==1));

r=randperm(length(tissue_genes));
disp('Check gene line up');
[['tissue';tissue_genes(r(1:20))] ['bulk';bulk_genes(r(1:20))] ['neuron';neuron_genes(r(1:20))]]



%% Put things into matrices:
Y=bulk_counts;Y_names=bulk_cells;
X=[neuron_counts tissue_counts];
X_names=[neuron_cells tissue_cells];
X_istissue=[zeros(1,length(neuron_cells)) ones(1,length(tissue_cells))];



%% Determine what to impute
S=double(X>imputation_cutoff); %% Good counts


%% log transform
Y=log1p(Y); % log - counts (bulk)
X=log1p(X); % log - counts (single)



%% find correspondence between bulk and single cell variables
A=char(Y_names');
B=isstrprop(A(:,1:4),'lower');
for i=1:size(A,1)
    [~,idx]=find(isstrprop(A(i,1:4),'lower')==1);
    Y_names_short{i}=A(i,1:(idx(1)-1));
end

for i=1:length(Y_names_short)
    if strcmpi(Y_names_short{i},'DD');
        match{i}=find(strcmpi(X_names,'VD_DD'));
    elseif strcmpi(Y_names_short{i},'IL2');
        match{i}=[find(strcmpi(X_names,'IL2_LR')) find(strcmpi(X_names,'IL2_DV'))];
    elseif strcmpi(Y_names_short{i},'VD');
        match{i}=find(strcmpi(X_names,'VD_DD'));
    else
        match{i}=find(strcmpi(X_names,Y_names_short{i}));
    end
end

M=zeros(length(X_names),length(Y_names));

for i=1:length(Y_names)
    M(match{i},i)=1;
end



%% Joint deconvolution and imputation - MAIN ROUTINE (algorithmic details inside the cons_nmf.m code)

[Z,U,Sigma,obj]=cons_nmf(Y,X,S,M,lambda,lambda2,lambda3,lambda4,hard_constraint,covariance_correction,iter);


%% Bulk subtraction
Y_clean = max(Y-Z*Sigma*(U-M.*U),0); %% subtract the contributions of all cells from bulk except the intended cell, threshold over zero for nicer visualization.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% visualizations

%% Figure - Bulk proportions (U variable)
set(0,'defaulttextInterpreter','none');
figure('units','normalized','outerposition',[0 0 1 1])
imagesc((Sigma*U)');colormap(othercolor('Blues9'));set(gca,'ytick',(1:length(Y_names)),'YTickLabel',Y_names,'xtick',(1:length(X_names)),'XTickLabel',X_names,'XTickLabelRotation',90,'FontWeight','bold');title('Bulk proportions');
hold on
for i=1:size(U,2)
    idx=find(M(:,i));
    for j=1:length(idx)
        plot(idx(j),i,'r.','MarkerSize',8);
    end
end
legend({'targets'});
try;export_fig('proportions.png');end

%% Figure - Imputed single cell reads
figure('units','normalized','outerposition',[0 0 1 1])
imagesc(log1p(Z));colormap(othercolor('Blues9'));set(gca,'xtick',(1:length(X_names)),'XTickLabel',X_names,'XTickLabelRotation',90,'FontWeight','bold');ylabel('Genes');title('Imputed Genes');
try;export_fig('imputed_single.png');end

%% Figure - Mean sample distance metric
figure('units','normalized','outerposition',[0 0 1 1])
X_names_double=[X_names X_names];
imagesc(log1p([Z X]));colormap(othercolor('Blues9'));set(gca,'xtick',(1:length(X_names_double)),'XTickLabel',[X_names_double],'XTickLabelRotation',90,'FontWeight','bold');ylabel('Genes');title('Imputed Genes / Actual genes');
try;export_fig('imputed_single_and_actual_single.png');end


%% Figure - imputation counts scatter vs actual
figure('units','normalized','outerposition',[0 0 1 1])
scatter(X(:),Z(:),'.');xlabel('Actual gene counts');ylabel('Imputed gene counts');
hold on
ux=unique(X(:));
for i=1:100
    plot(ux(i),quantile(Z(X==ux(i)),0.9),'ro');
end
legend({'Counts','Means'});
try;export_fig('imputed_vs_actual_scatter.png');end

%% Figure - Bulk - Single - Single (imputed) - Bulk (cleaned) next to each other
figure('units','normalized','outerposition',[0 0 1 1])
imagesc([log1p([bulk_counts tissue_counts neuron_counts]) Z Y_clean],[0 max(X(:))]);colormap(othercolor('Blues9'));
hold on
plot([size(bulk_counts,2)+0.5 size(bulk_counts,2)+0.5],[0.5 size(bulk_counts,1)+0.5],'r-','LineWidth',2);
plot([size([bulk_counts tissue_counts neuron_counts],2)+0.5 size([bulk_counts tissue_counts neuron_counts],2)+0.5],[0.5 size(bulk_counts,1)+0.5],'r-','LineWidth',2);
plot([size([bulk_counts tissue_counts neuron_counts Z],2)+0.5 size([bulk_counts tissue_counts neuron_counts Z],2)+0.5],[0.5 size(bulk_counts,1)+0.5],'r-','LineWidth',2);
text(10,1000,'Bulk','Color','k','FontWeight','bold','FontSize',15);
text(10+size(bulk_counts,2),1000,'Single','Color','k','FontWeight','bold','FontSize',15);
text(10+size([bulk_counts tissue_counts neuron_counts],2),1000,'Single (imputed)','Color','k','FontWeight','bold','FontSize',15);
text(10+size([bulk_counts tissue_counts neuron_counts Z],2),1000,'Bulk (cleaned)','Color','k','FontWeight','bold','FontSize',15);
set(gca,'xtick',(1:size(log1p([bulk_counts tissue_counts neuron_counts Z Y_clean]),2)),'XTickLabel',[bulk_cells tissue_cells neuron_cells X_names Y_names],'XTickLabelRotation',90);ylabel('genes');
set(gca,'TickLength',[0 0]);
try;export_fig('family_picture.png');end



%% Figure - Multidimensional scaling visualization before/after deconv.
figure('units','normalized','outerposition',[0 0 1 1])
K=Y_clean'*Y_clean;K=K./sqrt(diag(K)*diag(K)');
Y2=mdscale(K,2);

K0=Y'*Y;K0=K0./sqrt(diag(K0)*diag(K0)');
Y20=mdscale(K0,2);
minmax = @(x)((x-min(x(:)))./max(x(:)-min(x(:))));
Y20=minmax(Y20);
Y2=minmax(Y2);
colors=hsv(size(M,1));
for i=1:size(Y2,1)
    idx=find(M(:,i));
    subplot(1,2,2)
    hold on
    plot(Y2(i,1),Y2(i,2),'o','MarkerFaceColor',colors(idx(1),:),'MarkerEdgeColor','k','MarkerSize',5);
    text(Y2(i,1)-max(Y2(:))/100,Y2(i,2)-max(Y2(:))/100,Y_names_short{i},'Color','k','FontWeight','bold')
    axis square;grid on
    subplot(1,2,1)
    hold on
    plot(Y20(i,1),Y20(i,2),'o','MarkerFaceColor',colors(idx(1),:),'MarkerEdgeColor','k','MarkerSize',5);
    text(Y20(i,1)-max(Y2(:))/100,Y20(i,2)-max(Y2(:))/100,Y_names_short{i},'Color','k','FontWeight','bold')
    axis square;grid on
end
subplot(1,2,1);
title('Before deconv + subtraction bulk sample correlation MDS');
subplot(1,2,2);
title('After deconv + subtraction bulk sample correlation MDS');
try;export_fig('bulk_dispersion.png');end



%% Figure - Mean sample distance metric
clear B C
for i=1:size(M,1)
    idx=find(M(i,:));
    B(i,1)=0;
    B(i,2)=0;
    
    C(i,1)=0;
    C(i,2)=0;
    if length(idx)>1
        
        for j=1:length(idx)
            B(i,1)=B(i,1)+norm(mean(Y2(idx,:),1)-Y2(idx(j),:),'fro').^2 / length(idx);
            B(i,2)=B(i,2)+norm(mean(Y20(idx,:),1)-Y20(idx(j),:),'fro').^2 / length(idx);
            C(i,1)=C(i,1)+corr(mean(Y_clean(:,idx),2),Y_clean(:,idx(j)))/ length(idx);
            C(i,2)=C(i,2)+corr(mean(Y(:,idx),2),Y(:,idx(j)))/ length(idx);
        end
    end
end

figure('units','normalized','outerposition',[0 0 1 1])
hold on
for i=1:size(M,1)
    if and(B(i,1)~=0,B(i,2)~=0);
        plot(B(i,1),B(i,2),'o','MarkerFaceColor',colors(i,:),'MarkerEdgeColor','k','MarkerSize',5);
        text(B(i,1)-max(B(:))/100,B(i,2)-max(B(:))/100,X_names{i},'Color','k','FontWeight','bold');
    end
end
hold on;plot([min(B(:)) max(B(:))],[min(B(:)) max(B(:))],'k--','LineWidth',2);xlabel('After deconv correlation');ylabel('Before deconv correlation');set(gca,'FontWeight','bold');
axis equal;axis square;grid on
title('Within sample distance before/after deconv+subtraction');
try;export_fig('alexis_metric_1a.png');end


%% Figure - Mean sample correlation metric
figure('units','normalized','outerposition',[0 0 1 1])
hold on
for i=1:size(M,1)
    if and(C(i,1)~=0,C(i,2)~=0);
        plot(C(i,1),C(i,2),'o','MarkerFaceColor',colors(i,:),'MarkerEdgeColor','k','MarkerSize',5);
        text(C(i,1)-max(C(:))/1000,C(i,2)-max(C(:))/1000,X_names{i},'Color','k','FontWeight','bold');
    end
end
hold on;plot([min(C(C~=0)) max(C(C~=0))],[min(C(C~=0)) max(C(C~=0))],'k--','LineWidth',2);xlabel('After deconv correlation');ylabel('Before deconv correlation');set(gca,'FontWeight','bold');
axis equal;axis square;grid on
title('Within sample correlation before/after deconv+subtraction');
try;export_fig('alexis_metric_1b.png');end


%% Figure -  Correlation matrices before and after deconv.
figure('units','normalized','outerposition',[0 0 1 1])
imagesc([corr(Y(:,:))],[0 1]);axis square;colormap(othercolor('BuDRd_12'));colorbar;set(gca,'Xtick',(1:length(Y_names)),'Xticklabel',Y_names_short,'XTickLabelRotation',90,'ticklength',[0 0],'Ytick',(1:length(Y_names)),'Yticklabel',Y_names_short,'FontWeight','bold');title('Before deconv + subtraction');
try;export_fig('before_correlation.png');end
figure('units','normalized','outerposition',[0 0 1 1])
imagesc([corr(Y_clean(:,:))],[0 1]);axis square;colormap(othercolor('BuDRd_12'));colorbar;set(gca,'Xtick',(1:length(Y_names)),'Xticklabel',Y_names_short,'XTickLabelRotation',90,'ticklength',[0 0],'Ytick',(1:length(Y_names)),'Yticklabel',Y_names_short,'FontWeight','bold');title('After deconv + subtraction');
try;export_fig('after_correlation.png');end


%% Figure -  Maximum margin metric

S=double((M'*M)>0)-eye(size(M'*M));
Sp=(1-S) -eye(size(M'*M));
S(S==0)=nan;
Sp(Sp==0)=nan;

figure('units','normalized','outerposition',[0 0 1 1])
plot(vec(nanmin(corr(Y_clean).*S,[],2)-nanmax(corr(Y_clean).*Sp,[],2)),vec(nanmin(corr(Y).*S,[],2)-nanmax(corr(Y).*Sp,[],2)),'k.','MarkerSize',10);
tmp=[vec(nanmin(corr(Y_clean).*S,[],2)-nanmax(corr(Y_clean).*Sp,[],2)),vec(nanmin(corr(Y).*S,[],2)-nanmax(corr(Y).*Sp,[],2))];
hold on;plot([min(tmp(:)) max(tmp(:))],[min(tmp(:)) max(tmp(:))],'k--','LineWidth',2);xlabel('After deconv');ylabel('Before deconv');set(gca,'FontWeight','bold');
text(tmp(:,1),tmp(:,2),Y_names_short,'Color','k','FontWeight','bold');
axis equal;axis square;grid on
title('Minimum within sample correlation - maximum out of sample correlation before/after deconv+subtraction');
try;export_fig('maximum_margin_metric.png');end

%% Figure -  Mean margin metric
figure('units','normalized','outerposition',[0 0 1 1])
plot(vec(nanmean(corr(Y_clean).*S,2)-nanmean(corr(Y_clean).*Sp,2)),vec(nanmean(corr(Y).*S,2)-nanmean(corr(Y).*Sp,2)),'k.','MarkerSize',10);
tmp=[vec(nanmean(corr(Y_clean).*S,2)-nanmean(corr(Y_clean).*Sp,2)),vec(nanmean(corr(Y).*S,2)-nanmean(corr(Y).*Sp,2))];
hold on;plot([min(tmp(:)) max(tmp(:))],[min(tmp(:)) max(tmp(:))],'k--','LineWidth',2);xlabel('After deconv');ylabel('Before deconv');set(gca,'FontWeight','bold');
text(tmp(:,1),tmp(:,2),Y_names_short,'Color','k','FontWeight','bold');
axis equal;axis square;grid on
title('Mean within sample correlation - Mean out of sample correlation before/after deconv+subtraction');
try;export_fig('mean_margin_metric.png');end