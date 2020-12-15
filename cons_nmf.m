function [Z,U,C,obj]=cons_nmf(X,Y,S,M,lambda_1,lambda2,lambda3,lambda4,hard_constraint,covariance_correction,proportion_normalization,iter)

%% Constrained NMF -
%% Inputs: X - Bulk reads (genes x cells) [non-negative continuous]
%%         Y - Single cell reads (genes x cells) [non-negative continuous]
%%         S - Confident single cell reads (genes x cells) [0/1 binary]
%%         M - Single cell to bulk correspondence (cells x cells) [0/1 binary]
%%         lambda_1 - parameter to control single cell to bulk proportion sparsity (scalar) [non-negative continuous]
%%         lambda_2 - parameter to control imputation sparsity (scalar) [non-negative continuous]
%%         lambda_3 - parameter to control the purity of the single cell to bulk proportions  (scalar) [non-negative continuous]
%%         lambda_4 - parameter to control goodness of fit of the non-imputed single cell reads (scalar) [non-negative continuous]
%%         hard_constraint - switch to enforce hard constraint that confident single cell reads will remain the same as inputs (1) or will be imputed as well (0) [0/1 binary]
%%         covariance_correction - switch to utilize covariance structure of the single cell reads to decorrelate deconvolution proportions (0/1) [0/1 binary]
%%         proportion_normalization - switch to normalize proportions, 0 for off, 1 for rows, 2 for columns (0/1/2) 
%%         iter - iterations to run the algorithm (scalar) [non-negative integer]
%% Outputs:Z - "pure" single cell profiles (genes x cells) [non-negative continuous]
%%         U - single to bulk proportions (cells x cells) [non-negative continuous]
%%         C - single cell covariance estimate (cells x cells) [continuous]
%%         obj - objective function (iterations x 1) [non-negative continuous]
%%
%% This algorithm solves the problem: ||X - Z*U ||_F^2 + lambda_1 * ||U||_1 + lambda_2*||Z o (1-S)||_F^2 + lambda_3*||M o (U-M)||_F^2 + lambda_4*||S o (Y - Z*V) ||_F^2
%% Overall, the algorithm tries to deconvolve the bulk data into "pure" factors that capture single cell profiles and at the same time impute the single cell data such that the single cell reads match the "pure" profiles.
%% Here, the bulk data,X, is factorized into two terms, Z and U, where Z are fully captured single cell gene expression profiles and U are proportions of these that make up the bulk.
%% Simultaneously, the single cell data is imputed such that the non-confident reads (reads such that S==0) are going to be close to Z.
%% There are additional terms to enforce sparsity and purity of the terms.
%% -EV


visual=1;

%% initial covariance estimate (if switch is on)
if covariance_correction==1
    C=cov(Y);
else
    C=eye(size(Y,2));
end

%% initalization of terms
V=eye(size(Y,2));
Z=rand(size(Y));
U=rand(size(Y,2),size(X,2));


figure
for t=1:iter
    Z = (Z.*(X*U'*C' + lambda4*(Y.*S)*V'))./(Z*C*U*U'*C' + lambda4*((Z*V).*S)*V' + lambda2*(Z.*(1-S)) + realmin);%% Multiplicative update for Z
    if hard_constraint==1
        Z(S==1)=Y(S==1);
    end
    %     Z=Z./sum(Z,1); %% normalization that doens't seem to be necessary.
    U = (U.*(C'*Z'*X + lambda3*M))./(C'*Z'*Z*C*U + lambda3*(M.*U) + lambda_1 + realmin); %% Multiplicative update for U
    if proportion_normalization>0
        U=U./sum(U,proportion_normalization); %% normalization that doesn't seem to be necessary.
    end
    if covariance_correction==1
        C=cov(Z);
    else
        C=eye(size(Y,2));
    end
    obj(t,1)=norm(X-Z*C*U,'fro')^2 + lambda4*norm(S.*(Y-Z*V),'fro')^2 + lambda_1*sum(U(:)) + lambda2*norm((1-S).*Z,'fro')^2  + lambda3*norm(M.*(U-M),'fro')^2; %% objective update
    if visual==1
        plot(obj./max(obj,[],1));xlabel('iterations');ylabel('objective');title('Constrained NMF');set(gca,'FontWeight','bold');grid on
        drawnow
    end
end



