function[Outs]=ppl_ADMM_solver_8(MSI,HSI,opts)

if ~exist('opts', 'var')
    opts = [];
end    
if isfield(opts, 'nonzero_indices');     nonzero_indices = opts.nonzero_indices;      end  %% it is a set of matrices
if isfield(opts, 'tol');                 tol             = opts.tol;                  end
if isfield(opts, 'max_iter');            max_iter        = opts.max_iter;	          end
if isfield(opts, 'lambda1');             lambda1   	     = opts.lambda1;   	          end
if isfield(opts, 'lambda2');             lambda2   	     = opts.lambda2;   	          end
if isfield(opts, 'beta');                beta            = opts.beta;                 end
if isfield(opts, 'mu');                  mu       	     = opts.mu;                   end
if isfield(opts, 'S_initial');           S_initial       = opts.S_initial;            end
if isfield(opts, 'S_ideal');             S_ideal         = opts.S_ideal;              end
if isfield(opts, 'D1_spa');              D1_spa          = opts.D1_spa;               end
if isfield(opts, 'D2_spa');              D2_spa          = opts.D2_spa;               end
if isfield(opts, 'Dm_spe');              Dm_spe          = opts.Dm_spe;               end  %% it is a set of matrices
if isfield(opts, 'P1');                  P1              = opts.P1;     	          end
if isfield(opts, 'P2');                  P2              = opts.P2;     	          end
if isfield(opts, 'Pm');                  Pm              = opts.Pm;     	          end
if isfield(opts, 'b1');                  b1              = opts.b1;     	          end
if isfield(opts, 'b2');                  b2              = opts.b2;     	          end
if isfield(opts, 'b31');                 b31             = opts.b31;     	          end
if isfield(opts, 'b32');                 b32             = opts.b32;     	          end
if isfield(opts, 'b33');                 b33             = opts.b33;     	          end
if isfield(opts, 'b34');                 b34             = opts.b34;     	          end
if isfield(opts, 'b35');                 b35             = opts.b35;     	          end
if isfield(opts, 'b36');                 b36             = opts.b36;     	          end
if isfield(opts, 'b37');                 b37             = opts.b37;     	          end
if isfield(opts, 'b38');                 b38             = opts.b38;     	          end
if isfield(opts, 'UU');                  UU              = opts.UU;     	          end

dim         = size(S_initial);
S           = S_initial;
Z           = S;
Multiplier  = zeros(size(S));


%% Preparing for the iteration 
   M3 = tens2mat(MSI,3);
   H2 = tens2mat(HSI,2);
   H1 = tens2mat(HSI,1);
   %%%% spectral response function:      
    Dm1  = Dm_spe{1};
    Dm2  = Dm_spe{2};
    Dm3  = Dm_spe{3};
    Dm4  = Dm_spe{4};
    Dm5  = Dm_spe{5};
    Dm6  = Dm_spe{6};
    Dm7  = Dm_spe{7};
    Dm8  = Dm_spe{8};
    %%% row extraction vector:
    dm1  = row_extraction(b31);
    dm2  = row_extraction(b32);
    dm3  = row_extraction(b33);
    dm4  = row_extraction(b34);
    dm5  = row_extraction(b35);
    dm6  = row_extraction(b36);
    dm7  = row_extraction(b37);
    dm8  = row_extraction(b38);
    %%% nonzero_indices:
    nonzero_indices_m1 = nonzero_indices{1};
    nonzero_indices_m2 = nonzero_indices{2};
    nonzero_indices_m3 = nonzero_indices{3};
    nonzero_indices_m4 = nonzero_indices{4};
    nonzero_indices_m5 = nonzero_indices{5};
    nonzero_indices_m6 = nonzero_indices{6};
    nonzero_indices_m7 = nonzero_indices{7};
    nonzero_indices_m8 = nonzero_indices{8};
    %%% give the space of blur kernel to generate Pm:   
    pm1_long  = sum(Dm1,1);
    pm2_long  = sum(Dm2,1);
    pm3_long  = sum(Dm3,1);
    pm4_long  = sum(Dm4,1);
    pm5_long  = sum(Dm5,1);
    pm6_long  = sum(Dm6,1);
    pm7_long  = sum(Dm7,1);
    pm8_long  = sum(Dm8,1);

%% Run LADMM
 fprintf('============= Run Our ppl_ADMM =============\n');
 t1 = tic;
for iter = 1:max_iter
    Sk  = S;
    P1k = P1;
    P2k = P2;
    Pmk = Pm; 
    Zk  = Z;
    b1k = b1;
    b2k = b2;
    
    b31k = b31;
    b32k = b32;
    b33k = b33;
    b34k = b34;
    b35k = b35;
    b36k = b36;
    b37k = b37;
    b38k = b38;
    
   %% Updating S subproblem: [I,J,K]
   
    Hp1 = P1k'*P1k; 
    Hp2 = P2k'*P2k; 
    Hpm = Pmk'*Pmk;
    nabla1  =  tmprod(tmprod(Sk,Hp1,1),Hp2,2)-tmprod(tmprod(HSI,P1k',1),P2k',2);
    nabla2  =  tmprod(Sk,Hpm,3)-tmprod(MSI,Pmk',3);
    Nabla   = nabla1 + lambda1*nabla2;

    
    Ls      = max(1,norm_sparse(Hp1)*norm_sparse(Hp2)+lambda1*norm_sparse(Hpm));
    alpha   = 1/Ls;
    S1      = Sk - (Nabla+Multiplier)/alpha;
    S2      = Zk;
    C       =  (alpha/(alpha+beta))*S1+(beta/(alpha+beta))*S2;
    [S,tnn, trank(iter)] = prox_utnn(UU,C,lambda2./(alpha+beta));
    loss_f3 = lambda2*tnn;
   %% Updating Z subproblem: [I,J,K]
    Nu1     = S + Multiplier/beta;
    Nu12    = (beta./(beta+mu))*Nu1+(mu./(beta+mu))*Zk;
    Z       = Nu12;
    Z(Z<0)  = 0;   
   %% Updating the multipliers
    Multiplier_old = Multiplier;
    Multiplier = Multiplier + beta*(S-Z);  
   %% Updating P1 subproblem: [I',I]
    NSI_1 = tmprod(S,P2k,2);
    N1 = tens2mat(NSI_1,1);
    Nn = N1;
    Hn = H1;
    [b1,P1] = b_solver_L(Nn,Hn,D1_spa,b1k);
    
   %% Updating P2 subproblem: [J',J]
    NSI_2    = tmprod(S,P1,1);
    N2       = tens2mat(NSI_2,2);
    Nnn      = N2;
    Hnn      = H2;
    [b2,P2]  = b_solver_L(Nnn,Hnn,D2_spa,b2k);
    
   %% Updating Pm subproblem: [K',K]
    %%%% updating 1st band:
    S1        = tmprod(S,Dm1,3);
    S1n       = tens2mat(S1,3);
    S1nn      = S1n;
    M1n       = MSI(:,:,1);
    M1nn      = reshape(M1n,1,[]);
    [b31,pm1] = bm_solver_L(S1nn,M1nn,dm1,b31k);
    
    %%% updating 2nd band:
    S2         = tmprod(S,Dm2,3);
    S2n        = tens2mat(S2,3);
    S2nn       = S2n;
    M2n        = MSI(:,:,2);
    M2nn       = reshape(M2n,1,[]);
    [b32,pm2]  = bm_solver_L(S2nn,M2nn,dm2,b32k);
      
    %%% updating 3rd band:
    S3         = tmprod(S,Dm3,3);
    S3n        = tens2mat(S3,3);
    S3nn       = S3n;
    M3n        = MSI(:,:,3);
    M3nn       = reshape(M3n,1,[]);
    [b33,pm3]  = bm_solver_L(S3nn,M3nn,dm3,b33k);
    
    %%% updating 4th band:
    S4         = tmprod(S,Dm4,3);
    S4n        = tens2mat(S4,3);
    S4nn       = S4n;
    M4n        = MSI(:,:,4);
    M4nn       = reshape(M4n,1,[]);
    [b34,pm4]  = bm_solver_L(S4nn,M4nn,dm4,b34k);   
    %%% updating 5th band:
    S5         = tmprod(S,Dm5,3);
    S5n        = tens2mat(S5,3);
    S5nn       = S5n;
    M5n        = MSI(:,:,5);
    M5nn       = reshape(M5n,1,[]);
    [b35,pm5]  = bm_solver_L(S5nn,M5nn,dm5,b35k);
   
    %%% updating 6th band:
    S6         = tmprod(S,Dm6,3);
    S6n        = tens2mat(S6,3);
    S6nn       = S6n;
    M6n        = MSI(:,:,6);
    M6nn       = reshape(M6n,1,[]);
    [b36,pm6]  = bm_solver_L(S6nn,M6nn,dm6,b36k);
    
    %%% updating 7th band:
    S7         = tmprod(S,Dm7,3);
    S7n        = tens2mat(S7,3);
    S7nn       = S7n;
    M7n        = MSI(:,:,7);
    M7nn       = reshape(M7n,1,[]);
    [b37,pm7]  = bm_solver_L(S7nn,M7nn,dm7,b37k);
    
        %%% updating 8th band:
    S8         = tmprod(S,Dm8,3);
    S8n        = tens2mat(S8,3);
    S8nn       = S8n;
    M8n        = MSI(:,:,8);
    M8nn       = reshape(M8n,1,[]);
    [b38,pm8]  = bm_solver_L(S8nn,M8nn,dm8,b38k);
    
 
    
    %%% Upating pm:
    pm1_long(nonzero_indices_m1) = pm1;
    pm2_long(nonzero_indices_m2) = pm2;
    pm3_long(nonzero_indices_m3) = pm3;
    pm4_long(nonzero_indices_m4) = pm4;
    pm5_long(nonzero_indices_m5) = pm5;
    pm6_long(nonzero_indices_m6) = pm6;
    pm7_long(nonzero_indices_m7) = pm7;
    pm8_long(nonzero_indices_m8) = pm8;
    %%% Updating Pm matrix:
    Pm(1,:) = pm1_long;
    Pm(2,:) = pm2_long;
    Pm(3,:) = pm3_long;
    Pm(4,:) = pm4_long;
    Pm(5,:) = pm5_long;
    Pm(6,:) = pm6_long;
    Pm(7,:) = pm7_long;
    Pm(8,:) = pm8_long;
   %% Computing stop criterion   
   %%% for S:   
 %%% [S_inter,~, ~] = prox_our_tnn(S1,lambda2/alpha);
 %%5  eta_S(iter) = norm(S_inter(:)-S(:))/(1+norm(Multiplier_old(:))+norm(Sk(:)));
   %%% for Z 
 %%5  eta_Z(iter) = norm(S(:)-Z(:))/(1+norm(Z(:))+norm(S(:)));
 %%%  eta(iter)   = max(eta_S(iter),eta_Z(iter));
    %% fprintf all results
   Loss1(iter) = 0.5*frob(HSI-tmprod(tmprod(S,P1,1),P2,2))^2;
   Relerr1(iter) = (frob(HSI-tmprod(tmprod(S,P1,1),P2,2))^2)/(frob(HSI)^2);
   Loss2(iter) = (lambda1/2)*frob(MSI-tmprod(S,Pm,3))^2;
   Relerr2(iter) = (frob(MSI-tmprod(S,Pm,3))^2)/(frob(MSI)^2);
   Loss3(iter) = loss_f3;
   Aobj(iter)  = Loss1(iter)+ Loss2(iter)+ Loss3(iter)+inprod(Multiplier,S-Z)+(beta/2)*frob(S-Z)^2;
   Obj(iter)   = Loss1(iter)+ Loss2(iter)+ Loss3(iter);
      Psnr(iter) = PSNR(S,S_ideal);

   Relative_error(iter) = Relerr1(iter)+ lambda1*Relerr2(iter);
     if iter == 1 || mod(iter, 1) == 0
        fprintf('iter = %3.d, PSNR = %3.2f, relative error = %.4f, Augment Largrange obj = %.1f, Ori_obj = %.1f \n',iter,Psnr(iter),Relative_error(iter),Aobj(iter),Obj(iter));
     end  
     

   iter = iter+1;  
  % beta = 1.1*beta;
end
   t2 = toc(t1);
   fprintf('Our ppl_ADMM took %2.2f seconds to run.\n', t2);
    Outs.P1   = P1;
    Outs.P2   = P2;
    Outs.Pm   = Pm;
    Outs.SRI  = S;
    Outs.Z    = Z;
    Outs.Aobj = Aobj';
    Outs.Obj  = Obj';
    Outs.Relrr = Relative_error';
    Outs.Trank = trank';
 %   Outs.eta = eta';
 %   Outs.eta_S = eta_S';
 %   Outs.eta_Z = eta_Z';
    Outs.blur1_spa = b1;
    Outs.blur2_spa = b2;
    Outs.blur31_spe = b31;
    Outs.blur32_spe = b32;
    Outs.blur33_spe = b33;
    Outs.blur34_spe = b34;
    Outs.blur35_spe = b35;
    Outs.blur36_spe = b36;
    Outs.blur37_spe = b37;
    Outs.blur38_spe = b38;
end
    
    






