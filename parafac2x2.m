function [F,A,D,pvar,ssr] = parafac2x2(varargin)
%[F,A,D,pvar,ssr,timeOut] =
%parafac2x2(Xijkl,R,eps,maxiter,displ,animate,Bkli,Dkli,Akli,Bskli,Bili,Dili,Aili,Bsili)
%
%INPUTS
%Xijkl := 4 way tensor with drift across i and k.
%R := number of factors
%eps := convergence criterion (optional, default 1e-8)
%maxiter := maximum number of iterations (optional, default 1000)
%displ := [0,1] toggles the command line readout o
%animate := [0,1] toggles the animation display
%Bkli := \in R^[i,R,k*l] initialisation of the elution profiles of the ith mode
%Dkli := \in R^[R,R,k*l] initialisation of the quantitative components
%along the ith mode.
%Akli := \in R^[j,R] initialisation of the mass spectra along the klth
%intermediate of the data.
%Bskli := \in R^[j,R] initialisation of the Bs coupling matrix along the
%klth intermediate of the data.
%Bili := \in R^[k,R,i*l] initialisation of the elution profiles of the kth
%mode.
%Dili := \in R^[R,R,i*l] initialisation of the quantitative components
%along the kth mode.
%Aili := \in R^[j,R] initialisations of the mass spectra along the ilth
%intermediate of the data.
%Bsili := \in R^[R,R] initialisations of the coupling matrix along the ilth
%intermediate of the data.
%
%OUTPUTS
%F := \in R^[i,R,k,l] unique 2-dimensional score profiles for each sample
%A := \in R^[j,R] mass spectra
%D := \in R^[L,R] quantitative loadings
%pvar := percent variance explained in the model relative to the raw data
%ssr := descent of the relative sum of residual squares. If not decreasing at each iteration,
%indicates a problem with the ALS implementation of the algorithm.
%
%Algorithm for multi-linear decomposition of GCxGC-TOFMS data, based on
%alternating 2tr, 1tr frontal slabs to solve for elution scores. The two models are coupled via the mass spectral information which is used to solve for the sample loadings. Agnostic to
%drift in the first and second dimensions.
%
%The tensor must be organised properly for this algorithm to work. In case
%of failure, visually ensure the proper organisation of the modes according
%to: I X J X K X L where I are mass spectral acquisitions, J are the mass
%spectral ion channels, K are the modulations (per sample), and L are the
%individual samples themselves.
%
%In its current implementation, the only acceptable input is a tenosr, and
%not a multiway array. This will hopefully change in a future update, to
%handle more flexible regions of interest for chromatographic data.
%
%Bkli,Dkli,Akli,Bskli,Bili,Dili,Aili,Bsili are initialisation parameters.
%These are determined automatically if not otherwise indicated.
%
%(c) Michael Sorochan Armstrong, James J. Harynuk 2022
%
%This software has been stress-tested, and it is the hope of the developers
%that it works for its intended purposes. If it does not, please email the
%corresponding author: mdarmstr@ualberta.ca


if size(varargin,2) < 2
    error('Must declare input tensor, number of factors');
else
end

Xijkl = varargin{1};
R = varargin{2};
sz = size(Xijkl);

%Variable user input
switch size(varargin,2)
    case 14
        
        displ = varargin{5};
        
        if displ == 1
            disp('Utilising input initialisations for Bkl,Dkl,Akl,Bskl,Bil,Dil,Ail,Bsil');
        else
        end
        
        eps = varargin{3};
        maxiter = varargin{4};
        displ = varargin{5};
        animate = varargin{6};
        Bkli = varargin{7};
        Dkli = varargin{8};
        Akli = varargin{9};
        Bskli = varargin{10};
        Bili = varargin{11};
        Dili = varargin{12};
        Aili = varargin{13};
        Bsili = varargin{14};
        [F,A,D,pvar,ssr] = nnparafac2x2als(Xijkl,R,eps,maxiter,displ,animate,Bkli,Dkli,Akli,Bskli,Bili,Dili,Aili,Bsili);
    case 13
        displ = varargin{5};
        
        if displ == 1
            disp('Utilising input initialisations for Bkl,Dkl,Akl,Bskl,Bil,Dil,Ail');
        else
        end
        
        eps = varargin{3};
        maxiter = varargin{4};
        animate = varargin{6};
        Bsili = rand(R,R);
        Bkli = varargin{7};
        Dkli = varargin{8};
        Akli = varargin{9};
        Bskli = varargin{10};
        Bili = varargin{11};
        Dili = varargin{12};
        Aili = varargin{13};
        [F,A,D,pvar,ssr] = nnparafac2x2als(Xijkl,R,eps,maxiter,displ,animate,Bkli,Dkli,Akli,Bskli,Bili,Dili,Aili,Bsili);
    case 12
        displ = varargin{5};
        
        if displ == 1
            disp('Utilising input initialisations for Bkl,Dkl,Akl,Bskl,Bil,Dil');
        else
        end
        
        eps = varargin{3};
        maxiter = varargin{4};
        animate = varargin{6};
        Bsili = rand(R,R);
        Aili = rand(sz(2),R);
        Aili = Aili./vecnorm(Aili,2);
        Bkli = varargin{7};
        Dkli = varargin{8};
        Akli = varargin{9};
        Bskli = varargin{10};
        Bili = varargin{11};
        Dili = varargin{12};
        [F,A,D,pvar,ssr] = nnparafac2x2als(Xijkl,R,eps,maxiter,displ,animate,Bkli,Dkli,Akli,Bskli,Bili,Dili,Aili,Bsili);
    case 11
        displ = varargin{5};
        
        if displ ==1
            disp('Utilising input initialisations for Bkl,Dkl,Akl,Bskl,Bil');
        else
        end
        
        eps = varargin{3};
        maxiter = varargin{4};
        displ = varargin{5};
        animate = varargin{6};
        Bsili = rand(R,R);
        Aili = rand(sz(2),R);
        Aili = Aili./vecnorm(Aili,2);
        Dili = repmat(diag(ones(1,R)),[1,1,sz(1)*sz(4)]);
        Bkli = varargin{7};
        Dkli = varargin{8};
        Akli = varargin{9};
        Bskli = varargin{10};
        Bili = varargin{11};
        [F,A,D,pvar,ssr] = nnparafac2x2als(Xijkl,R,eps,maxiter,displ,animate,Bkli,Dkli,Akli,Bskli,Bili,Dili,Aili,Bsili);
    case 10
        displ = varargin{5};
        
        if displ == 1
            disp('Utilising input initialisations for Bkl,Dkl,Akl,Bskl');
        else
        end
        
        eps = varargin{3};
        maxiter = varargin{4};
        animate = varargin{6};
        Bsili = rand(R,R);
        Aili = rand(sz(2),R);
        Aili = Aili./vecnorm(Aili,2);
        Dili = repmat(diag(ones(1,R)),[1,1,sz(1)*sz(4)]);
        Bili = rand(sz(3),R,sz(1)*sz(4));
        Bkli = varargin{7};
        Dkli = varargin{8};
        Akli = varargin{9};
        Bskli = varargin{10};
        [F,A,D,pvar,ssr] = nnparafac2x2als(Xijkl,R,eps,maxiter,displ,animate,Bkli,Dkli,Akli,Bskli,Bili,Dili,Aili,Bsili);
    case 9
        displ = varargin{5};
        
        if displ == 1
            disp('Utilising input intialisations for Bkl, Dkl, Akl');
        else
        end
        
        eps = varargin{3};
        maxiter = varargin{4};
        animate = varargin{6};
        Bsili = rand(R,R);
        Aili = rand(sz(2),R);
        Aili = Aili./vecnorm(Aili,2);
        Dili = repmat(diag(ones(1,R)),[1,1,sz(1)*sz(4)]);
        Bili = rand(sz(3),R,sz(1)*sz(4));
        Bskli = rand(R);
        Bkli = varargin{7};
        Dkli = varargin{8};
        Akli = varargin{9};
        [F,A,D,pvar,ssr] = nnparafac2x2als(Xijkl,R,eps,maxiter,displ,animate,Bkli,Dkli,Akli,Bskli,Bili,Dili,Aili,Bsili);
    case 8
        displ = varargin{5};
        
        if displ == 1
            disp('Utilising input intialisations for Bkl, Dkl');
        else
        end
        
        eps = varargin{3};
        maxiter = varargin{4};
        animate = varargin{6};
        Bsili = rand(R,R);
        Aili = rand(sz(2),R);
        Aili = Aili./vecnorm(Aili,2);
        Dili = repmat(diag(ones(1,R)),[1,1,sz(1)*sz(4)]);
        Bili = rand(sz(3),R,sz(1)*sz(4));
        Bskli = rand(R);
        Akli = rand(sz(2),R);
        Akli = Akli./vecnorm(Akli,2);
        Bkli = varargin{7};
        Dkli = varargin{8};
        [F,A,D,pvar,ssr] = nnparafac2x2als(Xijkl,R,eps,maxiter,displ,animate,Bkli,Dkli,Akli,Bskli,Bili,Dili,Aili,Bsili);
    case 7
        displ = varargin{5};
        
        if displ == 1
            disp('Utilising input intialisation for Bkl');
        else
        end
        
        eps = varargin{3};
        maxiter = varargin{4};
        animate = varargin{6};
        Bsili = rand(R);
        Aili = rand(sz(2),R);
        Aili = Aili./vecnorm(Aili,2);
        Dili = repmat(diag(ones(1,R)),[1,1,sz(1)*sz(4)]);
        Bili = rand(sz(3),R,sz(1)*sz(4));
        Bskli = rand(R);
        Akli = rand(sz(2),R);
        Akli = Akli./vecnorm(Akli,2);
        Dkli = repmat(diag(ones(1,R)),[1,1,sz(3)*sz(4)]);
        Bkli = varargin{7};
        [F,A,D,pvar,ssr] = nnparafac2x2als(Xijkl,R,eps,maxiter,displ,animate,Bkli,Dkli,Akli,Bskli,Bili,Dili,Aili,Bsili);
    case 6
        displ = varargin{5};
        animate = varargin{6};
        eps = varargin{3};
        maxiter = varargin{4};
        
        if displ == 1
            disp('No input initialisations, parameter for maximum iteration number, or convergence threshold provided:');
            disp('Utilising best of 10 random initialisations');
        else
        end
        
        %preallocation of the different initialisation subsets tensors
        Bsili = zeros(R,R,10);
        Aili = zeros(sz(2),R,10);
        Dili = zeros(R,R,sz(1)*sz(4),10);
        Bili = zeros(sz(3),R,sz(1)*sz(4),10);
        Bskli = zeros(R,R,10);
        Akli = zeros(sz(2),R,10);
        Dkli = zeros(R,R,sz(3)*sz(4),10);
        Bkli = zeros(sz(1),R,sz(3)*sz(4));
        
        for re = 1:10
            if displ == 1
                fprintf('Initialisation %d out of 10...',re)
            else
            end
            Bsili(:,:,re) = rand(R,R);
            Bsili(:,:,re) = Bsili(:,:,re)./vecnorm(Bsili(:,:,re),2);
            Ailit = rand(sz(2),R);
            Aili(:,:,re) = Ailit./vecnorm(Ailit,2);
            Dili(:,:,:,re) = repmat(diag(ones(1,R)),[1,1,sz(1)*sz(4)]);
            Bili(:,:,:,re) = rand(sz(3),R,sz(1)*sz(4));
            Bskli(:,:,re) = rand(R,R);
            Bskli(:,:,re) = Bskli(:,:,re)./vecnorm(Bskli(:,:,re),2);
            Aklit = rand(sz(2),R);
            Akli(:,:,re) = Aklit./vecnorm(Aklit,2);
            Dkli(:,:,:,re) = repmat(diag(ones(1,R)),[1,1,sz(3)*sz(4)]);
            Bkli(:,:,:,re) = rand(sz(1),R,sz(3)*sz(4));
            [~,~,~,pvar,ssr(:,re)] = nnparafac2x2als(Xijkl,R,1e-16,15,0,0,Bkli(:,:,:,re),Dkli(:,:,:,re),Akli(:,:,re),Bskli(:,:,re),Bili(:,:,:,re),Dili(:,:,:,re),Aili(:,:,re),Bsili(:,:,re)); %#ok
            if displ == 1
                fprintf('complete \n')
            else
            end
        end
        
        [~,init_indx] = min(ssr(end,:));
        [F,A,D,pvar,ssr] = nnparafac2x2als(Xijkl,R,eps,maxiter,displ,animate,Bkli(:,:,:,init_indx),Dkli(:,:,:,init_indx),Akli(:,:,init_indx),Bskli(:,:,init_indx),Bili(:,:,:,init_indx),Dili(:,:,:,init_indx),Aili(:,:,init_indx),Bsili(:,:,init_indx));
        
    case 5
        displ = varargin{5};
        eps = varargin{3};
        maxiter = varargin{4};
        
        if displ == 1
            disp('No input initialisations, parameter for maximum iteration number, or convergence threshold provided:');
            disp('Utilising best of 10 random initialisations');
            disp('Animations disabled by default');
        else
        end
        
        %preallocation of the different initialisation subsets tensors
        Bsili = zeros(R,R,10);
        Aili = zeros(sz(2),R,10);
        Dili = zeros(R,R,sz(1)*sz(4),10);
        Bili = zeros(sz(3),R,sz(1)*sz(4),10);
        Bskli = zeros(R,R,10);
        Akli = zeros(sz(2),R,10);
        Dkli = zeros(R,R,sz(3)*sz(4),10);
        Bkli = zeros(sz(1),R,sz(3)*sz(4));
        
        for re = 1:10
            if displ == 1
                fprintf('Initialisation %d out of 10...',re)
            else
            end
            Bsili(:,:,re) = rand(R);
            Ailit = rand(sz(2),R);
            Aili(:,:,re) = Ailit./vecnorm(Ailit,2);
            Dili(:,:,:,re) = repmat(diag(ones(1,R)),[1,1,sz(1)*sz(4)]);
            Bili(:,:,:,re) = rand(sz(3),R,sz(1)*sz(4));
            Bskli(:,:,re) = rand(R);
            Aklit = rand(sz(2),R);
            Akli(:,:,re) = Aklit./vecnorm(Aklit,2);
            Dkli(:,:,:,re) = repmat(diag(ones(1,R)),[1,1,sz(3)*sz(4)]);
            Bkli(:,:,:,re) = rand(sz(1),R,sz(3)*sz(4));
            [~,~,~,pvar,ssr(:,re)] = nnparafac2x2als(Xijkl,R,1e-16,15,0,0,Bkli(:,:,:,re),Dkli(:,:,:,re),Akli(:,:,re),Bskli(:,:,re),Bili(:,:,:,re),Dili(:,:,:,re),Aili(:,:,re),Bsili(:,:,re)); %#ok
            if displ == 1
                fprintf('complete \n')
            else
            end
        end
        
        [~,init_indx] = min(ssr(end,:));
        [F,A,D,pvar,ssr] = nnparafac2x2als(Xijkl,R,eps,maxiter,displ,0,Bkli(:,:,:,init_indx),Dkli(:,:,:,init_indx),Akli(:,:,init_indx),Bskli(:,:,init_indx),Bili(:,:,:,init_indx),Dili(:,:,:,init_indx),Aili(:,:,init_indx),Bsili(:,:,init_indx));
        
    case 4
        
        disp('No input initialisations, parameter for maximum iteration number, or convergence threshold provided:');
        disp('Utilising best of 10 random initialisations');
        disp('Animations disabled by default; display enabled by default');
        
        %preallocation of the different initialisation subsets tensors
        Bsili = zeros(R,R,10);
        Aili = zeros(sz(2),R,10);
        Dili = zeros(R,R,sz(1)*sz(4),10);
        Bili = zeros(sz(3),R,sz(1)*sz(4),10);
        Bskli = zeros(R,R,10);
        Akli = zeros(sz(2),R,10);
        Dkli = zeros(R,R,sz(3)*sz(4),10);
        Bkli = zeros(sz(1),R,sz(3)*sz(4));
        
        for re = 1:10
            fprintf('Initialisation %d out of 10...',re)
            Bsili(:,:,re) = rand(R);
            Ailit = rand(sz(2),R);
            Aili(:,:,re) = Ailit./vecnorm(Ailit,2);
            Dili(:,:,:,re) = repmat(diag(ones(1,R)),[1,1,sz(1)*sz(4)]);
            Bili(:,:,:,re) = rand(sz(3),R,sz(1)*sz(4));
            Bskli(:,:,re) = rand(R);
            Aklit = rand(sz(2),R);
            Akli(:,:,re) = Aklit./vecnorm(Aklit,2);
            Dkli(:,:,:,re) = repmat(diag(ones(1,R)),[1,1,sz(3)*sz(4)]);
            Bkli(:,:,:,re) = rand(sz(1),R,sz(3)*sz(4));
            [~,~,~,pvar,ssr(:,re)] = nnparafac2x2als(Xijkl,R,1e-16,15,0,0,Bkli(:,:,:,re),Dkli(:,:,:,re),Akli(:,:,re),Bskli(:,:,re),Bili(:,:,:,re),Dili(:,:,:,re),Aili(:,:,re),Bsili(:,:,re)); %#ok
            fprintf('complete \n')
        end
        
        eps = varargin{3};
        maxiter = varargin{4};
        [~,init_indx] = min(ssr(end,:));
        [F,A,D,pvar,ssr] = nnparafac2x2als(Xijkl,R,eps,maxiter,1,0,Bkli(:,:,:,init_indx),Dkli(:,:,:,init_indx),Akli(:,:,init_indx),Bskli(:,:,init_indx),Bili(:,:,:,init_indx),Dili(:,:,:,init_indx),Aili(:,:,init_indx),Bsili(:,:,init_indx));
        
    case 3
        
        disp('No input initialisations, parameter for maximum iteration number, or convergence threshold provided:');
        disp('Utilising best of 10 random initialisations, maxiter = 1000');
        disp('Animations disabled by default; display enabled by default');
        
        %preallocation of the different initialisation subsets tensors
        Bsili = zeros(R,R,10);
        Aili = zeros(sz(2),R,10);
        Dili = zeros(R,R,sz(1)*sz(4),10);
        Bili = zeros(sz(3),R,sz(1)*sz(4),10);
        Bskli = zeros(R,R,10);
        Akli = zeros(sz(2),R,10);
        Dkli = zeros(R,R,sz(3)*sz(4),10);
        Bkli = zeros(sz(1),R,sz(3)*sz(4));
        
        for re = 1:10
            fprintf('Initialisation %d out of 10...',re)
            Bsili(:,:,re) = rand(R);
            Ailit = rand(sz(2),R,re);
            Aili(:,:,re) = Ailit(:,:,re)./vecnorm(Ailit,2);
            Dili(:,:,:,re) = repmat(diag(ones(1,R)),[1,1,sz(1)*sz(4)]);
            Bili(:,:,:,re) = rand(sz(3),R,sz(1)*sz(4));
            Bskli(:,:,re) = rand(R);
            Aklit = rand(sz(2),R);
            Akli(:,:,re) = Aklit./vecnorm(Aklit,2);
            Dkli(:,:,:,re) = repmat(diag(ones(1,R)),[1,1,sz(3)*sz(4)]);
            Bkli(:,:,:,re) = rand(sz(1),R,sz(3)*sz(4));
            [~,~,~,pvar,ssr(:,re)] = nnparafac2x2als(Xijkl,R,1e-16,15,0,0,Bkli,Dkli,Akli,Bskli,Bili,Dili,Aili,Bsili); %#ok
            fprintf('complete \n')
        end
        
        [~,init_indx] = min(ssr(end,:));
        
        eps = varargin{3};
        [F,A,D,pvar,ssr] = nnparafac2x2als(Xijkl,R,eps,1000,Bkli(:,:,:,init_indx),Dkli(:,:,:,init_indx),Akli(:,:,init_indx),Bskli(:,:,init_indx),Bili(:,:,:,init_indx),Dili(:,:,:,init_indx),Aili(:,:,init_indx),Bsili(:,:,init_indx));
        
    case 2
        disp('No input initialisations, parameter for maximum iteration number, or convergence threshold provided:');
        disp('Utilising best of 10 random initialisations, maxiter = 1000, and eps = 1e-8');
        disp('Animations disabled by default; display enabled by default');
        
        %preallocation of the different initialisation subsets tensors
        Bsili = zeros(R,R,10);
        Aili = zeros(sz(2),R,10);
        Dili = zeros(R,R,sz(1)*sz(4),10);
        Bili = zeros(sz(3),R,sz(1)*sz(4),10);
        Bskli = zeros(R,R,10);
        Akli = zeros(sz(2),R,10);
        Dkli = zeros(R,R,sz(3)*sz(4),10);
        Bkli = zeros(sz(1),R,sz(3)*sz(4));
        
        for re = 1:10
            fprintf('Initialisation %d out of 10...',re)
            Bsili(:,:,re) = rand(R,R);
            Ailit = rand(sz(2),R);
            Aili(:,:,re) = Ailit./vecnorm(Ailit,2);
            Dili(:,:,:,re) = repmat(diag(ones(1,R)),[1,1,sz(1)*sz(4)]);
            Bili(:,:,:,re) = rand(sz(3),R,sz(1)*sz(4));
            Bskli(:,:,re) = rand(R,R);
            Aklit = rand(sz(2),R);
            Akli(:,:,re) = Aklit./vecnorm(Aklit,2);
            Dkli(:,:,:,re) = repmat(diag(ones(1,R)),[1,1,sz(3)*sz(4)]);
            Bkli(:,:,:,re) = rand(sz(1),R,sz(3)*sz(4));
            [~,~,~,pvar,ssr(:,re)] = nnparafac2x2als(Xijkl,R,1e-16,15,0,0,Bkli(:,:,:,re),Dkli(:,:,:,re),Akli(:,:,re),Bskli(:,:,re),Bili(:,:,:,re),Dili(:,:,:,re),Aili(:,:,re),Bsili(:,:,re)); %#ok
            fprintf('complete \n')
        end
        
        [~,init_indx] = min(ssr(end,:));
        [F,A,D,pvar,ssr] = nnparafac2x2als(Xijkl,R,1e-8,1000,1,0,Bkli(:,:,:,init_indx),Dkli(:,:,:,init_indx),Akli(:,:,init_indx),Bskli(:,:,init_indx),Bili(:,:,:,init_indx),Dili(:,:,:,init_indx),Aili(:,:,init_indx),Bsili(:,:,init_indx));
        
end

end


function [F,A,D,pvar,ssr] = nnparafac2x2als(Xijkl,R,eps,maxiter,displ,animate,Bkl,Dkl,Akl,Bskl,Bil,Dil,Ail,Bsil)

%Internal settings: you shouldn't need to change these.
iter = 1;
coupit = 10;
nnls = 3; %1 for fnnls, 2 for fcnnls, 3 for cls
omg = 3;

YNorm = norm(Xijkl(:))^2;

sz = size(Xijkl);

%2tr frontal plane - tensor of observations
Xkl = reshape(Xijkl,[sz(1),sz(2),sz(3)*sz(4)]);

%1tr frontal plane - tensor of observations
Xilp = permute(Xijkl,[3 2 1 4]);
Xil = reshape(Xilp,[sz(3),sz(2),sz(1)*sz(4)]);

mkl = zeros(1,sz(3)*sz(4));
mhkl = mkl;
for kl = 1:sz(3)*sz(4)
    mkl(kl) = norm(Bkl(:,:,kl)*Dkl(:,:,kl)*Akl','fro')./norm(Bkl(:,:,kl),'fro');
end

mil = zeros(1,sz(1)*sz(4));
mhil = mil;
for il = 1:sz(1)*sz(4)
    mil(il) = norm(Bil(:,:,il)*Dil(:,:,il)*Ail','fro')./norm(Bil(:,:,il),'fro');
end

%Estimation of residuals, initial
XHkl = zeros(size(Xkl));
for kl = 1:sz(3)*sz(4)
    XHkl(:,:,kl) = Xkl(:,:,kl) - Bkl(:,:,kl)*(Dkl(:,:,kl))*Akl'; %+ mkl(kl)*norm(Bkl(:,:,kl) - Pkl(:,:,kl)*Bskl,'fro')^2;
end

XHil = zeros(size(Xil));
for il = 1:sz(1)*sz(4)
    XHil(:,:,il) = Xil(:,:,il) - Bil(:,:,il)*(Dil(:,:,il))*Ail';% + mil(il).*norm(Bil(:,:,il) - Pil(:,:,il)*Bsil,'fro')^2;
end

%Mass Spectral coupling constant
mA = (10^(omg))*0.5*(norm(XHkl(:),'fro')^2 + norm(XHil(:),'fro')^2)/norm(Akl,'fro')^2;

ssr1 = norm(XHkl(:)) + norm(XHil(:));
ssr2 = YNorm; %arbitrary value, algorithm will go for at least one iteration

if animate == 1
    h = figure;
    axis tight manual
    filename = strcat(date,'msdescent.gif');
else
end

%more preallocation, kl
Pkl = zeros(sz(1),R,sz(3)*sz(4));
Bstkl = zeros(R,R,sz(3)*sz(4));
BklDkl = Pkl;

%more preallocation, il
Pil = zeros(sz(3),R,sz(1)*sz(4));
Bstil = zeros(R,R,sz(1)*sz(4));
BilDil = Pil;

while abs(ssr1 - ssr2) > eps && abs(ssr1 - ssr2)/abs(ssr2) > eps && iter < maxiter
    
    %update error
    ssr1 = ssr2;
    
    %PARAFAC2 - Xkl
    KL = sz(3)*sz(4);
    
    %Pkl Estimation
    for kl = 1:KL
        [U,~,V] = svds(Bkl(:,:,kl)*Bskl,R);
        Pkl(:,:,kl) = U(:,:)*V(:,:)';
    end
    
    %Bskl estimation
    for kl = 1:KL
        Bstkl(:,:,kl) = mkl(kl).*Pkl(:,:,kl)'*Bkl(:,:,kl);
    end
    
    Bskl = 1/(sum(mkl))*(sum(Bstkl,3));
    
    for rr = 1:R
        Bskl(:,rr) = Bskl(:,rr)./norm(Bskl(:,rr));
    end
    
    %Akl Estimation
    for kl = 1:sz(3)*sz(4)
        BklDkl(:,:,kl) = Bkl(:,:,kl)*Dkl(:,:,kl);
    end
    
    BklDklIK = permute(BklDkl,[1,3,2]);
    BklDklIK = reshape(BklDklIK,size(BklDklIK,1)*size(BklDklIK,2),size(BklDklIK,3));
    
    Xjkil = permute(Xkl,[1,3,2]);
    Xjkil = reshape(Xjkil,size(Xjkil,1)*size(Xjkil,2),size(Xjkil,3));
    
    if nnls == 1
        for jj = 1:sz(2)
            Akl(jj,:) = fnnls(BklDklIK'*BklDklIK + mA.*eye(R),(mA.*Ail(jj,:) + Xjkil(:,jj)'*BklDklIK)');
        end
    elseif nnls == 2
        Akl = fcnnls(BklDklIK'*BklDklIK + mA.*eye(R),(mA.*Ail + Xjkil'*BklDklIK)')';
    elseif nnls == 3
        Akl = (pinv(BklDklIK'*BklDklIK + mA.*eye(R))*(mA.*Ail + Xjkil'*BklDklIK)')';
    end
    
    for rr = 1:R
        Akl(:,rr) = Akl(:,rr)./norm(Akl(:,rr));
    end
    
    Akl(isnan(Akl)) = 0;
    
    %Bkl Estimation
    for kl = 1:KL
        Bkl(:,:,kl) = (pinv(Dkl(:,:,kl)*(Akl'*Akl)*Dkl(:,:,kl) + mkl(kl)*eye(R))*(Xkl(:,:,kl)*Akl*Dkl(:,:,kl) + mkl(kl)*Pkl(:,:,kl)*Bskl)')';
    end
    
    for kl = 1:KL
        for rr = 1:R
            Btkl = Bkl(:,rr,kl);
            Bkl(:,rr,kl) = Btkl;
            Bkl(:,rr,kl) = Bkl(:,rr,kl)./norm(Bkl(:,rr,kl));
        end
    end
    
    Bkl(isnan(Bkl)) = 0;
    
    %Dkl Estimation (not non-negative)
    for kl = 1:KL
        Dkl(:,:,kl) = diag(diag(pinv(Akl'*Akl)*(Akl'*Xkl(:,:,kl)'*Bkl(:,:,kl))*pinv(Bkl(:,:,kl)'*Bkl(:,:,kl))));
    end
    
    %mkl update
    if iter == 1
        for kl = 1:sz(3)*sz(4)
            S = svds(Xkl(:,:,kl)-mean(Xkl(:,:,kl)),2);
            SNRkl = S(1)/S(2);
            mkl(kl) = 10^(-SNRkl/10)*sum(sum(sum(sqrt((Xkl(:,:,kl) - Bkl(:,:,kl)*Dkl(:,:,kl)*Akl').^2))))/sum(sum(sum(sqrt((Bkl(:,:,kl) - Pkl(:,:,kl)*Bskl).^2))));
        end
    else
        for kl = 1:sz(3)*sz(4)
            if iter < coupit
                mkl(kl) = min(mkl(kl)*1.02,1e12);
            end
        end
    end
    
    
    %%PARAFAC2 - Xil
    IL = sz(1)*sz(4);
    
    %Pil Estimation
    for il = 1:IL
        [U,~,V] = svds(Bil(:,:,il)*Bsil,R);
        Pil(:,:,il) = U(:,:)*V(:,:)';
    end
    
    %Bsil estimation
    for il = 1:IL
        Bstil(:,:,il) = mil(il).*Pil(:,:,il)'*Bil(:,:,il);
    end
    
    Bsil = 1/(sum(mil))*(sum(Bstil,3));
    
    for rr = 1:R
        Bsil(:,rr) = Bsil(:,rr)./norm(Bsil(:,rr));
    end
    
    %Ail estimation
    for il = 1:sz(1)*sz(4)
        BilDil(:,:,il) = Bil(:,:,il)*Dil(:,:,il);
    end
    
    BilDilIK = permute(BilDil,[1,3,2]);
    BilDilIK = reshape(BilDilIK,size(BilDilIK,1)*size(BilDilIK,2),size(BilDilIK,3));
    
    Xjilk = permute(Xil,[1,3,2]);
    Xjilk = reshape(Xjilk,size(Xjilk,1)*size(Xjilk,2),size(Xjilk,3));
    
    if nnls == 1
        for jj = 1:sz(2)
            Ail(jj,:) = fnnls(BilDilIK'*BilDilIK + mA.*eye(R),(mA.*Akl(jj,:) + Xjilk(:,jj)'*BilDilIK)');
        end
    elseif nnls == 2
        Ail = fcnnls(BilDilIK'*BilDilIK + mA.*eye(R),(mA.*Akl + Xjilk'*BilDilIK)')';
    elseif nnls == 3
        Ail = (pinv(BilDilIK'*BilDilIK + mA.*eye(R))*(mA.*Akl + Xjilk'*BilDilIK)')';
    end
    
    for rr = 1:R
        Ail(:,rr) = Ail(:,rr)./norm(Ail(:,rr));
    end
    
    Ail(isnan(Ail)) = 0;
    
    for il = 1:IL
        Bil(:,:,il) = (pinv(Dil(:,:,il)*(Ail'*Ail)*Dil(:,:,il) + mil(il)*eye(R))*(Xil(:,:,il)*Ail*Dil(:,:,il) + mil(il)*Pil(:,:,il)*Bsil)')';
    end
    
    for il = 1:IL
        for rr = 1:R
            Btil = Bil(:,rr,il);
            Bil(:,rr,il) = Btil;
            Bil(:,rr,il) = Bil(:,rr,il)./norm(Bil(:,rr,il));
        end
    end
    
    Bil(isnan(Bil)) = 0;
    
    %Dil Estimation
    for il = 1:sz(1)*sz(4)
        Dil(:,:,il) = diag(diag(pinv(Ail'*Ail)*(Ail'*Xil(:,:,il)'*Bil(:,:,il))*pinv(Bil(:,:,il)'*Bil(:,:,il))));
    end
    
    %mil Estimation
    if iter == 1
        for il = 1:sz(1)*sz(4)
            S = svds(Xil(:,:,il)-mean(Xil(:,:,il)),2);
            SNRil = S(1)/S(2);
            mil(il) = 10^(-SNRil/10)*sum(sum(sum(sqrt((Xil(:,:,il) - Bil(:,:,il)*Dil(:,:,il)*Ail').^2))))/sum(sum(sum(sqrt((Bil(:,:,il) - Pil(:,:,il)*Bsil).^2))));
        end
    else
        for il = 1:sz(1)*sz(4)
            if iter < coupit
                mil(il) = min(mil(il)*1.02,1e12);
            end
        end
    end
    
    
    for kl = 1:sz(3)*sz(4)
        Xhkl(kl) = norm(Xkl(:,:,kl) - Bkl(:,:,kl)*Dkl(:,:,kl)*Akl','fro')^2; %#ok
        mhkl(kl) = norm(Bkl(:,:,kl) - Pkl(:,:,kl)*Bskl,'fro')^2;
    end
    
   mhkl(isnan(mhkl)) = 0;
    
    for il = 1:sz(1)*sz(4)
        Xhil(il) = norm(Xil(:,:,il) - Bil(:,:,il)*Dil(:,:,il)*Ail','fro')^2; %#ok
        mhil(il) = norm(Bil(:,:,il) - Pil(:,:,il)*Bsil,'fro')^2;
    end
    
    mhil(isnan(mhil)) = 0;
    
    mhail = norm(Akl - Ail,'fro')^2;%/(norm(Akl,'fro'))^2;
    
    %Calculate the SSR as the norm of the residuals of the two models
    ssr2 = ((sum(Xhkl(:)) + sum(mhkl) + sum(Xhil(:)) + sum(mhil) + mhail)/(2*YNorm));
    ssr(iter) = ssr2; %#ok
    
    if iter == 1 && displ == 1
        varNames = {'Iteration','Absolute dError','Relative dError','SSR','mA'};
        fprintf(1,'\n%s\t\t%s\t\t%s\t\t%s\t\t%s\n',varNames{:})
        fprintf(1,' \t\t%d\t\t%e\t\t%e\t\t%e\t\t%e\n',[iter,abs(ssr2-ssr1),abs(ssr2-ssr1)/abs(ssr2),ssr(iter),mA]);
    elseif displ == 1
        fprintf(1,' \t\t%d\t\t%e\t\t%e\t\t%e\t\t%e\n',[iter,abs(ssr2-ssr1),abs(ssr2-ssr1)/abs(ssr2),ssr(iter),mA]);
    else
    end
    
    %Coupling constant determination
    if iter < coupit
        if iter == 1
            mA = (sum(Xhkl(:)) + sum(Xhil(:)))/norm(Akl,'fro')^2;
        end
        mA = min(1.03*mA,1e20);
    else
    end
    
    if animate == 1
        
        colmat = lines(R);
        
        subplot(2,1,1);
        
        for rr = 1:R
            if rr == 1
                plot(Akl(:,rr),'color',[colmat(rr,:),0.5],'LineWidth',1.5);
            else
                hold on;
                plot(Akl(:,rr),'color',[colmat(rr,:),0.5],'LineWidth',1.5);
            end
        end
        ylabel('A_{kl}','FontSize',14);
        title(sprintf('Relative Change in Error: %.4e', abs(ssr2-ssr1)/abs(ssr2)),'FontSize',14);
        hold off;
        
        subplot(2,1,2);
        for rr = 1:R
            if rr == 1
                plot(Ail(:,rr),'color',[colmat(rr,:),0.5],'LineWidth',1.5);
            else
                hold on;
                plot(Ail(:,rr),'color',[colmat(rr,:),0.5],'LineWidth',1.5);
            end
        end
        
        title(sprintf('Mass Spectral Coupling Constant: %.4e',mA),'FontSize',14);
        ylabel('A_{il}','FontSize',14);
        xlabel('m/z (Simulated)','FontSize',14);
        hold off;
        
        set(gcf,'color','w');
        
        drawnow;
        
        frame = getframe(h);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        
        if iter == 1
            imwrite(imind,cm,filename,'gif','Loopcount',3);
        else
            imwrite(imind,cm,filename,'gif','WriteMode','append');
        end
        
    end
    
    iter = iter + 1;
    
end

%Average the mass spectra, scores.
A1 = Ail + Akl;

for rr = 1:R
    A1(:,rr) = A1(:,rr)./norm(A1(:,rr));
end

A = A1;

%Create Scores Matrix, F

for kl = 1:sz(3)*sz(4)
    BklDkl(:,:,kl) = Bkl(:,:,kl)*Dkl(:,:,kl);
end

for il = 1:sz(1)*sz(4)
    BilDil(:,:,il) = Bil(:,:,il)*Dil(:,:,il);
end

Bklvec = 0:sz(3):sz(3)*sz(4);
Bilvec = 0:sz(1):sz(1)*sz(4);

for ll = 1:sz(4)
    %Extract one sample at a time - not efficient but whatever
    tr2p = BklDkl(:,:,(Bklvec(ll)+1):(Bklvec(ll+1)));
    tr1p = BilDil(:,:,(Bilvec(ll)+1):(Bilvec(ll+1)));
    
    tr1p = permute(tr1p,[1,3,2]);
    tr2p = permute(tr2p,[3,1,2]);
    
    F11(:,:,:,ll) = tr1p + tr2p; %#ok
    
end

F = permute(F11,[2,3,1,4]);

%Normalise scores matrices (averaged)
for ll = 1:sz(4)
    for rr = 1:R
        F(:,rr,:,ll) = F(:,rr,:,ll)./norm(squeeze(F(:,rr,:,ll)),'fro');
    end
end

D = zeros(R,sz(4));
for ll = 1:sz(4)
    %unfold the tensor
    X = permute(Xijkl(:,:,:,ll),[1 3 2]);
    X = reshape(X,[sz(1)*sz(3),sz(2)]);
    
    Ft = permute(F(:,:,:,ll),[1,3,2]);
    Ft = reshape(Ft,[sz(1)*sz(3),R]);
    
    D(:,ll) = diag(pinv(Ft'*Ft)*Ft'*X*pinv(A'));
end

%Calculate the percent variance explained in the final model
Xh = zeros(sz);
for kk = 1:sz(3)
    for ll = 1:sz(4)
        Xh(:,:,kk,ll) = F(:,:,kk,ll)*diag(D(:,ll))*A';
    end
end

pvar = 100*(1 - sum(sum(sum(sum((Xijkl - Xh).^2))))/sum(sum(sum(sum(Xijkl.^2)))));
if displ == 1
    disp('At least one of the conditions for convergence has been reached')
else
end

end
