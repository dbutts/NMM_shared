function nim_out = NMMfit_filters( nim, Robs, Xstims, Gmults, Uindx, silent, desired_optim_params, regmat_custom, targets, fit_thresh )
%
% Usage: nim_out = NMMfit_filters( nim, Robs, Xstims, <Gmults>, <Uindx>, <silent>, <desired_optim_params>, <regmat_custom>,<targets>,  <fit_thresh> )
%              or
% Usage: nim_out = NMMfit_filters( nim, Robs, Xstims, <Gmults>, <Uindx>, <fit_params> )
%   this is an alternate usage with struct 'fit_params'. Note fit_params must have field '.silent' set at minimum
%
% Optimizes the stimulus filters and threshold of threshlin terms
%
% INPUTS:
%       nim: model structure
%       Robs: binned spikes
%       Xstim: time-embedded stimulus mat
%       <Gmults>: multiplicative term to be combined with module output
%       <silent>: 0 to display optimization iterations, and 1 to suppress them
%       <desired_optim_params>: Struct of optimization parameters
%       <targets>: Vector of indices specifying which subunits to optimize.
%           (-1 = spk history filter, and -2 = offset only) Default is to optimize all elements
%       <fit_thresh>: 1 to fit internal thresholds, otherwise set
%	                    to 0 (default). Can also pass in array of length Nmods
%
%       <fit_params>: structure that contains common fit parameters:
%             .silent = 0,1 (see above)
%             .fit_thresh = 0,1 (up to # of mods -- see above
%             .opt_params (struct of optimization parameters)
%             .rescaleNLs (relevant if fitting nonparametric upstream nonlinearities (see NMMfit_upstreamNLs)
%             .regmat_custom
%             .targets

% OUTPUTS:
%       nim_out: output model struct

%% PROCESS INPUTS
Nmods = length(nim.mods);
if (nargin < 4) || (length(Gmults) < Nmods)
	Gmults{Nmods} = [];
end

if nargin < 5
	Uindx = [];
end

% Process Xstims (in case multiple Xstims)
if ~iscell(Xstims)
	tmp = Xstims;
	clear Xstims
	Xstims{1} = tmp;
end
if (nargin < 6) || isempty(silent)
	silent = 0;
end
if nargin < 7
	desired_optim_params = [];
end
if nargin < 8
	regmat_custom = [];
end
if nargin < 9
	targets = [];
end
if (nargin < 10) || isempty(fit_thresh)
	fit_thresh = 0;
end

if isfield(silent,'silent')
	fit_params = silent;
	clear silent
	silent = fit_params.silent;
	if isfield(fit_params,'optim_params')
		desired_optim_params = fit_params.optim_params;
	end
	if isfield(fit_params,'fit_thresh')
		fit_thresh = fit_params.fit_thresh;
	end
	if isfield(fit_params,'regmat_custom')
		regmat_custom = fit_params.regmat_custom;
	end
	if isfield(fit_params,'targets')
		targets = fit_params.targets;
	end
end

% Index X-matrices and Robs
RobsFULL = Robs;
XstimsFULL = Xstims;
if ~isempty(Uindx)
  for nn = 1:length(Xstims)
    Xstims{nn} = Xstims{nn}(Uindx,:);
  end
  Robs = RobsFULL(Uindx);
end

% Make sure Robs is a column vector
Robs = Robs(:);

if length(fit_thresh) ~= Nmods
	fit_thresh = ones(Nmods,1)*fit_thresh;  % expand thresh-fit to all modules
end
for nn = 1:Nmods
	[NT,filtLen] = size(Xstims{nim.mods(nn).Xtarget}); % stimulus dimensions
	if filtLen ~= prod(nim.stim_params(nim.mods(nn).Xtarget).stim_dims)
		error('Xstim dimensions dont match with stim_params')
	end
	% Look for modules where threshold-fits necessary
	if fit_thresh(nn) > 0
		if strcmp(nim.mods(nn).NLtype,'lin')
			fit_thresh(nn) = 0;
		elseif strcmp(nim.mods(nn).NLtype,'nonpar')
			fit_thresh(nn) = 0;
		end
	end			
end

spkhstlen = nim.spk_hist.spkhstlen;

if max(targets) > Nmods % check input targets
	error('Invalid target specified');
end
if isempty(targets) %default is to optimize all model components
	targets = 1:Nmods;
	if (spkhstlen > 0)
		targets = [targets -1]; %optimize spike hist filter
	end
elseif targets == -2
	targets = [];
end

Ntargets = sum(targets > 0); % number of targeted subunits
non_targets = setdiff([1:Nmods -1 -2],targets); % elements of the model held constant

if ismember(-1,targets) && spkhstlen == 0
	error('No spk history term initialized')
end

% Add spk NL constant if it isnt already there
if length(nim.spk_NL_params) < 4
	nim.spk_NL_params(4) = 0;
end

if ~silent
	fprintf( 'Filt-thresh optimization. Targets:  ' )
	disp(sprintf('%d ', targets ))
end

%% PARSE INITIAL PARAMETERS
% Compute initial fit parameters
initial_params = [];
sign_con = [];
for imod = targets(targets > 0)
	cur_kern = nim.mods(imod).filtK';
	if nim.mods(imod).Kcon ~= 0
    sign_con(length(initial_params)+(1:length(cur_kern))) = nim.mods(imod).Kcon;
	end
	
	% Add coefs to initial param vector, along with threshold
	if fit_thresh(imod)
		initial_params = [initial_params; cur_kern'; nim.mods(imod).NLx;];
	else
		initial_params = [initial_params; cur_kern';];
	end
end

% Add in spike history coefs
if ismember(-1,targets)
	initial_params = [initial_params; nim.spk_hist.coefs];
end

initial_params(end+1) = nim.spk_NL_params(1); % add constant offset
initial_params = initial_params(:);

%% COMPUTE L1 PENALTY IF APPLICABLE
lambda_L1 = zeros(size(initial_params));
cnt = 0;
for ii = 1:Ntargets
	filtLen = length(nim.mods(targets(ii)).filtK);
	cur_inds = cnt + (1:filtLen);
	lambda_L1(cur_inds) = nim.mods(targets(ii)).reg_params.lambda_L1;
	cnt = cnt + filtLen+fit_thresh(targets(ii));
end
lambda_L1 = lambda_L1/sum(Robs); % since we are dealing with LL/spk

%% PRECOMPUTE 'TENT-BASIS' DERIVATIVES OF UPSTREAM NLS IF NEEDED
if any(strcmp('nonpar',{nim.mods(targets(targets > 0)).NLtype}))
	for ii = 1:Ntargets
		if strcmp(nim.mods(targets(ii)).NLtype,'nonpar')
			NLx = nim.mods(targets(ii)).NLx;
			NL = nim.mods(targets(ii)).NLy;
            
			% Compute derivative of non-linearity
			fpr = zeros(1,length(NLx)-1);
			for n = 1:length(fpr)
				fpr(n) = (NL(n+1)-NL(n))/(NLx(n+1)-NLx(n));
			end
			fprimes{ii} = fpr;
		else
			fprimes{ii} = [];
		end
	end
else
	fprimes = [];
end

%% CREATE SPIKE HISTORY Xmat IF NEEDED
if nim.spk_hist.spkhstlen > 0
	Xspkhst = create_spkhist_Xmat( RobsFULL, nim.spk_hist.bin_edges );
  if ~isempty(Uindx)
    Xspkhst = Xspkhst(Uindx,:);
  end
else
	Xspkhst = [];
end

%% COMPUTE NET OUPTUT OF ALL NON-TARGET PREDICTORS
nt_gout = zeros(NT,1);
%Kmat = [nim.mods(:).filtK]; %filter matrix
%gint = Xstim*Kmat; %filter output of each subunit

for imod = non_targets(non_targets > 0) %for all subunits that aren't targeted
    
	fgint = Xstims{nim.mods(imod).Xtarget} * nim.mods(imod).filtK;
    
	% Process subunit g's with upstream NLs
	if strcmp(nim.mods(imod).NLtype,'nonpar')
		fgint = piecelin_process( fgint, nim.mods(imod).NLy, nim.mods(imod).NLx );
	elseif strcmp(nim.mods(imod).NLtype,'quad')
		%fgint = gint(:,imod).^2;
		fgint = (fgint-nim.mods(imod).NLx).^2;
	elseif strcmp(nim.mods(imod).NLtype,'lin')
		%fgint = gint(:,imod);
	elseif strcmp(nim.mods(imod).NLtype,'threshlin')
		fgint = fgint-nim.mods(imod).NLx;
		fgint(fgint < 0) = 0;
	elseif strcmp(nim.mods(imod).NLtype,'threshP')
		fgint = fgint-nim.mods(imod).NLx;
		fgint(fgint < 0) = 0;
		fgint = fgint.^(nim.mods(imod.NLy));
	else
		error('Invalid internal NL');
	end
    
	% Multiply by weight (and multiplier, if appl) and add to generating function
	if isempty(Gmults{imod})
		nt_gout = nt_gout + fgint*nim.mods(imod).sign;
	else
		nt_gout = nt_gout + (fgint.*Gmults{imod}) * nim.mods(imod).sign;
	end
end

if ismember(-1,non_targets) && spkhstlen > 0
	nt_gout = nt_gout + Xspkhst*nim.spk_hist.coefs(:);
end

%% IDENTIFY ANY CONSTRAINTS
use_con = 0;
A = []; Aeq = []; %initialize constraint matrices
LB = -Inf*ones(size(initial_params));
UB = Inf*ones(size(initial_params));

% Constrain any of the filters to be positive or negative
if ~isempty(sign_con)
  LB(sign_con == 1) = 0;
  UB(sign_con == -1) = 0;
 
  use_con = 1;
end

if spkhstlen > 0 && ismember(-1,targets) %if optimizing spk history term
	%negative constraint on spk history coefs
	if nim.spk_hist.negCon == 1
		%spkhist_inds = (Ntargets*filtLen + 1):(Ntargets*filtLen + spkhstlen);
		spkhist_inds = cnt + (1:spkhstlen);
		UB(spkhist_inds) = 0;
        
		use_con = 1;
	end
end

beq = zeros(size(Aeq,1),1);
b = zeros(size(A,1),1);

%% GENERATE REGULARIZATION MATRICES
L2_mats = create_L2_matrices_NMM( nim );
L2_mats.custom = regmat_custom;

%%
if max(lambda_L1) > 0 && use_con == 1
	disp('Cant use L1 with constrained optimization, aborting constraints');
	use_con = 0;
end

%%
optim_params.MaxFunEvals = 100*length(initial_params);
optim_params.MaxIter = 1e3;
optim_params.Display = 'off';
if silent == 0
    optim_params.Display = 'iter';
end
if use_con == 0 %if no constraints
    
    %if using L1 reg
    if max(lambda_L1) > 0
        if exist('L1General2_PSSas','file') == 2
            optim_params.optTol = 1e-4;
            optim_params.progTol = 1e-8;
            if silent == 0
                optim_params.verbose = 2;
            else
                optim_params.verbose = 0;
            end
            % Load in specified optimization parameters
            if ~isempty(desired_optim_params)
                spec_fields = fieldnames(desired_optim_params);
                for i = 1:length(spec_fields)
                    optim_params = setfield(optim_params,spec_fields{i},getfield(desired_optim_params,spec_fields{i}));
                end
            end
            
            [params] = L1General2_PSSas(@(K) LLfit_filters_internal(nim, K, Robs, Xstims,Xspkhst,Gmults,L2_mats,targets,nt_gout,fprimes,fit_thresh),...
                initial_params,lambda_L1,optim_params);
            %             [params] = L1General2_PSSas(@(K) NIM_fit_filters_internal(nim, K, Robs, Xstim,Xspkhst,XLin,L2_mats,targets,nt_gout,fprimes),...
            %                 initial_params,lambda_L1);
        else
            error('Need to install Mark Schmidts L1 optimization toolbox for using L1');
        end
    else % if not using L1 reg
        
        if exist('minFunc','file') == 2 %try to use Mark Schmidt's optimizer
            % if using Mark Schmidt's optimization, some differences in option parameters
            optim_params.optTol = 1e-4;
            optim_params.progTol = 1e-8;
            optim_params.Method = 'lbfgs';
            if silent == 0
                optim_params.verbose = 2;
            else
                optim_params.verbose = 0;
            end
            %load in specified optimization parameters
            if ~isempty(desired_optim_params)
                spec_fields = fieldnames(desired_optim_params);
                for i = 1:length(spec_fields)
                    optim_params = setfield(optim_params,spec_fields{i},getfield(desired_optim_params,spec_fields{i}));
                end
            end
            
            [params] = minFunc( @(K) LLfit_filters_internal( nim, K, Robs, Xstims, Xspkhst, Gmults, L2_mats, targets, nt_gout, fprimes, fit_thresh ), initial_params, optim_params);
            
        else %if using Matlab Optim toolbox:
            
            % Default optimization parameters
            optim_params.LargeScale = 'off';
            optim_params.TolFun = 1e-6;
            optim_params.TolX = 1e-6;
            optim_params.HessUpdate = 'bfgs';
            optim_params.GradObj = 'on';
            
            %load in specified optimization parameters
            if ~isempty(desired_optim_params)
                spec_fields = fieldnames(desired_optim_params);
                for i = 1:length(spec_fields)
                    optim_params = setfield(optim_params,spec_fields{i},getfield(desired_optim_params,spec_fields{i}));
                end
            end
            
            [params] = fminunc( @(K) LLfit_filters_internal( nim, K, Robs, Xstims, Xspkhst, Gmults, L2_mats, targets, nt_gout, fprimes, fit_thresh ), initial_params, optim_params);
            
        end
    end
else %if there are constraints
    
    % Try to use Mark Schmidt constrained optimizer
    if exist('minConf_TMP','file') == 2 && isempty(A) && isempty(Aeq)
        % if using Mark Schmidt's optimization, some differences in option parameters
        optim_params.optTol = 1e-4;
        optim_params.progTol = 1e-6;
        if silent == 0
            optim_params.verbose = 2;
        else
            optim_params.verbose = 0;
        end
        [params] = minConf_TMP( @(K) LLfit_filters_internal( nim, K, Robs, Xstims, Xspkhst, Gmults, L2_mats, targets, nt_gout, fprimes, fit_thresh ),...
            initial_params, LB,UB,optim_params);
    else
        % otherwise resort to matlab's
        optim_params.GradObj = 'on';
        optim_params.LargeScale = 'off';
        optim_params.Algorithm = 'active-set';
        optim_params.optTol = 1e-4;
        optim_params.progTol = 1e-6;
        [params] = fmincon( @(K) LLfit_filters_internal( nim, K, Robs, Xstims, Xspkhst, Gmults, L2_mats, targets, nt_gout, fprimes, fit_thresh ),...
            initial_params, A,b,Aeq,beq,LB,UB,[],optim_params);
    end
    
end

%% PARSE MODEL FIT
nim_out = nim;
nim_out.spk_NL_params(1) = params(end);

cnt = 0;
for ii = 1:Ntargets
	filtLen = length(nim.mods(targets(ii)).filtK);
	cur_kern = params( cnt + (1:filtLen) );
	nim_out.mods(targets(ii)).filtK = cur_kern(:);
	cnt = cnt + filtLen;
	if fit_thresh(ii)
		nim_out.mods(targets(ii)).NLx = params( cnt + 1 );  % pull new threshold
		cnt = cnt + 1;
	end
end
if ismember(-1,targets)
  nim_out.spk_hist.coefs = params(cnt + (1:spkhstlen));
end

[LL, nullLL, ~, G, gint, fgint, penLL] = NMMeval_model( nim_out, RobsFULL, XstimsFULL, Gmults, Uindx, regmat_custom );
nim_out.LL_seq = cat(1,nim_out.LL_seq,LL);
nim_out.penLL_seq = cat(1,nim_out.penLL_seq,penLL);
nim_out.opt_history = cat(1,nim_out.opt_history,{'filt'});

% Compute std dev of the output of each subunit
for n = 1:Nmods
	mod_norm = std(fgint(:,n));
	nim_out.mods(n).mod_norm = mod_norm;
end

end % End of main body of function

%%%% INTERNAL FUNCTIONS %%%%%%%

function [LL, LLgrad] = LLfit_filters_internal(nim, params, Robs, Xstims, Xspkhst, Gmults, L2_mats, targets, nt_gout, fprimes, fit_thresh )
%
% [LL, LLgrad] = LLfit_filters_internal(nim, params, Robs, Xstims, Xspkhst, Gmults, L2_mats, targets, nt_gout, fprimes, fit_thresh )
%
% Internal function for computing LL and LLgradient with respect to the
% stimulus filters

%% USEFUL PARAMETERS
Ntargets = sum(targets > 0);
%lin_dims = nim.stim_params.lin_dims;
spkhstlen = nim.spk_hist.spkhstlen;

%% ESTIMATE GENERATING FUNCTIONS (OVERALL AND INTERNAL)
theta = params(end); % offset
G = theta + nt_gout; % initialize overall generating function G

%kmat = reshape(params(1:Ntargets*filtLen),filtLen,Ntargets);
%gint = Xstim*kmat;

gint = nan(length(Robs),Ntargets);

NKtot = 0;  filtLen = zeros(Ntargets,1);  ks = cell(Ntargets,1);
threshs = zeros(Ntargets,1);
for ii = 1:Ntargets
    
	tar = targets(ii);
    
	% Pull out (potentially different-sized) filters from params
	filtLen(ii) = prod(nim.stim_params(nim.mods(tar).Xtarget).stim_dims);
	ks{ii} = params( NKtot + (1:filtLen(ii)) );
	NKtot = NKtot + filtLen(ii);
	if fit_thresh(tar)
		threshs(ii) = params( NKtot + 1 );
		NKtot = NKtot + 1;
	else
		if ~isempty(nim.mods(tar).NLx)
			threshs(ii) = nim.mods(tar).NLx;
		end
	end  
	gint(:,ii) = Xstims{nim.mods(tar).Xtarget} * ks{ii};
    
	% Process subunit g's with upstream NLs
	if strcmp(nim.mods(tar).NLtype,'nonpar')
		fgint = piecelin_process(gint(:,ii),nim.mods(tar).NLy,nim.mods(tar).NLx);
	elseif strcmp(nim.mods(tar).NLtype,'quad')
		fgint = (gint(:,ii)-threshs(ii)).^2;
	elseif strcmp(nim.mods(tar).NLtype,'lin')
		fgint = gint(:,ii);
	elseif strcmp(nim.mods(tar).NLtype,'threshlin')
		fgint = gint(:,ii)-threshs(ii);
		fgint(fgint < 0) = 0;
	elseif strcmp(nim.mods(tar).NLtype,'threshP')
		fgint = gint(:,ii)-threshs(ii);
		fgint(fgint < 0) = 0;
		fgint = fgint.^nim.mods(tar).NLy;
	else
		error('Invalid internal NL');
	end
    
	% Multiply by weight (and multiplier, if appl) and add to generating function
	if isempty(Gmults{tar})
		G = G + fgint*nim.mods(tar).sign;
	else
		G = G + (fgint.*Gmults{tar}) * nim.mods(tar).sign;
	end
    
end

% Add contribution from spike history filter
if spkhstlen > 0 && ismember(-1,targets)
	G = G + Xspkhst*params(NKtot + (1:spkhstlen));
end

%% Compute predicted firing rate
if strcmp(nim.spk_NL_type,'logexp')
    max_gbeta = 50; %to prevent numerical overflow
    bgint = G*nim.spk_NL_params(2); %g*beta
    expg = exp(bgint);
    too_large = (bgint > max_gbeta);
    r = nim.spk_NL_params(4) + nim.spk_NL_params(3)*log(1+expg); %alpha*log(1+exp(gbeta))
    r(too_large) = nim.spk_NL_params(4) + nim.spk_NL_params(3)*bgint(too_large); %log(1+exp(x)) ~ x in limit of large x
elseif strcmp(nim.spk_NL_type,'logistic')
    bgint = G*nim.spk_NL_params(2); %g*beta
    expg = exp(-bgint);
		r = nim.spk_NL_params(4) + 1./(1+expg); %1/(1+exp(-gbeta))  % note third param not used
elseif strcmp(nim.spk_NL_type,'exp')
    expg = exp(G);
    r = expg;
elseif strcmp(nim.spk_NL_type,'linear')
    r = G;    
else
    error('invalid spk nl');
end

% Enforce minimum predicted firing rate to avoid nan LLs
if ~strcmp(nim.spk_NL_type,'linear')
	min_pred_rate = 1e-50;
	if min(r) < min_pred_rate
		r(r < min_pred_rate) = min_pred_rate; %minimum predicted rate
	end
end

%% COMPUTE LL and LL gradient
if strcmp(nim.spk_NL_type,'linear') % use MSE as cost function 
  Nspks = length(Robs);
  LL = -sum( (Robs - r).^2 );
elseif strcmp(nim.spk_NL_type,'logistic')
	% Bernouli likelihood = robs log r + (1-robs) log (1-r)
	Nspks = sum(Robs);
	%spklocs = find(Robs == 1);
	%zerlocs = find(Robs == 0);
	%LL = sum(log(r(spklocs))) + sum(log(1-r(zerlocs)));
	LL = nansum( Robs.*log(r) + (1-Robs).*log(1-r) );
else  % Poisson likelihood
	Nspks = sum(Robs);
	LL = sum(Robs.* log(r) - r); %up to an overall constant
	%'residual' = (R/r - 1)*F'[] where F[.] is the spk NL
end

if strcmp(nim.spk_NL_type,'logexp')
	% 'residual' = (R/r - 1)*F'[] where F[.] is the spk NL
  residual = nim.spk_NL_params(3)*nim.spk_NL_params(2)*(Robs./r - 1) .* expg ./ (1+expg);
  residual(too_large) = nim.spk_NL_params(3)*nim.spk_NL_params(2)*(Robs(too_large)./r(too_large) - 1);
  residual(r == min_pred_rate) = 0; %for points beneath the lower threshold of the spk NL, take F'[.] = 0
elseif strcmp(nim.spk_NL_type,'logistic')
	% 'residual' = (R/r - (1-R)/(1-r))*F'[] where F[.] is the spk NL
	residual = nim.spk_NL_params(2)* (Robs./r - (1-Robs)./(1-r)) .* expg ./ (1+expg).^2;
elseif strcmp(nim.spk_NL_type,'exp')
	residual = Robs - r;
	residual(r == min_pred_rate) = 0;
elseif strcmp(nim.spk_NL_type,'linear')
	residual = 2*(Robs - r);    
else
	error('Unsupported spiking NL')
end

%initialize LL gradient
LLgrad = zeros(length(params),1);

% Calculate derivatives with respect to constant term (theta)
LLgrad(end) = sum(residual);

% Calculate derivative with respect to spk history filter
if spkhstlen > 0 && ismember(-1,targets)
	LLgrad(NKtot+(1:spkhstlen)) = residual'*Xspkhst;
end

% NOW COMPUTE LL grad WRT STIMULUS FILTERS (and thresholds)
% Calculate output of derivative module
chunk_size = 1000; %maximum chunk size for processing high-dimensional stim filters
if max(filtLen) <= chunk_size
	use_chunking = 0;
else
	use_chunking = 1;
	NChunks = ceil(filtLen/chunk_size); %divide stim filters into this many chunks for piecewise processing
end

placeholder = 0;
for ii = 1:Ntargets
	tar = targets(ii);
	if strcmp(nim.mods(tar).NLtype,'lin')
		% Check for multiplicative interactions
		if isempty(Gmults{tar})
			LLgrad(placeholder+(1:filtLen(ii))) = residual'*Xstims{nim.mods(tar).Xtarget} * nim.mods(tar).sign;
		else
			LLgrad(placeholder+(1:filtLen(ii))) = (Gmults{tar}.*residual')* Xstims{nim.mods(tar).Xtarget} * nim.mods(tar).sign;
		end
				
	else
		if strcmp(nim.mods(tar).NLtype,'nonpar')      
			fpg = piececonst_process(gint(:,ii),fprimes{ii}, nim.mods(tar).NLx);
		elseif strcmp(nim.mods(tar).NLtype,'quad')
			fpg = 2*(gint(:,ii)-threshs(ii));
		elseif strcmp(nim.mods(tar).NLtype,'threshlin')
			fpg = gint(:,ii) >= threshs(ii); 
		elseif strcmp(nim.mods(tar).NLtype,'threshP')
			fpg = gint(:,ii) >= threshs(ii);
			fpg = nim.mods(tar).NLy * fpg.*gint(:,ii) .^ (nim.mods(tar).NLy-1);
			if nim.mods(tar).NLy < 1
				% then there will be some nan
				fpg(isnan(fpg)) = 0;
			end
		else
			error('Unsupported NL type')
		end
		
		% Add gradient for threshold parameters
		if fit_thresh(tar)
			LLgrad(placeholder+filtLen(ii)+1) = -sum((fpg.*residual)) * nim.mods(tar).sign;
		end
		
		% Add for multiplicative interactions
		if ~isempty(Gmults{tar})
			fpg = fpg .* Gmults{tar};
		end
		LLgrad(placeholder+(1:filtLen(ii))) = (fpg.*residual)' * Xstims{nim.mods(tar).Xtarget} * nim.mods(tar).sign;
    
	end
		
	placeholder = placeholder + filtLen(ii) + fit_thresh(tar);
end

%% COMPUTE L2 PENALTIES AND ASSOCIATED CONTRIBUTIONS TO THE LL GRADIENT
smooth_penalty = zeros(Ntargets,1);
deriv_penalty = zeros(Ntargets,1);
ridge_penalty = zeros(Ntargets,1);
custom_penalty = zeros(Ntargets,1);

LLgrad_pen = zeros(size(LLgrad));
placeholder = 0;
for ii = 1:Ntargets
    
    tar = targets(ii);
    
    % Temporal regularization
    if nim.mods(tar).reg_params.lambda_dT > 0
        deriv_penalty(ii) = deriv_penalty(ii) + nim.mods(tar).reg_params.lambda_dT*sum((L2_mats.L2_dT{nim.mods(tar).Xtarget} * ks{ii}).^2);
        cur_grad_pen = 2*nim.mods(tar).reg_params.lambda_dT*(L2_mats.L2_dT{nim.mods(tar).Xtarget}' * L2_mats.L2_dT{nim.mods(tar).Xtarget} * ks{ii});
        LLgrad_pen(placeholder+(1:filtLen(ii))) = LLgrad_pen(placeholder+(1:filtLen(ii))) + cur_grad_pen;
    end
    if nim.mods(tar).reg_params.lambda_d2T > 0
        smooth_penalty(ii) = smooth_penalty(ii) + nim.mods(tar).reg_params.lambda_d2T*sum((L2_mats.L2_d2T{nim.mods(tar).Xtarget} * ks{ii}).^2);
        cur_grad_pen = 2*nim.mods(tar).reg_params.lambda_d2T*(L2_mats.L2_d2T{nim.mods(tar).Xtarget}' * L2_mats.L2_d2T{nim.mods(tar).Xtarget} * ks{ii});
        LLgrad_pen(placeholder+(1:filtLen(ii))) = LLgrad_pen(placeholder+(1:filtLen(ii))) + cur_grad_pen;
    end
    
    % Spatial (and spatiotemporal) regularization
    if nim.stim_params(nim.mods(tar).Xtarget).stim_dims(2) > 1
        % 	if nPix(1) > 1 % for spatial stimuli
        
        if nim.mods(tar).reg_params.lambda_dX > 0
            deriv_penalty(ii) = deriv_penalty(ii) + nim.mods(tar).reg_params.lambda_dX*sum((L2_mats.L2_dX{nim.mods(tar).Xtarget} * ks{ii}).^2);
            cur_grad_pen = 2*nim.mods(tar).reg_params.lambda_dX*(L2_mats.L2_dX{nim.mods(tar).Xtarget}' * L2_mats.L2_dX{nim.mods(tar).Xtarget} * ks{ii});
            LLgrad_pen(placeholder+(1:filtLen(ii))) = LLgrad_pen(placeholder+(1:filtLen(ii))) + cur_grad_pen;
        end
        if nim.mods(tar).reg_params.lambda_d2X > 0
            smooth_penalty(ii) = smooth_penalty(ii) + nim.mods(tar).reg_params.lambda_d2X*sum((L2_mats.L2_d2X{nim.mods(tar).Xtarget} * ks{ii}).^2);
            cur_grad_pen = 2*nim.mods(tar).reg_params.lambda_d2X*(L2_mats.L2_d2X{nim.mods(tar).Xtarget}' * L2_mats.L2_d2X{nim.mods(tar).Xtarget} * ks{ii});
            LLgrad_pen(placeholder+(1:filtLen(ii))) = LLgrad_pen(placeholder+(1:filtLen(ii))) + cur_grad_pen;
        end
        if nim.mods(tar).reg_params.lambda_d2XT > 0
            smooth_penalty(ii) = smooth_penalty(ii) + nim.mods(tar).reg_params.lambda_d2XT*sum((L2_mats.L2_d2XT{nim.mods(tar).Xtarget} * ks{ii}).^2);
            cur_grad_pen = 2*nim.mods(tar).reg_params.lambda_d2XT*(L2_mats.L2_d2XT{nim.mods(tar).Xtarget}' * L2_mats.L2_d2XT{nim.mods(tar).Xtarget} * ks{ii});
            LLgrad_pen(placeholder+(1:filtLen(ii))) = LLgrad_pen(placeholder+(1:filtLen(ii))) + cur_grad_pen;
        end
    end
    
    %for custom regularization
    if nim.mods(tar).reg_params.lambda_custom > 0
        custom_penalty(ii) = custom_penalty(ii) + nim.mods(tar).reg_params.lambda_custom*sum((L2_mats.custom * ks{ii}).^2);
        cur_grad_pen = 2*nim.mods(tar).reg_params.lambda_custom*(L2_mats.custom' * L2_mats.custom * ks{ii});
        LLgrad_pen(placeholder+(1:filtLen(ii))) = LLgrad_pen(placeholder+(1:filtLen(ii))) + cur_grad_pen;
    end
    
    if nim.mods(tar).reg_params.lambda_L2 > 0
        ridge_penalty(ii) = nim.mods(tar).reg_params.lambda_L2*(ks{ii}' * ks{ii});
        LLgrad_pen(placeholder+(1:filtLen(ii))) = LLgrad_pen(placeholder+(1:filtLen(ii))) + 2*nim.mods(tar).reg_params.lambda_L2*ks{ii};
    end
    
    placeholder = placeholder + filtLen(ii)+fit_thresh(tar);
end

LL = LL - sum(smooth_penalty) - sum(ridge_penalty) - sum(deriv_penalty) - sum(custom_penalty);
LLgrad = LLgrad - LLgrad_pen;


%% CONVERT TO NEGATIVE LLS AND NORMALIZE BY NSPKS

LL = -LL/Nspks;
LLgrad = -LLgrad/Nspks;

end
