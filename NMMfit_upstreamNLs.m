function nim_out = NMMfit_upstreamNLs( nim, Robs, Xstims, Gmults, Uindx, rescale_NLs, silent, desired_optim_params, regmat_custom, desired_targets )
%
% Usage: nim_out = NMMfit_upstreamNLs( nim, Robs, Xstims, <Gmults>, <Uindx>, <rescale_NLs>, <silent>, <desired_optim_params>, <regmat_custom>, <desired_targets> )
%              or
% Usage: nim_out = NMMfit_upstreamNLs( nim, Robs, Xstims, <Gmults{}>, <Uindx>, <fit_params> )
%   this is an alternate usage with struct 'fit_params'. Note fit_params must have field '.silent' set at minimum
%   see NMMfit_filters for more documentation about its use
%
% Optimizes the upstream NLs (in terms of tent-basis functions) (plus extra linear terms if desired) for
% given stimulus filters
%
% INPUTS:
%       nim: model structure
%       Robs: binned spikes
%       Xstim: time-embedded stimulus mat
%       <XLin>: Matrix specifying additional linear predictors
%       <desired_targets>: Vector of indices specifying which subunits to optimize.
%           (-1 = spk history filter, and -2 = extra linear filter) Default is to optimize all elements
%       <rescale_NLs>: set to 0 if you don't want to rescale upstream NLs after
%           estimating (otherwise set to 1)
%       <silent>: set to 0 if you want to turn on the optimization display
%       <desired_optim_params>: Struct of optimization parameters
% OUTPUTS:
%       nim_out: output model struct

%% USEFUL PARAMS
Nmods = length(nim.mods);
spkhstlen = nim.spk_hist.spkhstlen;

% Process Xstims (in case multiple Xstims)
if ~iscell(Xstims)
	tmp = Xstims;
	clear Xstims
	Xstims{1} = tmp;
end

for n = 1:Nmods
	[NT,filtLen] = size(Xstims{nim.mods(n).Xtarget}); % stimulus dimensions
	if filtLen ~= prod(nim.stim_params(nim.mods(n).Xtarget).stim_dims)
		error('Xstim dimensions dont match with stim_params')
	end
end

%% PARSE INPUTS
if (nargin < 4) || (length(Gmults) < Nmods)
  Gmults{Nmods} = [];
end

% Process Xstims (in case multiple Xstims)
if ~iscell(Xstims)
	tmp = Xstims;
	clear Xstims
	Xstims{1} = tmp;
end
if nargin < 5
	Uindx = [];
end
if (nargin < 6) || isempty(rescale_NLs)
	rescale_NLs = 1;
end
if (nargin < 7) || isempty(silent)
	silent = 0;
end
if (nargin < 8) || isempty(desired_optim_params)
	desired_optim_params = [];
end
if nargin < 9
	regmat_custom = [];
end
if nargin < 10
	desired_targets = [];
end

% Add any fit_params (if this was used in place of rescale_NLs
if isfield(rescale_NLs,'silent')
	fit_params = rescale_NLs;
	clear rescale_NLs
	if isfield(fit_params,'rescaleNLs')
		rescale_NLs = fit_params.rescaleNLs;
	elseif isfield(fit_params,'rescale_NLs')  % making sure not entered other way
		rescale_NLs = fit_params.rescale_NLs;
	else
		rescale_NLs = 1;
	end
	silent = fit_params.silent;

	if isfield(fit_params,'optim_params')
		desired_optim_params = fit_params.optim_params;
	end
	if isfield(fit_params,'regmat_custom')
		regmat_custom = fit_params.regmat_custom;
	end
	if isfield(fit_params,'targets')
		desired_targets = fit_params.targets;
	end
end

% Index X-matrices and Robs
RobsFULL = Robs;
if ~isempty(Uindx)
  for nn = 1:length(Xstims)
    Xstims{nn} = Xstims{nn}(Uindx,:);
  end
  Robs = RobsFULL(Uindx);
  NT = length(Robs);
end

% Make sure Robs is a column vector
Robs = Robs(:);

% Add spk NL constant if it isnt already there
if length(nim.spk_NL_params) < 4
	nim.spk_NL_params(4) = 0;
end

%% FIND SUBUNITS WITH NONPAR UPSTREAM NLS
poss_targets = [];
n_tbfs = [];
lambda_smooth = zeros(Nmods,1);
for imod = 1:Nmods
	if strcmp(nim.mods(imod).NLtype,'nonpar')
		lambda_smooth(imod) = nim.mods(imod).reg_params.lambda_NLd2;
		poss_targets = [poss_targets imod];
		n_tbfs = [n_tbfs; length(nim.mods(imod).NLx)];
	end
end
if isempty(poss_targets)
	disp('No subunits have nonpar NLs. Exiting function...');
	return
end
poss_targets = [poss_targets -1]; %add spk hist and lin filter to possible targets

if length(unique(n_tbfs)) ~= 1
	error('Have to have same number of tent-bases for each subunit');
end
n_tbfs = unique(n_tbfs);

% If targets are specified
if ~isempty(desired_targets)
	if any(~ismember(desired_targets,poss_targets))
		error('At least one desired target doesnt have non-parametric NL')
	end
	targets = desired_targets;
else %if no targets are specified use all possible targets (all subunits with nonpar NLs)
	targets = poss_targets;
end

non_targets = setdiff([1:Nmods -1],targets); %set of subunits that are not being optimized
Ntargets = sum(targets > 0); %number of target subunits

%% CREATE SPIKE HISTORY Xmat IF NEEDED
if spkhstlen > 0
	Xspkhst = create_spkhist_Xmat( RobsFULL, nim.spk_hist.bin_edges );
  if ~isempty(Uindx)
    Xspkhst = Xspkhst(Uindx,:);
  end
else
	Xspkhst = [];
end

%% COMPUTE NET OUPTUT OF ALL NON-TARGET SUBUNITS
nt_gout = zeros(NT,1);
%gint = Xstim*Kmat;
gint = nan(NT,Nmods);

for imod = 1:Nmods
	gint(:,imod) = Xstims{nim.mods(imod).Xtarget} * nim.mods(imod).filtK;
end

for imod = non_targets(non_targets > 0)

	% Process subunit g's with upstream NLs
	if strcmp(nim.mods(imod).NLtype,'nonpar')
		fgint = piecelin_process(gint(:,imod),nim.mods(imod).NLy,nim.mods(imod).NLx);
	elseif strcmp(nim.mods(imod).NLtype,'quad')
		fgint = (gint(:,imod)-nim.mods(imod).NLx).^2;
	elseif strcmp(nim.mods(imod).NLtype,'lin')
		fgint = gint(:,imod);
	elseif strcmp(nim.mods(imod).NLtype,'threshP')
		fgint = gint(:,imod)-nim.mods(imod).NLx;
		fgint(fgint < 0) = 0;
		fgint = fgint.^(nim.mods(imod.NLy));
	elseif strcmp(nim.mods(imod).NLtype,'threshlin')
		fgint = gint(:,imod)-nim.mods(imod).NLx;
		fgint(fgint < 0) = 0;
	else
		error('Invalid internal NL');
	end
    
	% Multiply by weight (and multiplier, if appl) and add to generating function
	if ~isempty(Gmults{imod})
		nt_gout = nt_gout + (fgint.*Gmults{imod}) * nim.mods(imod).sign;
	else
		nt_gout = nt_gout + fgint*nim.mods(imod).sign;
	end
end
if ismember(-1,non_targets) && spkhstlen > 0
	nt_gout = nt_gout + Xspkhst*nim.spk_hist.coefs(:);
end

%% SET TENT-BASIS X-AXES FOR TARGET SUBUNITS
nim = NMMinitialize_upstreamNLs(nim,Xstims,targets,nan,nan,nim.NL_tb_params.edge_p,nim.NL_tb_params.n_bfs,nim.NL_tb_params.space_type);
    
%% COMPUTE NEW X-MATRIX OUT OF TENT-BASIS OUTPUTS
XNL = zeros(NT,Ntargets*n_tbfs); %initialize X matrix which is for the NL BFs of each module
for ii = 1:Ntargets %for each module
		tar = targets(ii);
    % The output of the current model's internal filter projected onto the tent basis representation
		if isempty(Gmults{tar})
			tbf_out = nim.mods(tar).sign * tb_rep(gint(:,tar),nim.mods(tar).NLx);
		else
			%tbf_out = nim.mods(tar).sign * (Gmults{tar} .* tb_rep(gint(:,tar),nim.mods(tar).NLx));
			tbf_out = nim.mods(tar).sign * bsxfun(@times, tb_rep(gint(:,tar),nim.mods(tar).NLx), Gmults{tar} );
		end
    XNL(:,((ii-1)*n_tbfs + 1):(ii*n_tbfs)) = tbf_out; % assemble filtered NLBF outputs into X matrix
end

%% CREATE INITIAL PARAMETER VECTOR
% Compute initial fit parameters
initial_params = [];
for imod = targets(targets > 0)
	initial_params = [initial_params; nim.mods(imod).NLy']; %not incorporating the multiplier here because doing so messes with regularization
end

% Add in spike history coefs
if ismember(-1,targets) && spkhstlen > 0
	initial_params = [initial_params; nim.spk_hist.coefs];
end

initial_params(end+1) = nim.spk_NL_params(1); %add constant offset

%% CREATE L2 REG MATRICES IF USING REGULARIZATION
all_reg = zeros(Nmods,1);
for imod = 1:Nmods
	all_reg(imod) = nim.mods(imod).reg_params.lambda_NLd2;
end
if max(all_reg) > 0 % if using smoothness reg, generate the L2 matrix
	et = ones(n_tbfs,1);
	et([1 end]) = 0;
	L2_NLd2 = spdiags([et -2*et et], [-1 0 1], n_tbfs, n_tbfs)';
else
	L2_NLd2 = [];
end

%% PROCESS CONSTRAINTS
LB = []; UB = []; A = []; Aeq = []; % initialize constraint parameters

% Check for spike history coef constraints
if spkhstlen > 0 && ismember(-1,targets)
	% negative constraint on spk history coefs
	if nim.spk_hist.negCon == 1
		spkhist_inds = (Ntargets*n_tbfs + 1):(Ntargets*n_tbfs + spkhstlen);
		LB = -Inf*ones(size(initial_params));
		UB = Inf*ones(size(initial_params));
		UB(spkhist_inds) = 0;
	end
end

% Process NL monotonicity constraints, and constraints that the tent basis
% centered at 0 should have coefficient of 0 (eliminate y-shift degeneracy)
zvec = zeros(1,length(initial_params)); % indices of tent-bases centered at 0
for ii = 1:Ntargets
	cur_range = (ii-1)*n_tbfs + (1:n_tbfs);
    
	% For monotonicity constraint
	if nim.mods(targets(ii)).NLmon ~= 0
		for jj = 1:length(cur_range)-1
			cur_vec = zvec;
			cur_vec(cur_range([jj jj + 1])) = nim.mods(targets(ii)).NLmon*[1 -1];
			A = cat(1,A,cur_vec);
		end
	end
    
	% Constrain the 0-coefficient to be 0
	[~,zp] = find(nim.mods(targets(ii)).NLx == 0);
	if isempty(zp)
		error('Need one TB to be centered at 0')
	end
	cur_vec = zvec;
	cur_vec(cur_range(zp)) = 1;
	Aeq = cat(1,Aeq,cur_vec);
end
b = zeros(size(A,1),1);
beq = zeros(size(Aeq,1),1);

%% HANDLE OPTIMIZATION PARAMETERS
% default optimization parameters
if silent == 1
    optim_params.Display = 'off';
else
    optim_params.Display = 'iter';
end
optim_params.MaxFunEvals = 100*filtLen;
optim_params.MaxIter = 1e3;
optim_params.TolFun = 1e-6;
optim_params.TolX = 1e-6;
optim_params.HessUpdate = 'bfgs';
optim_params.GradObj = 'on';
optim_params.Algorithm = 'Active-set';

% Load in specified optimization parameters
if ~isempty(desired_optim_params)
	spec_fields = fieldnames(desired_optim_params);
	for i = 1:length(spec_fields)
		optim_params = setfield(optim_params,spec_fields{i},getfield(desired_optim_params,spec_fields{i}));
	end
end

% Constrained optimization
[params] = fmincon( @(K) NIMfit_upstreamNLs_internal(nim, K, Robs, XNL, Xspkhst, L2_NLd2, targets,nt_gout),...
    initial_params, A,b,Aeq,beq,LB,UB,[],optim_params);

%%
nim_out = nim;
nlmat = reshape(params(1:Ntargets*n_tbfs),n_tbfs,Ntargets); %take output K vector and restructure into a matrix of NLBF coefs, one for each module
nlmat_resc = nlmat;
for ii = 1:Ntargets;
    cur_pset = ((ii-1)*n_tbfs+1) : (ii*n_tbfs);
    thisnl = nlmat(:,ii); %NL coefs for current module
    
    cur_std = std(XNL(:,cur_pset)*thisnl);
    if rescale_NLs == 1
        %rescale so that the std dev of the subunit output is conserved
        thisnl = thisnl*nim.mods(targets(ii)).mod_norm/cur_std;
    else
        %otherwise adjust the model output std dev
        nim_out.mods(targets(ii)).mod_norm = cur_std; %
    end
    nim_out.mods(targets(ii)).NLy = thisnl';
    nlmat_resc(:,ii) = thisnl';
end

if ismember(-1,targets) && spkhstlen > 0
    nim_out.spk_hist.coefs = params((Ntargets*n_tbfs+1):(Ntargets*n_tbfs+spkhstlen));
end

%% IF RESCALING NLS, NEED TO RE-ESTIMATE OPTIMAL THETA
if rescale_NLs == 1
    resc_nlvec = nlmat_resc(:);
    
    new_g_out = XNL*resc_nlvec;
    G = nt_gout + new_g_out;
    if spkhstlen > 0 && ismember(-1,targets)
        G = G + Xspkhst*nim_out.spk_hist.coefs;
    end
    %if lin_dims > 0 && ismember(-2,targets)
    %    G = G + XLin*nim_out.kLin;
    %end
    
    initial_theta = params(end);
    opts.Display = 'off';
    opts.GradObj = 'on';
    opts.LargeScale = 'off';
    [best_const] = fminunc( @(K) internal_theta_opt(K,G,Robs,nim), initial_theta, opts);
    nim_out.spk_NL_params(1) = best_const;
else
    nim_out.spk_NL_params(1) = params(end);
end

%% COMPUTE FINAL LL Values
[LL, penLL,~,G] = NMMmodel_eval( nim_out, Robs, Xstims, Gmults, regmat_custom );
nim_out.LL_seq = cat(1,nim_out.LL_seq,LL);
nim_out.penLL_seq = cat(1,nim_out.penLL_seq,penLL);
nim_out.opt_history = cat(1,nim_out.opt_history,{'upstreamNLs'});

end

%% INTERNAL FUNCTIONS

function [LL, LLgrad] = NIMfit_upstreamNLs_internal(nim, params, Robs, Xnl, Xspkhst, L2_NLd2, targets, nt_gout)
%
% [LL, LLgrad] = NIMfit_upstreamNLs_internal(nim, params, Robs, Xnl, Xspkhst, L2_NLd2, targets, nt_gout)
%
% Internal function for computing LL and LLgradient with respect to the tent-basis coeffs

%% Useful params
Nmods = length(nim.mods);
Ntargets = sum(targets > 0);
n_tbfs = length(nim.mods(targets(1)).NLx);
NT = size(Xnl,1);
spkhstlen = nim.spk_hist.spkhstlen;

%% ESTIMATE GENERATING FUNCTIONS (OVERALL AND INTERNAL)
theta = params(end); %offset
G = theta + nt_gout;
G = G + Xnl*params(1:Ntargets*n_tbfs);

%add contribution from spike history filter
if spkhstlen > 0 && ismember(-1,targets)
	G = G + Xspkhst*params(Ntargets*n_tbfs + (1:spkhstlen));
end

%% Compute predicted firing rate
if strcmp(nim.spk_NL_type,'logexp')
    max_gbeta = 50; %to prevent numerical overflow
    bgint = G*nim.spk_NL_params(2); %g*beta
    expg = exp(bgint);
    too_large = (bgint > max_gbeta);
    r = nim.spk_NL_params(4) + nim.spk_NL_params(3)*log(1+expg); %alpha*log(1+exp(gbeta))
    r(too_large) = nim.spk_NL_params(4) + nim.spk_NL_params(3)*bgint(too_large); %log(1+exp(x)) ~ x in limit of large x
elseif strcmp(nim.spk_NL_type,'exp')
    expg = exp(G);
    r = expg;
		elseif strcmp(nim.spk_NL_type,'logistic')
    bgint = G*nim.spk_NL_params(2); %g*beta
    expg = exp(-bgint);
		r = nim.spk_NL_params(4) + 1./(1+expg); %1/(1+exp(-gbeta))  % note third param not used
elseif strcmp(nim.spk_NL_type,'linear')
    r = G;
else
    error('invalid spk nl');
end
% Enforce minimum predicted firing rate to avoid nan LLs
min_pred_rate = 1e-50;
if min(r) < min_pred_rate
    r(r < min_pred_rate) = min_pred_rate; %minimum predicted rate
end

%% COMPUTE LL and LL gradient
if strcmp(nim.spk_NL_type,'linear') % use MSE as cost function 
    Nspks = length(Robs);
    LL = -sum( (Robs - r).^2 );
elseif strcmp(nim.spk_NL_type,'logistic')
	% Bernouli likelihood = robs log r + (1-robs) log (1-r)
	Nspks = sum(Robs);
	LL = nansum( Robs.*log(r) + (1-Robs).*log(1-r) );
else
    Nspks = sum(Robs);
    LL = sum(Robs.* log(r) - r); %up to an overall constant
    %'residual' = (R/r - 1)*F'[] where F[.] is the spk NL
end

%'residual' = dLL/dr
if strcmp(nim.spk_NL_type,'logexp')
	residual = nim.spk_NL_params(3)*nim.spk_NL_params(2)*(Robs./r - 1) .* expg ./ (1+expg);
	residual(too_large) = nim.spk_NL_params(3)*nim.spk_NL_params(2)*(Robs(too_large)./r(too_large) - 1);
elseif strcmp(nim.spk_NL_type,'exp')
	residual = Robs - r;
elseif strcmp(nim.spk_NL_type,'logistic')
	% 'residual' = (R/r - (1-R)/(1-r))*F'[] where F[.] is the spk NL
	residual = nim.spk_NL_params(2)* (Robs./r - (1-Robs)./(1-r)) .* expg ./ (1+expg).^2;
elseif strcmp(nim.spk_NL_type,'linear')
	residual = 2*(Robs - r);    
else
	error('Unsupported spiking NL')
end

%initialize LL gradient
LLgrad = zeros(length(params),1);

LLgrad(1:Ntargets*n_tbfs) = residual'*Xnl;

% Calculate derivatives with respect to constant term (theta)
LLgrad(end) = sum(residual);

% Calculate derivative with respect to spk history filter
if spkhstlen > 0 && ismember(-1,targets)
    LLgrad((Ntargets*n_tbfs+1):(Ntargets*n_tbfs+spkhstlen)) = residual'*Xspkhst;
end

%% COMPUTE L2 PENALTIES AND GRADIENTS
smooth_penalty = zeros(Ntargets,1);
LLgrad_pen = zeros(size(LLgrad));
for n = 1:Ntargets
    cur_NLy = params((n-1)*n_tbfs + (1:n_tbfs));
    if nim.mods(targets(n)).reg_params.lambda_NLd2 > 0
        smooth_penalty(n) = nim.mods(targets(n)).reg_params.lambda_NLd2*sum((L2_NLd2 * cur_NLy).^2);
        
        cur_grad_pen = 2*nim.mods(targets(n)).reg_params.lambda_NLd2*(L2_NLd2' * L2_NLd2 * cur_NLy);
        LLgrad_pen((n-1)*n_tbfs + (1:n_tbfs)) = cur_grad_pen;
    end
end

LL = LL - sum(smooth_penalty);
LLgrad = LLgrad - LLgrad_pen;

%% CONVERT TO NEGATIVE LLS AND NORMALIZE BY NSPKS
LL = -LL/Nspks;
LLgrad = -LLgrad/Nspks;
end

%%
function [LL,grad] = internal_theta_opt(theta,G,Robs,nim)

G = G + theta;
if strcmp(nim.spk_NL_type,'logexp')
    max_gbeta = 50; %to prevent numerical overflow
    bgint = G*nim.spk_NL_params(2); %g*beta
    expg = exp(bgint);
    too_large = (bgint > max_gbeta);
    r = nim.spk_NL_params(3)*log(1+expg); %alpha*log(1+exp(gbeta))
    r(too_large) = nim.spk_NL_params(3)*bgint(too_large); %log(1+exp(x)) ~ x in limit of large x
elseif strcmp(nim.spk_NL_type,'exp')
    expg = exp(G);
    r = expg;
elseif strcmp(nim.spk_NL_type,'linear')
    r = G;    
else
    error('invalid spk nl');
end
%enforce minimum predicted firing rate to avoid nan LLs
min_pred_rate = 1e-50;
if min(r) < min_pred_rate
    r(r < min_pred_rate) = min_pred_rate; %minimum predicted rate
end

%% COMPUTE LL and LL gradient
if strcmp(nim.spk_NL_type,'linear') % use MSE as cost function 
    Nspks = length(Robs);
    LL = -sum( (Robs - r).^2 );
else
    Nspks = sum(Robs);
    LL = sum(Robs.* log(r) - r); %up to an overall constant
    %'residual' = (R/r - 1)*F'[] where F[.] is the spk NL
end

%'residual' = dLL/dr
if strcmp(nim.spk_NL_type,'logexp')
    residual = nim.spk_NL_params(3)*nim.spk_NL_params(2)*(Robs./r - 1) .* expg ./ (1+expg);
    residual(too_large) = nim.spk_NL_params(3)*nim.spk_NL_params(2)*(Robs(too_large)./r(too_large) - 1);
elseif strcmp(nim.spk_NL_type,'exp')
    residual = Robs - r;
elseif strcmp(nim.spk_NL_type,'linear')
    residual = 2*(Robs - r);        
else
    error('Unsupported spiking NL')
end

grad = sum(residual);


LL=-LL/Nspks;
grad=-grad'/Nspks; %'

end


