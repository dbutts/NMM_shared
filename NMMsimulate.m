function [pred_rate, spks, G, gint, fgint] = NMMsimulate( nim, Xstims, Gmults, Nreps)
%
% Usage: [pred_rate spks G gint fgint] = NMMsimulate( nim, Xstim, <Gmults>, <Nreps> )
%         or
%        [rt rSPK G gint fgint] = NMMsimulate( nim, Xstim, <Gmults>, Robs )
%
% If Robs is input in place of Nreps, then instead of generating spikes,
% this function will generate the probability of spiking in place of spikes
% give the observed firing
%
% Simulates spike trains of the specified model, also returns other useful outputs
%
% INPUTS:
%   nim: model structure
%   Robs: binned spikes
%   Xstim: time-embedded stimulus mat
%
% OUTPUTS:
%   pred_rate: predicted firing rate without spike history term. (unitless, so divide by dt to get Hz).
%       Note that if there is a spike history term, it does not take into account, and insteads 'spks' 
%       output should be used to estimate rate (PSTH-style)
%   spks: spike times generated from the simulation. Multiple repeats are stored in one list, with
%       a '-1' separating each repeat.
%   G: generating function (output of the model before the spk NL)
%   gint: TxNmods matrix of the output of each subunit's stimulus filter (before applying upstream NL)
%   fgint: TxNmods matrix of the output of each subunit (after applying upstream NL)
%					Note: if spikes input, last term of each is the spike-history	term output

if (nargin < 4) || isempty(Nreps)
	Nreps = 1;
else
	if length(Nreps) > 1  % then Robs instead of Nreps
		Robs = Nreps;
		Nreps = 1;
	else
		Robs = [];
	end
end

%% Process Xstims (in case multiple Xstims)
if ~iscell(Xstims)
	tmp = Xstims;
	clear Xstims
	Xstims{1} = tmp;
end

%% Key parameters
NT = size(Xstims{1},1);

Nmods = length(nim.mods);
spkhstlen = nim.spk_hist.spkhstlen;
dt = nim.stim_params.dt;

if (nargin < 3) || isempty(Gmults)
	Gmults{Nmods} = [];
end

%% CREATE L2 REGULARIZATION MATRICES
%L2_mats = create_L2_matrices_NMM( nim );

%% ESTIMATE GENERATING FUNCTIONS (OVERALL AND INTERNAL)
theta = nim.spk_NL_params(1); %offset
G = theta + zeros(NT,1); %initialize overall generating function G

%Kmat = [nim.mods(:).filtK];
%gint = Xstim*Kmat; %subunit generating functions

gint = nan(NT,Nmods);
fgint = nan(NT,Nmods);

for imod = 1:Nmods
	
	gint(:,imod) = Xstims{nim.mods(imod).Xtarget} * nim.mods(imod).filtK;
		
	% Process subunit g's with upstream NLs
	if strcmp(nim.mods(imod).NLtype,'nonpar')
		fgint(:,imod) = piecelin_process( gint(:,imod), nim.mods(imod).NLy, nim.mods(imod).NLx );
	elseif strcmp(nim.mods(imod).NLtype,'quad')
		fgint(:,imod) = (gint(:,imod)-nim.mods(imod).NLx).^2;
	elseif strcmp(nim.mods(imod).NLtype,'lin')
		fgint(:,imod) = gint(:,imod);
	elseif strcmp(nim.mods(imod).NLtype,'threshlin')
		fgint(:,imod) = gint(:,imod)-nim.mods(imod).NLx;
		fgint( fgint(:,imod)<0, imod ) = 0;
	elseif strcmp(nim.mods(imod).NLtype,'threshP')
		fgint = fgint-nim.mods(imod).NLx;
		fgint(fgint < 0) = 0;
		fgint = fgint.^(nim.mods(imod.NLy));
	else
		error('Invalid internal NL');
	end
    
	% Multiply by weight (and multiplier, if appl) and add to generating function
	if isempty(Gmults{imod})
		fgint(:,imod) = fgint(:,imod) * nim.mods(imod).sign;
	else
		fgint(:,imod) = (fgint(:,imod).*Gmults{imod}) * nim.mods(imod).sign;
	end
	G = G + fgint(:,imod);
end

% Calculate predicted rate (without spike history)
%logexp = 0;
if strcmp(nim.spk_NL_type,'logexp')
	max_gbeta = 50; %to prevent numerical overflow
	bgint = G*nim.spk_NL_params(2); %g*beta
	%expg = exp(bgint);
	too_large = (bgint > max_gbeta);
	pred_rate = nim.spk_NL_params(3)*log(1+exp(bgint)); %alpha*log(1+exp(gbeta))
	pred_rate(too_large) = nim.spk_NL_params(3)*bgint(too_large); %log(1+exp(x)) ~ x in limit of large x
	%logexp = 1;
elseif strcmp(nim.spk_NL_type,'logistic')
	pred_rate = nim.spk_NL_params(4) + 1./(1+exp(-G*nim.spk_NL_params(2))); %1/(1+exp(-gbeta))  % note third param not used
elseif strcmp(nim.spk_NL_type,'exp')
	%expg = exp(G);
	pred_rate = exp(G);
elseif strcmp(nim.spk_NL_type,'linear')
	pred_rate = G;
else
	error('invalid spk nl');
end

% Simulate Poisson spikes if Nreps
if isempty(Robs)
	spks = [];
	if spkhstlen == 0
		% no spike history term, just need Poisson generator
		for rep = 1:Nreps
			spkstemp = find(rand(NT,1) < pred_rate)*dt - dt/2; 
			spks = [spks' spkstemp' -1]';
			% spks = [spks (find(rand(NT,1) < pred_rate)'*dt - dt/2 -1];
		end	
	else
	
		% then generating function affected by generated spikes
		Lh = nim.spk_hist.bin_edges(end); 
		h = zeros(1,Lh); % spike-history term
		for n = 1:nim.spk_hist.spkhstlen
			h(nim.spk_hist.bin_edges(n):(nim.spk_hist.bin_edges(n+1)-1)) = nim.spk_hist.coefs(n);
		end
	
		% Simulate over time (all reps at once)
		spkstemp = zeros(NT+Lh,Nreps);  % add buffer at beginning for spike history
		for t = 1:NT
			Gspkhist = ones(1,Nreps) * G(t) + h * spkstemp(Lh+t-(1:Lh),:);

			%if logexp > 0
			if strcmp(nim.spk_NL_type,'logexp')
				r = nim.spk_NL_params(3)*log(1+exp(Gspkhist*nim.spk_NL_params(2)));
			elseif strcmp(nim.spk_NL_type,'logistic')
				r = nim.spk_NL_params(4) + 1./(1+exp(-Gspkhist*nim.spk_NL_params(2))); %1/(1+exp(-gbeta))  % note third param not used
			elseif strcmp(nim.spk_NL_type,'exp')
				r = exp(Gspkhist);
			end
			spkstemp(t+Lh,:) = rand(1,Nreps) < r;
		end
	
		for n = 1:Nreps
			spks = [spks' ((find(spkstemp(:,n) > 0)'-Lh)*dt - dt/2) -1]';
		end
	end
	
	if Nreps == 1
		spks = spks(1:end-1); % take the -1 off the end if only one rep
	end

else
	
	% Generate actual firing rate
	if spkhstlen == 0
		spks = pred_rate;
	else
		
		% Add spike-history term's effect on G and recalculate rate
		Gspkhist = create_spkhist_Xmat( Robs, nim.spk_hist.bin_edges ) * nim.spk_hist.coefs(:);
		gint(:,end+1) = Gspkhist;
		fgint(:,end+1) = Gspkhist;
		G = G + Gspkhist;

		% Calculate predicted rate (without spike history)
		%logexp = 0;
		if strcmp(nim.spk_NL_type,'logexp')
			max_gbeta = 50; %to prevent numerical overflow
			bgint = G*nim.spk_NL_params(2); %g*beta
			%expg = exp(bgint);
			too_large = (bgint > max_gbeta);
			spks = nim.spk_NL_params(3)*log(1+exp(bgint)); %alpha*log(1+exp(gbeta))
			spks(too_large) = nim.spk_NL_params(3)*bgint(too_large); %log(1+exp(x)) ~ x in limit of large x
			%logexp = 1;
		elseif strcmp(nim.spk_NL_type,'logistic')
			spks = nim.spk_NL_params(4) + 1./(1+exp(-G*nim.spk_NL_params(2))); %1/(1+exp(-gbeta))  % note third param not used
		elseif strcmp(nim.spk_NL_type,'exp')
			%expg = exp(G);
			spks = exp(G);
		elseif strcmp(nim.spk_NL_type,'linear')
			spks = G;
		end
	end
end
 