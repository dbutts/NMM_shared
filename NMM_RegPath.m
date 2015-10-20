function [bestfit,L2best] = NMM_RegPath( fit0, Robs, Xs, Uindx, XVindx, targets, L2s, lambdaID, fitparams )
%
% Usage: [bestift L2best] = NMM_RegPath( fit0, Robs, Xs, Uindx, XVindx, <targets>, <L2s>, <lambdaID>, <fitparams> )
%
% Rudimentary function to fit filters or non-parametric nonlinearities for range of reg values. This function
% is not necessarily robust or tested outside of particular situations.
%
% lambda-ID refers to which lamdba to do regularization path: default = 'lambda_d2T'
%   -- other main options include 'lambda_NLd2' and 'lambda_L1'

if nargin < 6
  targets = [];
end
%if (nargin < 7) || isempty(L2s)
%  L2s = [0 0.5 1 2 4 8 12 16 20 40 80 120 160 200 400 800 1200 2000 4000];
%end
if nargin < 7
  L2s = [];
end
if (nargin < 8) || isempty(lambdaID)
  lambdaID = 'lambda_d2T';
end

if nargin < 9
  fitparams = [];
	silent = 1;
else
	if isfield(fitparams,'silent')
		silent = fitparams.silent;
	else
		silent = 0;
	end
end
silent = 0;
fitparams.silent = 1;

Nreg = length(L2s);

if isempty(targets)
  targets = 1:length(fit0.mods);
end
fitparams.targets = targets;
fitparams.optim_params.optTol = 1e-6;
fitparams.optim_params.progTol = 1e-10;


if isempty(L2s)
  
  %% Do order-of-mag reg first
  L2s = [0 0.1 1.0 10 100 1000 10000 1e5];
  if ~silent
    fprintf( 'Order-of-magnitude L2 reg path: targets =' )
    disp(sprintf( ' %d', targets ))
  end
  LLregs = zeros(length(L2s),1);
  for nn = 1:length(L2s)
    regfit = NMMadjust_regularization( fit0, targets, lambdaID, L2s(nn) );
    if strcmp( lambdaID, 'lambda_NLd2' )
      regfit = NMMfit_upstreamNLs( regfit, Robs, Xs, [], Uindx, 1, 1 );
    else
      regfit = NMMfit_filters( regfit, Robs, Xs, [], Uindx, fitparams );
    end
    fitsaveM{nn} = regfit;
    [LL,LLnull] = NMMeval_model( regfit, Robs, Xs, [], XVindx );
    LLregs(nn) = LL-LLnull;
    if ~silent
      fprintf( '  %7.1f: %f\n', L2s(nn), LLregs(nn) )
    end
    if nn > 2
      if (LLregs(nn) < LLregs(nn-1)) && (LLregs(nn) < LLregs(nn-2))
        nn = length(L2s)+1;
      end
    end
  end
  [~,bestord] = sort(LLregs);
  % Check to see if contiguous
  if abs(bestord(end)-bestord(end-1)) == 1
    loweredge = min(bestord(end-1:end));
  else  % otherwise go the right direction on the best
    if bestord(end) > bestord(end-1)
      loweredge = bestord(end)-1;
    else
      loweredge = bestord(end);
    end
  end
  mag = L2s(loweredge);
  LLbounds = LLregs(loweredge+[0 1]); 

  % Zoom in on best regularization
  if mag == 0
    L2s = [0 0.01 0.02 0.04 0.1];
  else
    L2s = mag*[1 2 4 6 8 10];
  end
  Nreg = length(L2s);
  LLregs = zeros(Nreg,1);
  LLregs(1) = LLbounds(1);    LLregs(end) = LLbounds(2);
  fitsave{1} = fitsaveM{loweredge};
  fitsave{Nreg} = fitsaveM{loweredge+1};

  if ~silent
    fprintf( 'Zooming in on L2 reg path (%0.1f-%0.1f):\n', mag, mag*10 )
  end

  for nn = 2:(Nreg-1)
    regfit = NMMadjust_regularization( fit0, targets, lambdaID, L2s(nn) );
    if strcmp( lambdaID, 'lambda_NLd2' )
      regfit = NMMfit_upstreamNLs( regfit, Robs, Xs, [], Uindx, 1, 1 );
    else
      regfit = NMMfit_filters( regfit, Robs, Xs, [], Uindx, fitparams );
    end
    fitsave{nn} = regfit;

    [LL,LLnull] = NMMeval_model( regfit, Robs, Xs, [], XVindx );
    LLregs(nn) = LL-LLnull;
    if ~silent
      fprintf( '  %6.1f: %f\n', L2s(nn), LLregs(nn) )
    end
    if nn > 2
      if (LLregs(nn) < LLregs(nn-1)) && (LLregs(nn) < LLregs(nn-2))
        nn = length(L2s)+1;
      end
    end
  end
  
else
  
  %% Use L2 list estalished in function call
  if ~silent
    fprintf( 'L2 reg path (%d): targets =', Nreg )
    disp(sprintf( ' %d', targets ))
  end
  LLregs = zeros(Nreg,1);

  for nn = 1:Nreg
    regfit = NMMadjust_regularization( fit0, targets, lambdaID, L2s(nn) );
    if strcmp( lambdaID, 'lambda_NLd2' )
      regfit = NMMfit_upstreamNLs( regfit, Robs, Xs, [], Uindx, 1, 1 );
    else
      regfit = NMMfit_filters( regfit, Robs, Xs, [], Uindx, fitparams );
    end
    fitsave{nn} = regfit;
    [LL,LLnull] = NMMeval_model( regfit, Robs, Xs, [], XVindx );
    LLregs(nn) = LL-LLnull;
  
    if ~silent
      fprintf( '  %5.1f: %f\n', L2s(nn), LLregs(nn) )
    end
  end

end

[~,bestnn] = max(LLregs);
L2best = L2s(bestnn);
bestfit = fitsave{bestnn};
