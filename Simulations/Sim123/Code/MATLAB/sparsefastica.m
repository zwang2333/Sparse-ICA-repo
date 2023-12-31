function [A, W, time] = sparsefastica(X, whiteningMatrix, dewhiteningMatrix, approach, ...
			numOfIC, finetune, myy, stabilization, ...
			epsilon1, epsilon2, maxNumIterations, maxFinetune, initState, ...
			guess, sampleSize, displayMode, displayInterval, ...
			s_verbose, sigma, rou, func, round_2nd)
% sparsefastica: Main algorithm of SparseFastICA.
% This code is modified based on Hugo et al.'s fpica code, please see more details from http://research.ics.aalto.fi/ica/fastica/
% Copyright (c) State Key Laboratory of Cognitive Neuroscience and Learning @ Beijing Normal University
%	            Written by Ruiyang Ge 
% Mail to Authors: ruiyangge@hotmail.com
%% Parameters
% Perform independent component analysis using Hyvarinen's fixed point
% algorithm. Outputs an estimate of the mixing matrix A and its inverse W.
%
% whitesig                              :the whitened data as row vectors
% whiteningMatrix                       :can be obtained with function whitenv
% dewhiteningMatrix                     :can be obtained with function whitenv
% approach      [ 'symm' | 'defl' ]     :the approach used (deflation or symmetric); we use symmetric way in SparseFastICA.
% numOfIC       [ 0 - Dim of whitesig ] :number of independent components estimated
% finetune      [same as g + 'off']     :the nonlinearity used in finetuning.
% myy                                    :step size in stabilized algorithm
% stabilization [ 'on' | 'off' ]        :if myy < 1 then automatically on
% epsilon1                               :stopping criterion for the first step of SparseFastICA
% epsilon2                               :stopping criterion for the second step of SparseFastICA
% maxNumIterations                      :maximum number of iterations in the first step of SparseFastICA
% maxFinetune                           :maximum number of iteretions for finetuning
% initState     [ 'rand' | 'guess' ]    :initial guess or random initial state. See below
% guess                                 :initial guess for A. Ignored if initState = 'rand'
% sampleSize    [ 0 - 1 ]               :percentage of the samples used in one iteration
% displayMode   [ 'signals' | 'basis' | :plot running estimate
%                 'filters' | 'off' ]
% displayInterval                       :number of iterations we take between plots
% s_verbose       [ 'on' | 'off' ]        :report progress in text format
% sigma                                 :original value of sigma
% rou                                   :decreasing parameter for rou at each loop
% fun           [ 'pow3' | 'tanh' |]    :the nonlinearity used; ; we use pow3 and tanh in SparseFastICA.
% round_2nd                             :the maximum number of iterations in the second step of SparseFastICA

%% EXAMPLE
%        [A, W, time] = fastica_sparse(whitesig, whiteningMatrix, dewhiteningMatrix, approach,
%        numOfIC, finetune,  mu, stabilization, epsilon1, epsilon2, 
%        maxNumIterations, maxFinetune, initState, guess, sampleSize,
%        displayMode, displayInterval, verbose, sigma, rou, func, round_2nd);
% 
%
% 
%  

% @(#)$Id: fpica.m,v 1.7 2005/06/16 12:52:55 jarmo Exp $
% 

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Global variable for stopping the ICA calculations from the GUI
global g_FastICA_interrupt;
if isempty(g_FastICA_interrupt)
  clear global g_FastICA_interrupt;
  interruptible = 0;
else
  interruptible = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Default values

if nargin < 3, error('Not enough arguments!'); end
[vectorSize, numSamples] = size(X);
if nargin < 22, round_2nd = 20; end
if nargin < 21, func = 'pow3'; end
if nargin < 20, rou = 0.965; end
if nargin < 19, sigma = 1; end
if nargin < 18, s_verbose = 'on'; end
if nargin < 17, displayInterval = 1; end
if nargin < 16, displayMode = 'on'; end
if nargin < 15, initState = 'rand'; end
if nargin < 14, sampleSize = 1; end
if nargin < 13, guess = 1; end
if nargin < 12, maxFinetune = 100; end
if nargin < 11, maxNumIterations = 1000; end
if nargin < 10, epsilon2 = 0.000001; end
if nargin < 9, epsilon1 = 0.0001; end
if nargin < 8, stabilization = 'on'; end
if nargin < 7, myy = 1; end
if nargin < 6, finetune = 'off'; end
if nargin < 5, numOfIC = vectorSize; end     % vectorSize = Dim
if nargin < 4, approach = 'defl'; end

a1 = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Checking the data

if ~isreal(X)
  error('Input has an imaginary part.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Checking the value for verbose

switch lower(s_verbose)
 case 'on'
  b_verbose = 1;
 case 'off'
  b_verbose = 0;
 otherwise
  error(sprintf('Illegal value [ %s ] for parameter: ''verbose''\n', s_verbose));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Checking the value for approach

switch lower(approach)
 case 'symm'
  approachMode = 1;
 case 'defl'
  approachMode = 2;
 otherwise
  error(sprintf('Illegal value [ %s ] for parameter: ''approach''\n', approach));
end
if b_verbose, fprintf('Used approach [ %s ].\n', approach); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Checking the value for numOfIC

if vectorSize < numOfIC
  error('Must have numOfIC <= Dimension!');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Checking the sampleSize
if sampleSize > 1
  sampleSize = 1;
  if b_verbose
    fprintf('Warning: Setting ''sampleSize'' to 1.\n');
  end  
elseif sampleSize < 1
  if (sampleSize * numSamples) < 1000
    sampleSize = min(1000/numSamples, 1);
    if b_verbose
      fprintf('Warning: Setting ''sampleSize'' to %0.3f (%d samples).\n', ...
	      sampleSize, floor(sampleSize * numSamples));
    end  
  end
end
if b_verbose
  if  b_verbose & (sampleSize < 1)
    fprintf('Using about %0.0f%% of the samples in random order in every step.\n',sampleSize*100);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Checking the value for nonlinearity.

switch lower(func)
 case 'pow3'
  gOrig = 10;
 case 'tanh'
  gOrig = 20;
 otherwise
  error(sprintf('Illegal value [ %s ] for parameter: ''g''\n', func));
end
if sampleSize ~= 1
  gOrig = gOrig + 2;
end
if myy ~= 1
  gOrig = gOrig + 1;
end

if b_verbose,
  fprintf('Used nonlinearity [ %s ].\n', func);
end

finetuningEnabled = 1;
switch lower(finetune)
 case 'pow3'
  gFine = 10 + 1;
 case 'tanh'
  gFine = 20 + 1;
 case 'off'
  if myy ~= 1
    gFine = gOrig;
  else 
    gFine = gOrig + 1;
  end
  finetuningEnabled = 0;
 otherwise
  error(sprintf('Illegal value [ %s ] for parameter: ''finetune''\n', ...
		finetune));
end

if b_verbose & finetuningEnabled
  fprintf('Finetuning enabled (nonlinearity: [ %s ]).\n', finetune);
end

switch lower(stabilization)
 case 'on'
  stabilizationEnabled = 1;
 case 'off'
  if myy ~= 1
    stabilizationEnabled = 1;
  else
    stabilizationEnabled = 0;
  end
 otherwise
  error(sprintf('Illegal value [ %s ] for parameter: ''stabilization''\n', ...
		stabilization)); 
end

if b_verbose & stabilizationEnabled
  fprintf('Using stabilized algorithm.\n');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Some other parameters
myyOrig = myy;
% When we start fine-tuning we'll set myy = myyK * myy
myyK = 0.01;
% How many times do we try for convergence until we give up.
failureLimit = 5;


usedNlinearity = gOrig;
stroke = 0;
notFine = 1;
long = 0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Checking the value for initial state.

switch lower(initState)
 case 'rand'
  initialStateMode = 0;
 case 'guess'
  if size(guess,1) ~= size(whiteningMatrix,2)
    initialStateMode = 0;
    if b_verbose
      fprintf('Warning: size of initial guess is incorrect. Using random initial guess.\n');
    end
  else
    initialStateMode = 1;
    if size(guess,2) < numOfIC
      if b_verbose
	fprintf('Warning: initial guess only for first %d components. Using random initial guess for others.\n', size(guess,2)); 
      end
      guess(:, size(guess, 2) + 1:numOfIC) = ...
					     rand(vectorSize,numOfIC-size(guess,2))-.5;
    elseif size(guess,2)>numOfIC
      guess=guess(:,1:numOfIC);
      fprintf('Warning: Initial guess too large. The excess column are dropped.\n');
    end
    if b_verbose, fprintf('Using initial guess.\n'); end
  end
 otherwise
  error(sprintf('Illegal value [ %s ] for parameter: ''initState''\n', initState));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Checking the value for display mode.

switch lower(displayMode)
 case {'off', 'none'}
  usedDisplay = 0;
 case {'on', 'signals'}
  usedDisplay = 1;
  if (b_verbose & (numSamples > 10000))
    fprintf('Warning: Data vectors are very long. Plotting may take long time.\n');
  end
  if (b_verbose & (numOfIC > 25))
    fprintf('Warning: There are too many signals to plot. Plot may not look good.\n');
  end
 case 'basis'
  usedDisplay = 2;
  if (b_verbose & (numOfIC > 25))
    fprintf('Warning: There are too many signals to plot. Plot may not look good.\n');
  end
 case 'filters'
  usedDisplay = 3;
  if (b_verbose & (vectorSize > 25))
    fprintf('Warning: There are too many signals to plot. Plot may not look good.\n');
  end
 otherwise
  error(sprintf('Illegal value [ %s ] for parameter: ''displayMode''\n', displayMode));
end

% The displayInterval can't be less than 1...
if displayInterval < 1
  displayInterval = 1;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if b_verbose, fprintf('Starting ICA calculation...\n'); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SYMMETRIC APPROACH
if approachMode == 1,

  % set some parameters more...
  usedNlinearity = gOrig;
  stroke = 0;
  notFine = 1;
  long = 0;
  
  A = zeros(vectorSize, numOfIC);  % Dewhitened basis vectors.
  if initialStateMode == 0
    % Take random orthonormal initial vectors.
     B = orth (randn (vectorSize, numOfIC));
  elseif initialStateMode == 1
    % Use the given initial vector as the initial state
    B = whiteningMatrix * guess;
  end
  
  BOld = zeros(size(B));
  BOld2 = zeros(size(B));
  
  tic; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Timing begin
  %epsilon1 = epsilon*100;
    
  % This is the actual fixed-point iteration loop.
  for round = 1:maxNumIterations + 1,
    if round == maxNumIterations + 1,
      fprintf('No convergence after %d steps\n', maxNumIterations);
      fprintf('Note that the plots are probably wrong.\n');
      if ~isempty(B)
	% Symmetric orthogonalization.
	B = B * real(inv(B' * B)^(1/2));

	W = B' * whiteningMatrix;
	A = dewhiteningMatrix * B;
      else
	W = [];
	A = [];
      end
      return;
    end
    
    if (interruptible & g_FastICA_interrupt)
      if b_verbose 
        fprintf('\n\nCalculation interrupted by the user\n');
      end
      if ~isempty(B)
	W = B' * whiteningMatrix;
	A = dewhiteningMatrix * B;
      else
	W = [];
	A = [];
      end
      return;
    end
    
    
    % Symmetric orthogonalization.
    B = B * real(inv(B' * B)^(1/2));
    
    % Test for termination condition. Note that we consider opposite
    % directions here as well.
    minAbsCos = min(abs(diag(B' * BOld)));
    minAbsCos2 = min(abs(diag(B' * BOld2)));
    
    if (1 - minAbsCos < epsilon1)
      if finetuningEnabled & notFine
        if b_verbose, fprintf('Initial convergence, fine-tuning: \n'); end;
        notFine = 0;
        usedNlinearity = gFine;
        myy = myyK * myyOrig;
        BOld = zeros(size(B));
        BOld2 = zeros(size(B));
	
      else
        if b_verbose, fprintf('Convergence after %d steps\n', round); end
	
        % Calculate the de-whitened vectors.
        A = dewhiteningMatrix * B;
        break;
      end
    elseif stabilizationEnabled
      if (~stroke) & (1 - minAbsCos2 < epsilon1) % Note from Zihang: change epsilon to epsilon1
	if b_verbose, fprintf('Stroke!\n'); end;
	stroke = myy;
	myy = .5*myy;
	if mod(usedNlinearity,2) == 0
	  usedNlinearity = usedNlinearity + 1;
	end
      elseif stroke
	myy = stroke;
	stroke = 0;
	if (myy == 1) & (mod(usedNlinearity,2) ~= 0)
	  usedNlinearity = usedNlinearity - 1;
	end
      elseif (~long) & (round>maxNumIterations/2)
	if b_verbose, fprintf('Taking long (reducing step size)\n'); end;
	long = 1;
	myy = .5*myy;
	if mod(usedNlinearity,2) == 0
	  usedNlinearity = usedNlinearity + 1;
	end
      end
    end
    
    BOld2 = BOld;
    BOld = B;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Show the progress...
    if b_verbose
      if round == 1
        fprintf('Step no. %d\n', round);
      else
        fprintf('Step no. %d, change in value of estimate: %.3g \n', round, 1 - minAbsCos);
      end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Also plot the current state...
    switch usedDisplay
     case 1
      if rem(round, displayInterval) == 0,
	% There was and may still be other displaymodes...
	% 1D signals
	icaplot('dispsig',(X'*B)');
	drawnow;
      end
     case 2
      if rem(round, displayInterval) == 0,
	% ... and now there are :-)
	% 1D basis
	A = dewhiteningMatrix * B;
	icaplot('dispsig',A');
	drawnow;
      end
     case 3
      if rem(round, displayInterval) == 0,
	% ... and now there are :-)
	% 1D filters
	W = B' * whiteningMatrix;
	icaplot('dispsig',W);
	drawnow;
      end
     otherwise
    end
    
    switch usedNlinearity
      % pow3
     case 10
       B = (X * (( X' * B) .^ 3)) / numSamples - 3 * B;
      
     case 11
      % optimoitu - epsilonin kokoisia eroja
      % t?m?on optimoitu koodi, katso vanha koodi esim.
      % aikaisemmista versioista kuten 2.0 beta3
      Y = X' * B;
      Gpow3 = Y .^ 3;
      Beta = sum(Y .* Gpow3);
      D = diag(1 ./ (Beta - 3 * numSamples));
      B = B + myy * B * (Y' * Gpow3 - diag(Beta)) * D;     
      
     case 12
      Xsub=X(:, getSamples(numSamples, sampleSize));
      B = (Xsub * (( Xsub' * B) .^ 3)) / size(Xsub,2) - 3 * B;
     case 13
      % Optimoitu
      Ysub=X(:, getSamples(numSamples, sampleSize))' * B;
      Gpow3 = Ysub .^ 3;
      Beta = sum(Ysub .* Gpow3);
      D = diag(1 ./ (Beta - 3 * size(Ysub', 2)));
      B = B + myy * B * (Ysub' * Gpow3 - diag(Beta)) * D;
      
      % tanh
     case 20
      hypTan = tanh(a1 * X' * B);
      B = X * hypTan / numSamples - ...
	  ones(size(B,1),1) * sum(1 - hypTan .^ 2) .* B / numSamples * ...
	  a1;  

     case 21
      % optimoitu - epsilonin kokoisia 
      Y = X' * B;
      hypTan = tanh(a1 * Y);
      Beta = sum(Y .* hypTan);
      D = diag(1 ./ (Beta - a1 * sum(1 - hypTan .^ 2)));
      B = B + myy * B * (Y' * hypTan - diag(Beta)) * D;
     case 22
      Xsub=X(:, getSamples(numSamples, sampleSize));
      hypTan = tanh(a1 * Xsub' * B);
      B = Xsub * hypTan / size(Xsub, 2) - ...
	  ones(size(B,1),1) * sum(1 - hypTan .^ 2) .* B / size(Xsub, 2) * a1;
     case 23
      % Optimoitu
      Y = X(:, getSamples(numSamples, sampleSize))' * B;
      hypTan = tanh(a1 * Y);
      Beta = sum(Y .* hypTan);
      D = diag(1 ./ (Beta - a1 * sum(1 - hypTan .^ 2)));
      B = B + myy * B * (Y' * hypTan - diag(Beta)) * D;

     otherwise
      error('Code for desired nonlinearity not found!');
    end
  end
  
  
%%%%%%%%%%%%%%%%%%%% 2nd loop
  %epsilon2 = epsilon;
%   minAbsCos = 0;
for round = 1:round_2nd
    % Symmetric orthogonalization.
    B = B * real(inv(B' * B)^(1/2));
    % Test for termination condition. Note that we consider opposite directions here as well.
%     minAbsCos = min(abs(diag(B' * BOld)));  

if (1 - minAbsCos > epsilon2)  
% if (1 - minAbsCos > 0)    
    
    BOld = B;  
    
if round ==1
    sigma = sigma * 1.000;%%%%%%% sigma is the initial value
else
    sigma = sigma * rou; %%%%%%%% rou is the decreasing parameter for rou at each loop
end

switch usedNlinearity
    
    case 10 %%%%%%%%%%% pow3
      U = X' * B;
      Usquared=U .^ 2;
      ex = exp(-1 * Usquared / (2*sigma.^2));
      gauss = -1*U.*ex/sigma.^2;
      dGauss = (Usquared/sigma.^4 - 1/sigma.^2).*ex;
     % B = (X * (( X' * B) .^ 3)) / numSamples - 3 * B +...
     %     lambda * (X * gauss - ones(size(B,1),1) * sum(dGauss).* B)/ numSamples;%%%%%%% SparseFast-pow3 with finite lambda
      B = (X * gauss - ones(size(B,1),1) * sum(dGauss).* B)/ numSamples;%%%%%%% SparseFast-pow3 with infinite lambda



    case 20  %%%%%%%%%%%%%%%%%% tanh
      U = X' * B;
      Usquared=U .^ 2;
      ex = exp(-1 * Usquared / (2*sigma.^2));
      gauss = -1*U.*ex/sigma.^2;
      dGauss = (Usquared/sigma.^4-1/sigma.^2).*ex;

      hypTan = tanh(a1 * X' * B);
     % B = X * hypTan / numSamples - ...
	 % ones(size(B,1),1) * sum(1 - hypTan .^ 2) .* B / numSamples * a1 +...
     %     lambda * (X * gauss - ones(size(B,1),1) * sum(dGauss).* B)/ numSamples;%%%%%%% 20121211 LS0; 
      B = (X * gauss - ones(size(B,1),1) * sum(dGauss).* B)/ numSamples;%%%%%%% 20121211 LS0; 
 

end
      fprintf('Step no. %d, change in value of estimate: %3f; sigma=%f\n', round, 1 - minAbsCos,sigma);
      
    % Symmetric orthogonalization.
    B = B * real(inv(B' * B)^(1/2));
    
    % Test for termination condition. Note that we consider opposite
    % directions here as well.
    minAbsCos = min(abs(diag(B' * BOld)));    
end

end

%%%%%%%%%%%%%%%%%%%% 2nd loop
  
  time = toc;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Timing end
  
  % Calculate ICA filters.
  W = B' * whiteningMatrix;
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Also plot the last one...
  switch usedDisplay
   case 1 
    % There was and may still be other displaymodes...
    % 1D signals
    icaplot('dispsig',(X'*B)');
    drawnow;
   case 2
    % ... and now there are :-)
    % 1D basis
    icaplot('dispsig',A');
    drawnow;
   case 3
    % ... and now there are :-)
    % 1D filters
    icaplot('dispsig',W);
    drawnow;
   otherwise
  end
  
else
    return;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% In the end let's check the data for some security
if ~isreal(A)
  if b_verbose, fprintf('Warning: removing the imaginary part from the result.\n'); end
  A = real(A);
  W = real(W);
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Subfunction
% Calculates tanh simplier and faster than Matlab tanh.
function y=tanh(x)
y = 1 - 2 ./ (exp(2 * x) + 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Samples = getSamples(max, percentage)
Samples = find(rand(1, max) < percentage);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculates norm_2 simplier and faster than Matlab norm(Data,2).
function y=norm2(x)
y = sqrt(sum(sum((x.^2))));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Subfunction
% Calculates norm_1 simplier and faster than Matlab norm(Data,1).
function y=norm1(x)
y = sum(sum(abs(x)));