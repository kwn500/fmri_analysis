%% SIMULATE VOXEL DATA ANALYSIS
% In this script, I am going to:
% 1) Generate a hemodynamic response function
% 2) Simulate a BOLD response from a single voxel by convolving an event
%    sequence with an HRF
% 3) Compute responses to individual events by:
%    a) Deconvolution (gives an amplitude response per event)
%    b) GLM (gives a separate beta weight for each event)
% 4) Finally I will assess how adding noise to the voxel response affects
%    the reliability of 3a and 3b.

clear all;
close all

ynicInit spm8

%% 1: GENERATE HRF

% Make an HRF using the SPM double-gamma
tr = 3; 
% Like in our EoO data - TR is 3s
[h,p] = spm_hrf(tr);
% h is the HRF. We are using the defaults for p (the parameters of the 
% gamma that respresent things like delay in seconds etc. Check spm_hrf for
% more details.)
cm = GetConvolutionMatrix(h',112);
% Multiplying the convolution matrix by the event sequence does the same 
% thing as if we use the function conv to convolve h with the event 
% sequence - in theory.........

figure(1)
subplot(2,1,1)
plot(h)
xlabel('TRs')
ylabel('Amplitude')
title('SPM double-gamma')
subplot(2,1,2)
imagesc(cm)
colormap hot
title('Convolution matrix')

%% 2: GENERATE SYNTHETIC FMRI DATA

% Make a sequence of events
eventSeq = zeros(112,1);
% Specify the timing of left and right eye events
lTimes = [10,16,25,32,70,80,100];
rTimes = [2,12,18,42,52,64,98];
% Make each event have different amplitudes
eventSeq(lTimes) = 2;
eventSeq(rTimes) = 1;

% Simulate a scan by convolving the event sequence with the HRF
tSeries = conv(eventSeq,h);
% Or do the same thing using the convolution matrix
cDat = (eventSeq')*cm;
 
figure(2)
plot(tSeries,'r')
hold on
plot(cDat,'b')
plot(eventSeq,'g')
legend({'tSeries using conv','cDat using convolution matrix','Events'})
xlabel('TRs')
ylabel('Response amplitude')
title('Convolved event sequence')
% I think the reason why these look different is because the hrf does not
% begin at '1' but at '2' (something with matlab not doing zero indexing).
% The conv matrix also doesn't extend beyond the number of specified TRs so
% it doesn't add a 'tail' like conv does. I could fix this but clearly the
% convlution matrix does what we want it to do anyway :)

%% 3a: DECONVOLUTION

% Here's how we would deconvolve this timecourse using the convolution
% matrix (note: we can only do this if we know the shape of the HRF. We are
% going to assume the double-gamma is good enough for us. Of course for 
% this example we know that it is seeing as this is what we used in the 
% first place to generate the convolved timecourse.)

% here i think cDat is representing an example fMRI response, and cm is the
% convoluted design matrix (hrf and design matrix). 
dcSeries = cDat/cm;
% and this gives us the pure event sequence.


% The result is the pure events sequence
figure(3)
plot(dcSeries,'b');
hold on 
plot(eventSeq, 'r');
xlabel('TRs')
ylabel('Response amplitude')
title('Deconvolved timecourse')

% We can check this using the function norm
res = norm(dcSeries' - eventSeq);
disp(sprintf('\n The difference between the deconvolved timecourse and the original event sequence is %d',res));

% Now that we have this, we want to get back the amplitudes of our
% different events. We know what the timing of these events was (in lTimes
% and rTimes)... So we can just use those times to index into the dcSeries
% variable and pull out the amplitude value at that point.
lAmps = dcSeries(:,lTimes);
rAmps = dcSeries(:,rTimes);

% We know what these amplitudes should be, because we set them ourselves.
% So let's compare these vectors.

% i think here, when we're saying amplitudes, we mean the numerical code of
% the condition in the event file? 
origLAmps = eventSeq(lTimes,:);
origRAmps = eventSeq(rTimes,:);

resL = norm(lAmps-origLAmps');
resR = norm(rAmps-origRAmps');
disp(sprintf('\n The difference between the original amplitudes and the deconvolved amplitudes is %d for left events and %d for right events', resL, resR));

%% 3b: GLM WITH CONV

% Let's try the same kind of thing but with a GLM analysis. We want to end
% up with separate betas for each event - regardless of event type.

% First, get the design matrix
eventTimes = sort(horzcat(lTimes,rTimes));
timePoints = size(cDat,2);

% Make a design matrix
events = zeros(timePoints,length(eventTimes));
for thisEvent = 1:length(eventTimes)
    events(eventTimes(thisEvent),thisEvent) = 1;
end

designMatrix=[ones(timePoints,1),events];

figure(4)
imagesc(designMatrix)
title('Design Matrix')

% Convolve this with our HRF - do this is the same way as how we generated
% cDat! - not using matlab 'conv' but multiply by convolution matrix...
% although here i'm pretty sure they're using the matlab function convolve.
for thisEv = 1:size(designMatrix,2)
    designMatrix(:,thisEv) = conv(designMatrix(:,thisEv),h,'same');
end

designMatrix(:,1)=1;

% the design matrix now models the hdr for each event. 
figure(5)
imagesc(designMatrix)
title('Convolved Design Matrix - conv')

% We can now solve the GLM to get back estimated beta values for each point
% b=((X'X)^-1)X'*Y, where X is the design matrix and Y is the fMRI response

% is beta how well the modelled hdr fits the fmri response? 

X = designMatrix;
Y = conv(eventSeq,h,'same'); % We need to make sure we are convolving the 
% timeseries data in the same way as how we convolved the design matrix
% with the HRF. So I've used 'conv' again rather than multiplying by cm as
% I did in section 2. The 'same' flag makes sure that the dimensions don't
% change.

figure(6)
plot(Y,'b')
hold on
plot(cDat,'r')
legend({'Y','cDat'})
% WEIRD - I don't think we should use 'conv'

% i'm a bit lost at this point- what do the betas represent?
betas = inv((X'*X))*X'*Y % YAY. By the way - this is different from just 
% doing the pseudo-inverse because we also *Y at the end... But the
% inv(X'*X)*X is the deconvolution part. You are finding the inverse of
% matrix X (i.e. solving for X) and then multiplying that by the
% timecourse. This returns a value for each predictor in Y.

% Get original amplitudes that we set ourselves for different events
origEVs = eventSeq(find(eventSeq));

%How do they compare?
resBeta = norm(betas(2:end)-origEVs);
disp(sprintf('\n The difference between the original amplitudes and the betas is %d ', resBeta));
resBetaDAmps = norm(betas(2:end)-dcSeries(eventTimes)');
disp(sprintf('\n The difference between the betas and the deconvolved amps is %d ', resBetaDAmps));

%% 3b: ALTERNATIVE - GLM WITH CONVOLUTION MATRIX

% Exactly the same as above, except now I'm going to do the whole thing
% without using 'conv'.

eventTimes = sort(horzcat(lTimes,rTimes));
timePoints = size(cDat,2);
% Re-set the design matrix...
events = zeros(timePoints,length(eventTimes));
for thisEvent = 1:length(eventTimes)
    events(eventTimes(thisEvent),thisEvent) = 1;
end
designMatrix=[ones(timePoints,1),events];
% Convolve
for thisEv = 1:size(designMatrix,2)
    designMatrix(:,thisEv) = designMatrix(:,thisEv)'*cm;
end
figure(7)
imagesc(designMatrix)
title('Convolved Design Matrix - cm')
% Get betas
X = designMatrix;
Y = cDat';
betas = inv((X'*X))*X'*Y % YAY!
% Compare
origEVs = eventSeq(find(eventSeq));
resBeta = norm(betas(2:end)-origEVs);
disp(sprintf('\n The difference between the original amplitudes and the betas is %d ', resBeta));
resBetaDAmps = norm(betas(2:end)-dcSeries(eventTimes)');
disp(sprintf('\n The difference between the betas and the deconvolved amps is %d ', resBetaDAmps));

%% 4: SIMULATE NOISE

% Let's see how deconvolution and beta estimates do when we add increasing
% levels of noise to our timeseries.

% First, add noise
noiseSD = .25; % Set the standard deviation of the noise
rnoise = cDat+noiseSD*randn(1,timePoints);

figure(8)
plot(cDat,'b')
hold on
plot(rnoise,'r')
legend({'Original response','Noisy response'})
title('Adding noise, SD = .25')

% Deconvolved amps
dcSeries = rnoise/cm;
deconAmps = dcSeries(:,eventTimes);
resDAmps = norm(deconAmps'-origEVs);
disp(sprintf('\n The difference between the deconvolved amps and orginal amplitudes is %d ', resDAmps));

% Event betas
X = designMatrix; % Using the one we just made by *cm
Y = rnoise';
betas = inv((X'*X))*X'*Y
resBetas = norm(betas(2:end)-origEVs);
disp(sprintf('\n The difference between the betas and orginal amplitudes is %d ', resBetas));

% Let's do that in a loop and iteratively increase the noise levels
noiseLevels = linspace(0,1,100);

for thisNoiseLevel = 1:length(noiseLevels)
    
    rnoise = cDat+noiseLevels(thisNoiseLevel)*randn(1,timePoints);

    dcSeries = rnoise/cm;
    deconAmps = dcSeries(:,eventTimes);
    resDAmps(thisNoiseLevel) = norm(deconAmps'-origEVs);
    
    X = designMatrix; % Using the one we just made by *cm
    Y = rnoise';
    betas = inv((X'*X))*X'*Y;
    resBetas(thisNoiseLevel) = norm(betas(2:end)-origEVs);
    
end % Next noise level

figure(9)
plot(noiseLevels,resDAmps,'ro')
hold on
plot(noiseLevels,resBetas,'bx')
legend({'Deconvolved amplitudes','Event betas'})
title('Effect of increasing noise on recovered event responses')
xlabel('Standard deviation of noise')
ylabel('Error')
    