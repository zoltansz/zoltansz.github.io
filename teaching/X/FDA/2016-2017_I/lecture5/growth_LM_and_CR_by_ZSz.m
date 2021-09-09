%Landmark and continuous registration on the growth dataset.

%clear start:
  clear all;
  close all;

%%%%%%%%%%%%%%%%%%%%%%%
%Landmark registration:
%%%%%%%%%%%%%%%%%%%%%%%

%load data, monotone smoothing ([2]: Sec 5.4.2.2 adapted to Matlab):
N = 54; n = 31;
fid = fopen('hgtf.dat','rt');
height = reshape(fscanf(fid,'%f'),[n,N]); % 31x54 = number_of_observations x number_of_curves
fclose(fid);
age = [1:0.25:2, 3:8, 8.5:0.5:18]' % length(age) = 31
rng_age = [age(1), age(end)]; m = 6; B = length(age) + m - 2;

wbasis = create_bspline_basis(rng_age,B,m,age);
growfdPar = fdPar(fd(zeros(B,N),wbasis), 3, 10^(-0.5));
disp('**************Monotone smoothing****************');
[Wfd,betaf,hgtfhatfd] = smooth_monotone(age, height, growfdPar);


%unregistered acceleration curves and their mean:
accelfdUN = deriv_fd(hgtfhatfd, 2); %'UN' = unregistered
accelmeanfdUN = mean(accelfdUN);

%select the age of the center of the pubertal growth spurt for each subject; 
%=age (the last one) where acceleration curve crosses zero with negative slope.

%-----------------------
%manual labelling, save:

if 1 %set it to zero once the landmarks are marked
%agefine = linspace(age(1), age(end), 101); %if finer resolution is needed
PGSctr = zeros(N,1); %1 <-> age
for icase = 1 : N
   %if finer resolution is not needed:
      plot(accelfdUN(icase)); 
   %if finer resolution is needed:   
      %acc_i = eval_fd(agefine, accelfdUN(icase));
      %plot(agefine, acc_i); grid on
   %xlabel, ylabel, title:
      title(['Subject-',num2str(icase)]);
      xlabel('age');
      ylabel('height acceleration');
   [age_i, temp] = ginput(1); %'1' <-> get 1 point
   PGSctr(icase) = age_i; 
end
PGSctrmean = mean(PGSctr);
save('PGS.mat','PGSctr','PGSctrmean');

end %of 'if 0'

%load
load('PGS.mat','PGSctr','PGSctrmean');
%-----------------------

%set up the warping functions (4-order spline, minimal smoothing):
wbasisLM = create_bspline_basis(rng_age, 4, 3, [age(1),PGSctrmean,age(end)]); %3: does not matter here
WfdLM = fd(zeros(4,1), wbasisLM);
WfdParLM = fdPar(WfdLM,1,1e-12);

%do the landmark registration:
disp('*******************Landmark registration**************')
[accelfdLM, warpfdLM, WfdLM] = landmarkreg(accelfdUN, PGSctr, PGSctrmean, WfdParLM, true); %true <-> h_i (warping function): strictly monotone
accelmeanfdLM = mean(accelfdLM); %'LM' <-> LandMark

%comparsion of the unregistered and the landmark registered curves:
figure; 
   plot(accelfdUN); 
   axis([1,18,-4,4]);  %'axis' to focus on the relevant part
   xlabel('age'); ylabel('acceleration'); title('unregistered curves');  
figure; 
   plot(accelfdLM); 
   axis([1,18,-4,4]); title('LM-registered curves');
   xlabel('age'); ylabel('acceleration'); 

%visualization of the warping functions:
figure; plot(warpfdLM); title('warping functions'); xlabel('age'), ylabel('own age');

%visualize the _original_ curves:
heightLM = zeros(length(age),N);
for icase = 1:N
     hi_at_age = eval_fd(age, warpfdLM(icase));    
     heightLM(:,icase) = eval_fd(hi_at_age, hgtfhatfd(icase));%a few NaNs since hi can map age outside of [1,18], i.e, the domain of 'hgtfhatfd'
end
figure; 
   plot(hgtfhatfd);    
   axis([1,18,60,190]); %'axis' to focus on the relevant part
   title('height curves');
   xlabel('age'); ylabel('height');
figure; 
   plot(age,heightLM); 
   axis([1,18,60,190]); %'axis' to focus on the relevant part
   title('height curves after acceleration registration'); 
   xlabel('age'); ylabel('height');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Continuous registration (initialized with the result of the landmark registration):
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

wbasisCR  = create_bspline_basis([1,18], 15, 5);  %'CR'<-> Continuous Registration
WfdParCR = fdPar(fd(zeros(15,N),wbasisCR), 1, 1); %'1,1' <->L = D^1, \lambda = 1
[accelfdCR, warpfdCR, WfdCR] = register_fd(mean(accelfdLM), accelfdLM, WfdParCR);

%visualize the CR-registered acceleration curves:
figure; 
   plot(accelfdCR); 
   axis([1,18,-4,4]);  %'axis' to focus on the relevant part
   xlabel('age'); ylabel('height'); title('LM-registered curves'); 

%compare the unregistered, LM-registered and LM+CR-registered mean curves:
figure;
accelmeanfdCR = mean(accelfdCR);
plot(age,eval_fd(age,accelmeanfdUN),'b', age,eval_fd(age,accelmeanfdLM),'g', age,eval_fd(age,accelmeanfdCR),'r');
%title, legend:
   title('mean curves'); legend({'unregistered mean','LM-registered mean','LM+CR-registered mean'})
   xlabel('age'); ylabel('acceleration');
