%% AERO 557: FINAL PROJECT

    % Jason Beals
    % 03/12/2022
    
    clear
    close all
    clc
    
    addpath 'C:\Users\jason\OneDrive - Cal Poly\MATLAB\School\Functions\Math'
    
    % Create a genetic algorithim to solve Ackley's function:
    
%% Ackley's Funciton:

dSpace = 0.1;
[X,Y,Z] = ackley2Dplot(dSpace);

figure(1)
title('Ackley''s Function Visualization')
surfc(X,Y,Z)
grid on

%% Initial Population Creation: 

seed = input(sprintf('Please input the spawn seed:\n'));
rng(seed)
D = [-32.7,32.7]; % Define domain

% Set initial Population Size
N0 = 50; % Initial population size:
Ni = N0;
xi = D(1) + (D(2)-D(1)).*rand(Ni,1); % rand inside of D
yi = D(1) + (D(2)-D(1)).*rand(Ni,1); % " "

% Visualize Initial population:
figure(2)
contour(X,Y,Z)
title(sprintf('Initial Population, N = %d; Seed = %d.',N0,seed))
hold on
scatter(xi,yi,'r.','LineWidth',16)
hold off

%% Settings:

% Initial Population Size is defined above. This code can be edited to lose
% population counts over time, or gain, but at the moment it is configueed
% to maintain a constant population and is not able to modified with a
% single setting.

topX       = 10;        % How many will mate and be copied?
Xratio     = 0.5;       % How many will surive (Xratio = [0 1])
pc         = 0.75;      % How successful is mating? (pc = [0 1])
mc         = 0.1;      % How many will mutate?
crosspoint = 3;         % How many bits will stay the same?

maxit      = 50;        % How many iteration is acceptable?
fitTOL     = 50;       % How fit must the generation be? (50-120 is pretty tight).

 movie      = 'On';       % Toggle movie setings
% movie      = 'Off';

CrossoverSetting = 'Random';
% CrossoverSetting = 'Fixed'; % Uncomment to fix a crossover point as defined above by "crosspoint"

%% Enter Main Loop:
fitness = 5000;
it      = 0;
while fitness >= fitTOL
%% Iterator & Setup:
it = it + 1;
if it > maxit; break; end

% While Loop Selection:
if it == 1
%% Initial Selection:
% Use Natural Selection: ==> rank fitness from best to worst. 

%{
% global Xratio
% Xratio = 0.5; % Typical value, start with this.
%}

% Evaluate Initial Population's Fitness (Ackley's Function is def here). 
for i = 1:length(xi)
    in    = [xi(i),yi(i)];
    fi(i) = ackley(in);
end

Nkeep = Ni*Xratio;
% How much are we replenishing?:
Nr   = Ni - Nkeep; % Nreplenish

% Here we want to find the lowest cost, and we are trying to go to a global
% minimum so we can use mink:
fit      = zeros(Nkeep,3);
fit(:,3) = mink(fi,Nkeep); % Find the fittest members within the kept pop size
for i = 1:length(fit)
    % Find the population member corresponding to that fit score
    idx = find(fi==fit(i,3),1);
    fit(i,1) = xi(idx);
    fit(i,2) = yi(idx);
end
  
else
% Establish fit again:
fi       = NewPop(:,3); 
fit(:,3) = mink(NewPop(:,3),Nkeep); % Find the fittest members within the kept pop size
for i = 1:length(fit)
    % Find the population member corresponding to that fit score
    idx = find(fi==fit(i,3),1);
    fit(i,1) = NewPop(idx,1);
    fit(i,2) = NewPop(idx,2);
end
end % End initial iteration if loop

  
%% Crossover (for the discarded population members)

%{
% Inputs / Settings / Globals:
global topX
topX = 5; % Here, if crossover fails this is the memebrs that will be exactly copied.
          % This also sets which membera will be used as the "Males"
          % meaning that they mate with ALL of the remaining "Female"
          % members.

global pc crosspoint
pc = 0.75; % Probabliity that crossover is successful
crosspoint = 3; % Point in the strong of bits you cross to
CrossoverSetting = 'Random';
% CrossoverSetting = 'Fixed';
%}
% Convert to binary:
bl = length(dec2bin(max(abs(D)))); 
Xbin = dec2bin(abs(fit(:,1)),bl); LbinX = length(Xbin(1,:));
Ybin = dec2bin(abs(fit(:,2)),bl); LbinY = length(Xbin(1,:));
% Necessary for sign logic:
SIGN = fit(:,1:2) >= 0; % 1 indicates a positive value

for q = 1:Nkeep
    fitbin(q,:) = append(Xbin(q,:),Ybin(q,:));
end


% Now enter mating pool:
mult = 0;
for k = 1:Nr
  %% Continually loop through fit members to keep static population:
    if k>length(fit(:,1))
       mult = mult + 1;
    end
    % Index going through fit population:
    i = k - mult*k; 
    
  % Setup male and female:
    m  = randi(5);
    f  = randi(length(fit(:,1)));
    if m==f % If they are the same, choose a different female (guaranteed):
        f = topX + randi(length(fit(:,1))-topX);
    end
  %% Will mating occur?
   rng('shuffle')
   r = abs(rand(1));
    if r <= pc
      %% It's a child!:
          switch CrossoverSetting
              case 'Fixed'
                  % This will crossover at the specficied crosspoint
                  newbin{k,1}   = append(fitbin(f,1:crosspoint),fitbin(m,crosspoint+1:end)); 
              case 'Random'
                  % This randomizes the amount of genes shared during each
                  % succ
%                   rng('shuffle')
                  crosspoint = randi(length(fitbin(1,:))-1);
                  newbin{k,1}   = append(fitbin(f,1:crosspoint),fitbin(m,crosspoint+1:end)); 
          end         
     % Mating did not occur
    else % Then keep one of the top 5
        newbin{k,1}   = fitbin(m,:);  
    end
  %% Sign and Decimal Logic:
    % The bits only handled integers.... Now for sign and decimal logic:
    dec         = rand(1); % easy enough
    % Sign logic:
    momsign = SIGN(f,:); dadsign = SIGN(m,:); % 1 is a positive sign. 
    for z = 1:2
    if momsign(z) + dadsign(z) == 2 % If signs agree, take it. 
        sign(z) = 1;
    elseif momsign(z) + dadsign(z) == 0
        sign(z) = -1;
    else % Take dads sign if they don't agree
        if dadsign(z) == 1; sign(z) = 1; else; sign(z) = -1; end
    end
    
    % Throw in mutation
     r_sign = rand(1);
     if r_sign <= mc; sign(z) = -1*sign(z); end
    end
    % output to newbin:
        newbin{k,2} = dec;
        newbin{k,3} = sign; 
    
end

%% Now for Mutation:

%{
global mc
mc = 0.05; % Percentage of populations genes that will randomly mutate:
%}

Nmut = floor(mc*(Nr+Nkeep));

rng('shuffle')
for k = 1:Nmut
    
    i = randi(Nr+Nkeep);            % Population member to mutate
    j = randi(length(fitbin(1,:)));  % Bit to flip
    
    if i <= Nkeep % We are mutating a fit gene:
%        disp("Unmutated member"); disp(fitbin(i,:));
       bit = fitbin(i,j);
       switch bit
           case '0'
              fitbin(i,j) = '1';
           case '1'
              fitbin(i,j) = '0';
       end
%        disp("Mutated member"); disp(fitbin(i,:));

    else          % We are mutating a new gene:
       i = i - Nkeep; % make sure to reset index for the new vector
%        disp("Unmutated member"); disp(newbin{i,1}); 
       bit = newbin{i,1}(1,j);
       switch bit
           case '0'
              newbin{i,1}(1,j) = '1';
           case '1'
              newbin{i,1}(1,j) = '0';
       end
%        disp("Mutated member"); disp(newbin{i,1}); 
    end
end

%% Re-evalute new population:

% First birth the population members:
% decimals and signs their epigentics, hm? ;)
offspring = zeros(Nr,3);
for i = 1:Nr
    newXbin(i,:) = newbin{i,1}(1,1:LbinX);
    newYbin(i,:) = newbin{i,1}(1,LbinX+1:end);
    int          = [bin2dec(newXbin(i,:)),bin2dec(newYbin(i,:))];
    int          = [int(1)+newbin{i,2},int(2)+newbin{i,2}].*newbin{i,3};
    % Allocate to the offspring
    offspring(i,1) = int(1); offspring(i,2) = int(2);
    % Evaluate fitness of the new offspring:
    offspring(i,3) = ackley(int);
end

%% Save necessary variables:
NewPop          = [fit;offspring];
Gen(it).NewPop  = [fit;offspring];
if it == 1
Gen(it).Compare = [mink(NewPop(:,3),50) mink(fi,50)'];
else
Gen(it).Compare = [mink(NewPop(:,3),50) mink(fi,50)];
end
% Re-entry criteria:
fitness = sum(NewPop(:,3));
Gen(it).Fitness = fitness; 

% figure(2+it)
% contour(X,Y,Z)
% hold on
% scatter(NewPop(:,1),NewPop(:,2),'r.')
% title(sprintf('Generation = %d',it))
% hold off
% 
fprintf('Generation %d complete. Population fitness, J = %4.1f. Fittest member, Min = %4.1f.\n',it,fitness,min(NewPop(:,3)))



end % end main loop

figure(2+it)
contour(X,Y,Z)
hold on
scatter(NewPop(:,1),NewPop(:,2),'r.')


switch movie
    case 'On'
v = VideoWriter('Genetic_Algorithim.mp4') ;
open(v)

for q = 1:10
    figure(9)
    title(sprintf('Generation %d',q))
    contour(X,Y,Z)
    hold on
    scatter(Gen(q).NewPop(:,1),Gen(q).NewPop(:,2),'r.')
    hold off
    F(i) = getframe;
    writeVideo(v,F(i))
end
    case 'Off'
    figure(3)
    contour(X,Y,Z)
    hold on
    scatter(NewPop(:,1),NewPop(:,2),'r.')
    title("Final Generation")
end

% repeat = 0.1;
% movie(F,repeat,20)









function [X,Y,Z] = ackley2Dplot(dSpace,a,b,c)

% funciton handle:
if (nargin < 4)
    c = 2*pi;
end
if (nargin < 3)
    b = 0.2;
end
if (nargin < 2)
    a = 20;
end

% Create equal grid:
x = -32.7:dSpace:32.7;
y = x;
[X,Y] = meshgrid(x,y);

% Find vlaue at all spots (this is in 2D)
Z = zeros(size(X));
for i = 1:length(x)
for j = 1:length(y)
    in = [X(i,j),Y(i,j)]; 
    Z(i,j) = ackley(in,a,b,c);
end
end

% Ackley's in interesating fashion:
function [y] = ackley(xx, a, b, c)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% ACKLEY FUNCTION
%
% Authors: Sonja Surjanovic, Simon Fraser University
%          Derek Bingham, Simon Fraser University
% Questions/Comments: Please email Derek Bingham at dbingham@stat.sfu.ca.
%
% Copyright 2013. Derek Bingham, Simon Fraser University.
%
% THERE IS NO WARRANTY, EXPRESS OR IMPLIED. WE DO NOT ASSUME ANY LIABILITY
% FOR THE USE OF THIS SOFTWARE.  If software is modified to produce
% derivative works, such modified software should be clearly marked.
% Additionally, this program is free software; you can redistribute it 
% and/or modify it under the terms of the GNU General Public License as 
% published by the Free Software Foundation; version 2.0 of the License. 
% Accordingly, this program is distributed in the hope that it will be 
% useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
% of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
% General Public License for more details.
%
% For function details and reference information, see:
% http://www.sfu.ca/~ssurjano/
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% INPUTS:
%
% xx = [x1, x2, ..., xd]
% a = constant (optional), with default value 20
% b = constant (optional), with default value 0.2
% c = constant (optional), with default value 2*pi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

d = length(xx);

if (nargin < 4)
    c = 2*pi;
end
if (nargin < 3)
    b = 0.2;
end
if (nargin < 2)
    a = 20;
end

sum1 = 0;
sum2 = 0;
for ii = 1:d
	xi = xx(ii);
	sum1 = sum1 + xi^2;
	sum2 = sum2 + cos(c*xi);
end

term1 = -a * exp(-b*sqrt(sum1/d));
term2 = -exp(sum2/d);

y = term1 + term2 + a + exp(1);

end

end

function [y] = ackley(xx, a, b, c)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% ACKLEY FUNCTION
%
% Authors: Sonja Surjanovic, Simon Fraser University
%          Derek Bingham, Simon Fraser University
% Questions/Comments: Please email Derek Bingham at dbingham@stat.sfu.ca.
%
% Copyright 2013. Derek Bingham, Simon Fraser University.
%
% THERE IS NO WARRANTY, EXPRESS OR IMPLIED. WE DO NOT ASSUME ANY LIABILITY
% FOR THE USE OF THIS SOFTWARE.  If software is modified to produce
% derivative works, such modified software should be clearly marked.
% Additionally, this program is free software; you can redistribute it 
% and/or modify it under the terms of the GNU General Public License as 
% published by the Free Software Foundation; version 2.0 of the License. 
% Accordingly, this program is distributed in the hope that it will be 
% useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
% of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
% General Public License for more details.
%
% For function details and reference information, see:
% http://www.sfu.ca/~ssurjano/
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% INPUTS:
%
% xx = [x1, x2, ..., xd]
% a = constant (optional), with default value 20
% b = constant (optional), with default value 0.2
% c = constant (optional), with default value 2*pi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

d = length(xx);

if (nargin < 4)
    c = 2*pi;
end
if (nargin < 3)
    b = 0.2;
end
if (nargin < 2)
    a = 20;
end

sum1 = 0;
sum2 = 0;
for ii = 1:d
	xi = xx(ii);
	sum1 = sum1 + xi^2;
	sum2 = sum2 + cos(c*xi);
end

term1 = -a * exp(-b*sqrt(sum1/d));
term2 = -exp(sum2/d);

y = term1 + term2 + a + exp(1);

end


%% Graveyard:

% best = find(fi<=fit(topX),topX); % Find the topX population members.
% good = find(fi>fit(topX)&&fi<=max(fit),abs(length(fit)-topX));
% % Males will mate with females and each other.
% male(:,1) = xi(best,:); male(:,2) = yi(best,:);
% female(:,1) = xi(good,:); female(:,2) = yi(good,:);
