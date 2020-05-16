%% Plot single traces at boundary position

clear % clear all variables in the workspace
close all % close any figure windows that might be open
clc

addpath(genpath('./lib'))
%% Part 0: Initialization
% Specify the dataset to analyze

% We need to specify the prefix of the project we wish to analyze
%Prefix = '2020-02-23-optoknirps_new_embryo6'; % WT
Prefix = '2020-02-21-optoknirps_new_embryo5'; % with export
%Prefix = '2020-02-23-optoknirps_new_embryo7'; % Y->N
%Prefix = '2020-03-04-optoknirps_new_embryo10'; % N->Y
%% 
% Variables depending on the data
%% 
% * Start/Final frames

% WT: embryo6
%start_frame = 40;
%final_frame = 149;

% with export: embryo5
start_frame = 20;
final_frame = 124;

% Y->N: embryo7
%start_frame = 55;
%final_frame = 157;

% N->Y: embryo 10
%start_frame = 10;
%final_frame = 120;

APbinRange = 35:47;
%% Part 1: Add working path and load the data
% We will initialize the paths that are used later

% Make figure path
FigPath = ['/Users/jiaxi.zhao/Dropbox/OptoDropbox/#analysis/202004 - perturbation analysis/with_export/figure/' Prefix];

mkdir(FigPath)

% Specify the position of DynamicsResults folder
DynamicsResultsFolder = '/Users/jiaxi.zhao/Dropbox/OptoDropbox';
%% 
% Now, load all the data

% Now we'll load the CompiledParticles dataset
% This generates a string that gives the path to CompiledParticles
LoadParticlePath = [DynamicsResultsFolder '/' Prefix '/CompiledParticles.mat']; 
LoadNucleiPath = [DynamicsResultsFolder '/' Prefix '/' Prefix '_lin.mat'];
%LoadComNucleiPath = [DynamicsResultsFolder '/' Prefix '/CompiledNuclei.mat']; 

%data.particle = load(LoadParticlePath);
%data.nuclei = load(LoadNucleiPath);
load(LoadNucleiPath)
load(LoadParticlePath);

ElapsedTime = ElapsedTime-ElapsedTime(start_frame);
Time = ElapsedTime(start_frame:final_frame);
%APbinID = data.particle.APbinID;
%% Part 2: Quality check on nuclei tracking

% Quality check and save the nuclei that passed the test
% Requires the nuclei to have continuous trace until lateral movement
X_pass = [];
Y_pass = [];

X_fail = [];
Y_fail = [];

sch_pass = [];

for i = 1:size(schnitzcells,2)
    index_s = find(schnitzcells(i).frames == start_frame);
    index_f = find(schnitzcells(i).frames == final_frame);
    fluo_temp = max(schnitzcells(i).Fluo(index_s:index_f,:),[],2);
    result = sum(isnan(fluo_temp));
    % quality check
    if ~isempty(index_s) && ~isempty(index_f) && (index_f-index_s == final_frame-start_frame)  ...
            && (result == 0)
        sch_pass = [sch_pass, i];
        X_pass = [X_pass, schnitzcells(i).cenx(index_s)];
        Y_pass = [Y_pass, schnitzcells(i).ceny(index_s)];
    else
        X_fail = [X_fail, schnitzcells(i).cenx(index_s)];
        Y_fail = [Y_fail, schnitzcells(i).ceny(index_s)];
    end
end
figure(1)
plot(X_pass,Y_pass,'o')
hold on
plot(X_fail,Y_fail,'o')
axis equal
xlim([0 1024])
ylim([0 256])
title('Nuclei that passed quality check')
%% Part 3: Compile the schnitzcells together

% Initialize storage
processed_data(1).xcoord = [];
processed_data(1).ycoord = [];
processed_data(1).APpos = [];

processed_data(1).schnitznum = [];
processed_data(1).NuclearFluo = [];
processed_data(1).SpotFluo = [];

% Compile all the nuclei in each frame and assign basic info
for i = 1:length(sch_pass)
    sch_num = sch_pass(i);
    index_s = find(schnitzcells(sch_num).frames == start_frame);
    index_f = find(schnitzcells(sch_num).frames == final_frame);
    
    for j = index_s:index_f
        sch_now = sch_pass(i);
        frame_now = schnitzcells(sch_now).frames(j);
        x_coord = schnitzcells(sch_now).cenx(j);
        y_coord = schnitzcells(sch_now).ceny(j);
        AP_pos = schnitzcells(sch_now).APpos(j);
        fluo = max(schnitzcells(sch_now).Fluo(j,:));
        try
            processed_data(frame_now).xcoord = [processed_data(frame_now).xcoord, x_coord];
            processed_data(frame_now).ycoord = [processed_data(frame_now).ycoord, y_coord];
            processed_data(frame_now).APpos = [processed_data(frame_now).APpos, AP_pos];
            processed_data(frame_now).schnitznum = [processed_data(frame_now).schnitznum, sch_now];
            processed_data(frame_now).SpotFluo = [processed_data(frame_now).SpotFluo, 0];
            processed_data(frame_now).NuclearFluo = [processed_data(frame_now).NuclearFluo fluo];  
        catch
            processed_data(frame_now).xcoord = x_coord;
            processed_data(frame_now).ycoord = y_coord;
            processed_data(frame_now).APpos = AP_pos;
            processed_data(frame_now).schnitznum = sch_now;
            processed_data(frame_now).SpotFluo = 0;
            processed_data(frame_now).NuclearFluo = fluo;
        end
    end
end

% Assign particle info to all the nuclei
for i = 1:size(CompiledParticles{1,1},2)
    schnitz_num = CompiledParticles{1,1}(i).schnitz;
    for j = 1:size(CompiledParticles{1,1}(i).Frame,2)
        frame = CompiledParticles{1,1}(i).Frame(j);
        if frame<=final_frame
            num = find(processed_data(frame).schnitznum==schnitz_num);
            processed_data(frame).SpotFluo(num) = CompiledParticles{1,1}(i).Fluo(j);
        end
    end
end

%% 
% Let's try to plot the data

frame_plot = 100;

xpos = processed_data(frame_plot).xcoord;
ypos = processed_data(frame_plot).ycoord;

% try to plot and test a bit...
pts = [xpos' ypos'];

fig = figure(2);
[v,c] = voronoin(double(pts));

for i = 1:length(c)
    if all(c{i}~=1)
    x = v(c{i},1);
    y = v(c{i},2);
    %a = processed_data(frame_plot).SpotFluo(i);
    a = processed_data(frame_plot).NuclearFluo(i);
    patch(x,y,a);
    colorbar
    caxis([0 9E5])
    end
end

axis equal
xlim([0 1024])
ylim([0 256])
% Convert the format: Store them into individual nuclei traces
nuclei_fluo_traces = zeros(length(sch_pass),final_frame);
spot_fluo_traces = zeros(length(sch_pass),final_frame);
x_pos = zeros(length(sch_pass),final_frame);
y_pos = zeros(length(sch_pass),final_frame);
AP_pos = zeros(length(sch_pass),final_frame);

for i = 1:length(sch_pass)
    for j = start_frame:final_frame
        nuclei_fluo_traces(i,j) = processed_data(j).NuclearFluo(i);
        spot_fluo_traces(i,j) = processed_data(j).SpotFluo(i);
        x_pos(i,j) = processed_data(j).xcoord(i);
        y_pos(i,j) = processed_data(j).ycoord(i);
        AP_pos(i,j) = processed_data(j).APpos(i);
    end
end
%% Part 4: Process traces

% We would like to combine everything together and interpolate values
% Define time point to interpolate
time_point_final = ElapsedTime(start_frame):1/6:ElapsedTime(final_frame);
x = ElapsedTime(1:final_frame);

nuclei_fluo_ins = zeros(length(sch_pass),length(time_point_final));
spot_fluo_ins= zeros(length(sch_pass),length(time_point_final));

nuclei_fluo_avr = zeros(length(sch_pass),length(time_point_final));
spot_fluo_avr = zeros(length(sch_pass),length(time_point_final));

for i = 1:length(sch_pass)
    
    spot_temp = spot_fluo_traces(i,:);
    nuclei_temp = nuclei_fluo_traces(i,:);

    % Interpolate traces
    spot_interp = interp1(x,spot_temp,time_point_final,'pchip');
    nuclei_interp = interp1(x,nuclei_temp,time_point_final,'pchip');

    % pmoving average
    %nuclei_interp_avr = movmean(nuclei_interp,13);
    %spot_interp_avr = movmean(spot_interp,13);
    nuclei_interp_avr = nuclei_interp;
    spot_interp_avr = spot_interp;

    spot_fluo_ins(i,:) = spot_interp;
    nuclei_fluo_ins(i,:) = nuclei_interp;

    % use integrated mRNA level as input
    spot_fluo_avr(i,:) = spot_interp_avr;
    nuclei_fluo_avr(i,:) = nuclei_interp_avr;

    %{
    %spot_temp = spot_fluo_output(i,:);
    if sum(spot_temp) > 0
        TraceFig = figure;
        subplot(1,2,1)
        plot(x,nuclei_temp,'-o',time_point_final,nuclei_interp,'-x');
        hold on
        plot(time_point_final,nuclei_interp_avr,'-','LineWidth',2)
        xlabel('Time (min)')
        ylabel('Fluorescence (AU)')
        xlim([ElapsedTime(start_frame),ElapsedTime(final_frame)])
        ylim([0 9E5])
        title('Input: Protein Concentration')
        hold off

        subplot(1,2,2)
        plot(x,spot_temp,'-o',time_point_final,spot_interp,'-x');
        hold on
        plot(time_point_final,spot_interp_avr,'-','LineWidth',2)
        xlabel('Time (min)')
        ylabel('Fluorescence (AU)')
        xlim([ElapsedTime(start_frame),ElapsedTime(final_frame)])
        ylim([0 5E5])
        title('Output: Spot FLuorescence')
        legend('Spot Fluorescence','Integrated Spot Fluorescence')
        hold off

        saveas(TraceFig,[FigPath,'/traces_all/Particle_',num2str(i),'.jpg'])
    end
    %}
end
%% Plot traces for nuclei in different AP bin

test_frame = 70;
nuclei_plot = [];
APmin = APbinID(min(APbinRange));
APmax = APbinID(max(APbinRange));
X = [];

% Get nuclei within AP range
for i = 1:length(sch_pass)
    if (AP_pos(i,test_frame)<APmax) && (AP_pos(i,test_frame)>APmin)
        nuclei_plot = [nuclei_plot i];
        X = [X AP_pos(i,test_frame)];
    end
end
%% 
% Assign each nuclei to corresponding AP bin

edges = APbinID(APbinRange);
Y = discretize(X,edges);
%% 
% test

for i = 1:length(nuclei_plot)
    bin_plot = APbinRange(Y(i));
    spot_temp = spot_fluo_traces(nuclei_plot(i),:);
    nuclei_temp = nuclei_fluo_traces(nuclei_plot(i),:);

    if sum(spot_temp) > 0
        %{
        fig = figure;
        subplot(1,2,1)
        plot(ElapsedTime(1:final_frame),nuclei_temp,'-','LineWidth',2);
        xlabel('Time (min)')
        ylabel('Fluorescence (AU)')
        xlim([ElapsedTime(start_frame),ElapsedTime(final_frame)])
        ylim([0 1.1E6])
        %ylim([0 9E5])
        title('Input: Protein Concentration')
        hold off

        subplot(1,2,2)
        plot(ElapsedTime(1:final_frame),spot_temp,'- .','MarkerSize',12,'LineWidth',2)
        xlabel('Time (min)')
        ylabel('Fluorescence (AU)')
        xlim([ElapsedTime(start_frame),ElapsedTime(final_frame)])
        ylim([0 5E5])
        title('Output: Spot FLuorescence')
        legend('Spot Fluorescence')
        hold off
        
        sgtitle(['x = ',num2str(APbinID(bin_plot))]);

        FigAPPath = [FigPath,'/traces_AP/bin_',num2str(bin_plot),'_x_',num2str(APbinID(bin_plot))];
        mkdir(FigAPPath);
        %}
        
        fig = figure;
        yyaxis left
        plot(ElapsedTime(1:final_frame),nuclei_temp,'-','LineWidth',2);
        xlim([ElapsedTime(start_frame),ElapsedTime(final_frame)])
        %ylim([0 1.3E6])
        ylim([0 9E5])
        hold on

        yyaxis right
        plot(ElapsedTime(1:final_frame),spot_temp,'- .','MarkerSize',12,'LineWidth',2)
        xlabel('Time (min)')
        ylabel('Fluorescence (AU)')
        xlim([ElapsedTime(start_frame),ElapsedTime(final_frame)])
        ylim([0 3.5E5])
        hold off
        

        xlabel('Time (min)')
        ylabel('Fluorescence (AU)')
        legend('Input: Protein Concentration','Output: Spot FLuorescence')
        

        sgtitle(['x = ',num2str(APbinID(bin_plot))]);

        FigAPPath = [FigPath,'/traces_AP/bin_',num2str(bin_plot),'_x_',num2str(APbinID(bin_plot))];
        mkdir(FigAPPath);
        
        saveas(fig,[FigAPPath,'/Particle_',num2str(i),'.jpg'])
    end

end