% This code is intended for testing the boundary extraction process

clear all
close all
clc

%% Part 1: Read and display data
load('./data/2020-02-23-optoknirps_new_embryo6_lin.mat');
load('./data/CompiledParticles.mat')

start_frame = 40;
final_frame = 149;

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
        fluo = max(schnitzcells(sch_now).Fluo(j,:));
        try
            processed_data(frame_now).xcoord = [processed_data(frame_now).xcoord, x_coord];
            processed_data(frame_now).ycoord = [processed_data(frame_now).ycoord, y_coord];
            processed_data(frame_now).schnitznum = [processed_data(frame_now).schnitznum, sch_now];
            processed_data(frame_now).SpotFluo = [processed_data(frame_now).SpotFluo, 0];
            processed_data(frame_now).NuclearFluo = [processed_data(frame_now).NuclearFluo fluo];  
        catch
            processed_data(frame_now).xcoord = x_coord;
            processed_data(frame_now).ycoord = y_coord;
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

%{
frame_plot = 149;
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
%}

% Convert the format: Store them into individual nuclei traces
nuclei_fluo_traces = zeros(length(sch_pass),final_frame);
spot_fluo_traces = zeros(length(sch_pass),final_frame);
x_pos = zeros(length(sch_pass),final_frame);
y_pos = zeros(length(sch_pass),final_frame);

for i = 1:length(sch_pass)
    for j = start_frame:final_frame
        nuclei_fluo_traces(i,j) = processed_data(j).NuclearFluo(i);
        spot_fluo_traces(i,j) = processed_data(j).SpotFluo(i);
        x_pos(i,j) = processed_data(j).xcoord(i);
        y_pos(i,j) = processed_data(j).ycoord(i);
    end
end

%% Part 4: Process traces
% We would like to combine everything together and interpolate values


% Define time point to interpolate
time_point_final = 12:1/6:48;
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

    %spot_temp = spot_fluo_output(i,:);
    if sum(spot_temp) > 0
        %{
        fig = figure(5);
        subplot(1,2,1)
        plot(x,nuclei_temp,'-o',time_point_final,nuclei_interp,'-x');
        hold on
        plot(time_point_final,nuclei_interp_avr,'-','LineWidth',2)
        xlabel('Time (min)')
        ylabel('Fluorescence (AU)')
        xlim([12 48])
        ylim([0 9E5])
        title('Input: Protein Concentration')
        hold off

        subplot(1,2,2)
        plot(x,spot_temp,'-o',time_point_final,spot_interp,'-x');
        hold on
        plot(time_point_final,spot_interp_avr,'-','LineWidth',2)
        xlabel('Time (min)')
        ylabel('Fluorescence (AU)')
        xlim([12 48])
        ylim([0 5E5])
        title('Output: Spot FLuorescence')
        legend('Spot Fluorescence','Integrated Spot Fluorescence')
        hold off

        saveas(fig,['./figure/traces/Particle_',num2str(i),'.jpg'])
        %}
    end

end


%% Part 5: Visualize input and output
%{
factor = 1E4;

%for i = 140
for i = start_frame:final_frame
    % Plot voronoi first
    xpos = processed_data(i).xcoord;
    ypos = processed_data(i).ycoord;

    % try to plot and test a bit...
    pts = [xpos' ypos'];

    fig = figure(2);
    clf
    [v,c] = voronoin(double(pts));

    for j = 1:length(c)
        if all(c{j}~=1)
        x = v(c{j},1);
        y = v(c{j},2);
        a = processed_data(i).NuclearFluo(j);
        patch(x,y,a);
        colorbar
        caxis([1E5 8E5])
        colormap('jet');
        end
    end

    axis equal
    xlim([0 1024])
    ylim([0 256])
    hold on
    
    % Plot Spots
    for j = 1:length(processed_data(i).xcoord)
        x = processed_data(i).xcoord(j);
        y = processed_data(i).ycoord(j);
        marker_size = processed_data(i).SpotFluo(j)/factor;
        if marker_size>0
            plot(x, y, '.', 'MarkerSize',marker_size,'Color','#A2142F')
        end
    end
    hold off
    title(['Time: ',num2str(ElapsedTime(i)),' min'])
    
    saveas(fig,['./figure/test/Frame_',num2str(i),'.jpg'])
    
end
%}

%% Part 6: Plot input-output function for each stripe
% initialize some parameters
test_frame = 149;

nuclei_stripe4 = [];
nuclei_stripe6 = [];

edges = (1.5:1:7.5)*1E5;
len = length(edges)-1;
xmid = 0.5*(edges(1:end-1)+edges(2:end));

% Store calculated results
y4_result = zeros(length(sch_pass),length(xmid));
y6_result = zeros(length(sch_pass),length(xmid));
err4_result = zeros(length(sch_pass),length(xmid));
err6_result = zeros(length(sch_pass),length(xmid));

% Get nuclei for stripe 4
for i = 1:length(sch_pass)
    if (x_pos(i,test_frame)<475) && (x_pos(i,test_frame)>350)
        nuclei_stripe4 = [nuclei_stripe4 i];
    end
end

% Get nuclei for stripe 6
for i = 1:length(sch_pass)
    if (x_pos(i,test_frame)<740) && (x_pos(i,test_frame)>625)
        nuclei_stripe6 = [nuclei_stripe6 i];
    end
end

for i = 1:length(time_point_final)
    
    x4 = nuclei_fluo_avr(nuclei_stripe4,i);
    y4 = spot_fluo_avr(nuclei_stripe4,i);
    
    x6 = nuclei_fluo_avr(nuclei_stripe6,i);
    y6 = spot_fluo_avr(nuclei_stripe6,i);
    
    [~,~,loc4]=histcounts(x4,edges);
    yplot4 = accumarray(loc4(loc4>0),y4(loc4>0),[len 1])./accumarray(loc4(loc4>0),1,[len 1]);
    err4 = sqrt(accumarray(loc4(loc4>0),(y4(loc4>0)).^2,[len 1])./accumarray(loc4(loc4>0),1,[len 1]))./sqrt(accumarray(loc4(loc4>0),1,[len 1]));
    
    [~,~,loc6]=histcounts(x6,edges);
    yplot6 = accumarray(loc6(loc6>0),y6(loc6>0),[len 1])./accumarray(loc6(loc6>0),1,[len 1]);
    err6 = sqrt(accumarray(loc6(loc6>0),(y6(loc6>0)).^2,[len 1])./accumarray(loc6(loc6>0),1,[len 1]))./sqrt(accumarray(loc6(loc6>0),1,[len 1]));
    
    y4_result(i,:) = yplot4;
    y6_result(i,:) = yplot6;
    err4_result(i,:) = err4;
    err6_result(i,:) = err6;
    
    %{
    fig = figure(5);
    plot(x4,y4,'.','MarkerSize',15)
    hold on
    plot(x6,y6,'.','MarkerSize',15)
    boundedline(xmid, yplot4, err4, '-b*','nan', 'gap','alpha');
    boundedline(xmid, yplot6, err6, '-r*','nan', 'gap','alpha');
    %plot(xmid,y_plot,'-o','LineWidth',2)
    xlim([1.5E5 7E5])
    ylim([0 3.5E5])
    title(['Time: ',num2str(time_point_final(i)),' min'])
    xlabel('Knirps Protein Input (AU)')
    ylabel('mRNA Production Rate (AU)')
    legend('stripe 4', 'stripe 6','FontSize',16)
    saveas(fig,['./figure/input_output_stripe_4/Frame_',num2str(i),'.jpg'])
    hold off
    %}
end

%% Part 7: Plot input-output function for different time points

time_plot = [97,157,217];

cmap_temp = [[0 0.4470 0.7410];[0.8500 0.3250 0.0980];[0.9290 0.6940 0.1250];[0.4940 0.1840 0.5560]];

% Stripe 4
fig = figure(6);
%errorbar(xmid, y4_result(time_plot(1),:), err4_result(time_plot(1),:),'-*','LineWidth',2);
%hold on
%errorbar(xmid, y4_result(time_plot(2),:), err4_result(time_plot(2),:), '-*','LineWidth',2);
%errorbar(xmid, y4_result(time_plot(3),:), err4_result(time_plot(3),:), '-*','LineWidth',2);
%errorbar(xmid, y4_result(time_plot(4),:), err4_result(time_plot(4),:), '-*','LineWidth',2);

boundedline(xmid, y4_result(time_plot(1),:), err4_result(time_plot(1),:),'-*','nan', 'gap','alpha','cmap', cmap_temp(1,:));
hold on
boundedline(xmid, y4_result(time_plot(2),:), err4_result(time_plot(2),:), '-*','nan', 'gap','alpha','cmap', cmap_temp(2,:));
boundedline(xmid, y4_result(time_plot(3),:), err4_result(time_plot(3),:), '-*','nan', 'gap','alpha','cmap', cmap_temp(3,:));
%boundedline(xmid, y4_result(time_plot(4),:), err4_result(time_plot(4),:), '-r*','nan', 'gap','alpha','cmap', cmap_temp(4,:));

xlim([1.5E5 7E5])
ylim([0 3.5E5])

xlabel('Knirps Protein Input (AU)')
ylabel('mRNA Production Rate (AU)')
title('Input-Output Function (Stripe 4)')

legend('28 min','error','38 min','error','48 min','error','FontSize',16)
%legend('28 min','38 min','48 min','FontSize',16)
saveas(fig,'./figure/input_output_time_stripe_4.pdf')

% Stripe 6
fig = figure(7);
%errorbar(xmid, y6_result(time_plot(1),:), err6_result(time_plot(1),:),'-*','LineWidth',2);
%hold on
%errorbar(xmid, y6_result(time_plot(2),:), err6_result(time_plot(2),:), '-*','LineWidth',2);
%errorbar(xmid, y6_result(time_plot(3),:), err6_result(time_plot(3),:), '-*','LineWidth',2);
%errorbar(xmid, y6_result(time_plot(4),:), err6_result(time_plot(4),:), '-*','LineWidth',2);

boundedline(xmid, y6_result(time_plot(1),:), err6_result(time_plot(1),:),'-*','nan', 'gap','alpha','cmap', cmap_temp(1,:));
hold on
boundedline(xmid, y6_result(time_plot(2),:), err6_result(time_plot(2),:), '-*','nan', 'gap','alpha','cmap', cmap_temp(2,:));
boundedline(xmid, y6_result(time_plot(3),:), err6_result(time_plot(3),:), '- *','nan', 'gap','alpha','cmap', cmap_temp(3,:));
%boundedline(xmid, y6_result(time_plot(4),:), err6_result(time_plot(4),:), '-r*','nan', 'gap','alpha','cmap', cmap_temp(4,:));

xlim([1.5E5 7E5])
ylim([0 3.5E5])

xlabel('Knirps Protein Input (AU)')
ylabel('mRNA Production Rate (AU)')
title('Input-Output Function (Stripe 6)')

%legend('28 min','38 min','48 min','FontSize',16)
legend('28 min','error','38 min','error','48 min','error','FontSize',16)
saveas(fig,'./figure/input_output_time_stripe_6.pdf')

%% Plot traces for stripe 4 and stripe 6

%{
% For stripe 4
for i = 1:length(nuclei_stripe4)
    
    spot_temp = spot_fluo_traces(nuclei_stripe4(i),:);
    nuclei_temp = nuclei_fluo_traces(nuclei_stripe4(i),:);

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

    %spot_temp = spot_fluo_output(i,:);
    if sum(spot_temp) > 0
        
        fig = figure(5);
        subplot(1,2,1)
        plot(x,nuclei_temp,'-o',time_point_final,nuclei_interp,'-x');
        hold on
        plot(time_point_final,nuclei_interp_avr,'-','LineWidth',2)
        xlabel('Time (min)')
        ylabel('Fluorescence (AU)')
        xlim([12 48])
        ylim([0 9E5])
        title('Input: Protein Concentration')
        hold off

        subplot(1,2,2)
        plot(x,spot_temp,'-o',time_point_final,spot_interp,'-x');
        hold on
        plot(time_point_final,spot_interp_avr,'-','LineWidth',2)
        xlabel('Time (min)')
        ylabel('Fluorescence (AU)')
        xlim([12 48])
        ylim([0 5E5])
        title('Output: Spot FLuorescence')
        legend('Spot Fluorescence')
        hold off

        saveas(fig,['./figure/stripe4_traces/Particle_',num2str(i),'.jpg'])
    end

end
%}

% For stripe 6
for i = 1:length(nuclei_stripe6)
    
    spot_temp = spot_fluo_traces(nuclei_stripe6(i),:);
    nuclei_temp = nuclei_fluo_traces(nuclei_stripe6(i),:);

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

    %spot_temp = spot_fluo_output(i,:);
    if sum(spot_temp) > 0
        
        fig = figure(5);
        subplot(1,2,1)
        plot(x,nuclei_temp,'-o',time_point_final,nuclei_interp,'-x');
        hold on
        plot(time_point_final,nuclei_interp_avr,'-','LineWidth',2)
        xlabel('Time (min)')
        ylabel('Fluorescence (AU)')
        xlim([12 48])
        ylim([0 9E5])
        title('Input: Protein Concentration')
        hold off

        subplot(1,2,2)
        plot(x,spot_temp,'-o',time_point_final,spot_interp,'-x');
        hold on
        plot(time_point_final,spot_interp_avr,'-','LineWidth',2)
        xlabel('Time (min)')
        ylabel('Fluorescence (AU)')
        xlim([12 48])
        ylim([0 5E5])
        title('Output: Spot FLuorescence')
        legend('Spot Fluorescence')
        hold off

        saveas(fig,['./figure/stripe6_traces/Particle_',num2str(i),'.jpg'])
    end

end

