clear
close all

load('S:\Jake\Dropbox\ProcessedData\eve_opto_v1\hmm_input_output.mat','hmm_input_output')
FigPath = 'S:\Jake\Dropbox\LocalEnrichmentFigures\PipelineOutput\eve_opto_v1\hmm_plots\';

%%
index_vec = 1:numel(hmm_input_output);
n_samples = 100;
rng(123)
plot_indices = randsample(index_vec,n_samples,false);

for p = plot_indices(10)
    ParticleID = hmm_input_output(p).ParticleID;
%     r = hmm_input_output(p).r_inf;
    fluo_fig = figure;
    hold on
    plot(hmm_input_output(p).Time,hmm_input_output(p).SpotFluorescence,'LineWidth',1.5)
    plot(hmm_input_output(p).Time,hmm_input_output(p).FluoViterbi,'--','LineWidth',1.5)
    
    ylabel('spot fluorescence (au)')
    yyaxis right
    
    plot(hmm_input_output(p).Time,hmm_input_output(p).NuclearKnirps,'-','Color','g','LineWidth',1.5)
    ylabel('knirps level (au)')
    
    xlabel('time')
    
    set(gca,'Fontsize',14)
    
    saveas(fluo_fig,[FigPath 'particle' num2str(ParticleID) '_fluo_plot.png'])
    
        
    v_fig = figure;
    hold on
    plot(hmm_input_output(p).Time,hmm_input_output(p).InitiationRate,'LineWidth',1.5)
    stairs(hmm_input_output(p).Time,hmm_input_output(p).PromoterState,'-','LineWidth',1.5)
    
    ylabel('initiation rate')
    yyaxis right
    
    plot(hmm_input_output(p).Time,hmm_input_output(p).NuclearKnirps,'-','Color','g','LineWidth',1.5)
    ylabel('knirps level (au)')
    
    xlabel('time')
    
    set(gca,'Fontsize',14)
    
    saveas(v_fig,[FigPath 'particle' num2str(ParticleID) '_v_plot.png'])
end    
    