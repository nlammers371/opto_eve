clear 
close all

load("S:\Jake\Dropbox\ProcessedData\eve_opto_v1\hmm_input_output_w7_K3_f2D_p2D.mat");
load("S:\Jake\Dropbox\ProcessedData\eve_opto_v1\nucleus_struct.mat");

%%
particle_id_index = unique([hmm_input_output.ParticleID]);
nc_particle_id_index = [nucleus_struct.ParticleID];
n_unique = numel(particle_id_index);
hmm_input_output = hmm_input_output(1:n_unique);

% add average protein levels back into hmm structure
for i = 1:numel(hmm_input_output)
    ParticleID = hmm_input_output(i).ParticleID;
    nc_ind = nc_particle_id_index == ParticleID;
    raw_protein = nucleus_struct(nc_ind).raw_nc_protein;
    % interpolate raw protein and add to structure
    raw_time = nucleus_struct(nc_ind).time;
    interp_time = hmm_input_output(i).time;
    hmm_input_output(i).nuclear_protein = interp1(raw_time,raw_protein,interp_time);
    hmm_input_output(i).source_path = nucleus_struct(nc_ind).source_path;
    hmm_input_output(i).SetID = floor(hmm_input_output(i).ParticleID);
    hmm_input_output(i).z_viterbi = hmm_input_output(i).r_inf(hmm_input_output(i).z_viterbi);
end
%%
keep_fields = {'time','ParticleID','SetID','source_path','fluo','r_vec','z_viterbi','f_viterbi','spot_protein','nuclear_protein','apPosParticle'};
rename_fields = {'Time','ParticleID','SetID','SourcePath','SpotFluorescence','InitiationRate','PromoterState','FluoViterbi','LocalKnirps','NuclearKnirps','APPosition'};
hmm_out = struct;
for i = 1:numel(hmm_input_output)
    for f = 1:numel(rename_fields)
        hmm_out(i).(rename_fields{f}) = hmm_input_output(i).(keep_fields{f});
    end
end
hmm_input_output = hmm_out;
    
    
save('S:\Jake\Dropbox\ProcessedData\eve_opto_v1\hmm_input_output.mat','hmm_input_output')
    
