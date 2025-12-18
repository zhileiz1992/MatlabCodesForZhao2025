% perform statistical tests
clear; close all;


%% General folder setting
fd_z4 = '/mnt/z4';
fd_base = fullfile(fd_z4, 'zz367', 'EphysMONAO', 'Analyzed');
birdIDs = {'pair5RigCCU29', 'pair4RigACU68', 'pair4RigBCU53', 'pair2RigBCU25'};
pairIDs = {'pair5CU29CU55', 'pair4CU68CU53', 'pair4CU68CU53', 'pair2CU20CU25'};
pretty_ids = {'M1', 'M2', 'M3', 'M4'};

bi=1;
birdID = birdIDs{bi};
pairID = pairIDs{bi};

% neuronID = '20240829-ch9';
neuronID = '20240912-ch11';


% load XGB result
fd_res = fullfile(fd_base, 'Figures', pairID, 'CompoundSyl', birdID, 'ensembleModel', 'XGB_res', neuronID);
fn_res = fullfile(fd_res, sprintf('%s.%s.metrics_all.mat', birdID, neuronID));
load(fn_res);



%% perform tests on compound syllables
y1 = metrics_all.Compound_Full.improvement_over_null;
y2 = metrics_all.Compound_Acoustic.improvement_over_null;
y3 = metrics_all.Compound_Time.improvement_over_null;
p12 = signrank(y1, y2);
p13 = signrank(y1, y3);
% compare the difference
p23 = signrank(y1-y2, y1-y3);


%% perform tests on calls
y1 = metrics_all.Call_Full.improvement_over_null;
y2 = metrics_all.Call_Acoustic.improvement_over_null;
y3 = metrics_all.Call_Time.improvement_over_null;
p12 = signrank(y1, y2);
p13 = signrank(y1, y3);
% compare the difference
p23 = signrank(y1-y2, y1-y3);




