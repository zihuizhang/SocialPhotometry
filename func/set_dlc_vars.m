function [dlc_vars,chamber_names] = set_dlc_vars(behav_type)
% set dlc variables depending on the session type
% will be used to parse csv columns
% need to be exactly the same as DLC config file
dlc_vars = {};
switch behav_type
    case 'threechamber'
        dlc_vars = {'snout', 'leftear','rightear','tailbase','implant',...
            'cupA','cupB','upperleftA','lowerleftA','upperrightA','lowerrightA',...
            'upperleftB','lowerleftB','upperrightB','lowerrightB'};
    case 'openfield'
        dlc_vars = {'snout','leftear','rightear','tailbase','implant',...
            'upperleft','lowerleft','upperright','lowerright'};
    case 'openfieldTwoMice'
        dlc_vars = {'snout','leftear','rightear','tailbase','implant',...
            'juv_snout','juv_leftear','juv_rightear','juv_tailbase',...
            'upperleft','lowerleft','upperright','lowerright'};
    case 'homecageSingleMouse'
        dlc_vars = {'snout','leftear','rightear','tailbase','implant',...
            'upperleft','lowerleft','upperright','lowerright'};
    case 'homecageFrootloop'
        dlc_vars = {'snout','leftear','rightear','tailbase','implant','frootloop'};
    case 'openfieldFrootloop'
        dlc_vars = {'snout','leftear','rightear','tailbase','implant','frootloop'};
    case {'HomecageIntruderBox61','JuvenileInteractionBox61','HomecageTwoMice(Inscopics)'}
        dlc_vars = { 'Ear_left_1', 'Ear_right_1','Nose_1', 'Center_1', 'Lateral_left_1', 'Lateral_right_1', 'Tail_base_1', 'Tail_end_1',... % resident
            'Ear_left_2', 'Ear_right_2','Nose_2', 'Center_2', 'Lateral_left_2', 'Lateral_right_2', 'Tail_base_2', 'Tail_end_2'}; % intruder
    case 'zoom150_homecagejuvenile'
        dlc_vars = { 'mouse_nose1', 'mouse_body1','mouse_butt1', 'juv_nose1', 'juv_body1', 'juv_butt1',...
            'mouse_nose2', 'mouse_body2','mouse_butt2', 'juv_nose2', 'juv_body2', 'juv_butt2'}; % two boxes
    case 'Juvenile'
        dlc_vars = { 'mouse_nose1', 'mouse_body1','mouse_butt1', 'juv_nose1', 'juv_body1', 'juv_butt1',...
            'mouse_nose2', 'mouse_body2','mouse_butt2', 'juv_nose2', 'juv_body2', 'juv_butt2'}; % two boxes
    case 'zoom150_homecagetoymouse'
        dlc_vars = { 'mouse_nose1', 'mouse_body1','mouse_butt1', 'juv_nose1', 'juv_body1', 'juv_butt1',...
            'mouse_nose2', 'mouse_body2','mouse_butt2', 'juv_nose2', 'juv_body2', 'juv_butt2'}; % two boxes
    case 'zoom150_homecagefrootloop'
        dlc_vars = {'mouse_nose1', 'mouse_body1','mouse_butt1', 'frootloop1',...
            'mouse_nose2', 'mouse_body2','mouse_butt2', 'frootloop2'}; % two boxes
    case 'zoom150_threecham'
        dlc_vars = { 'mouse_nose1', 'mouse_body1','mouse_butt1', 'box1_side1_a', 'box1_side1_b', 'box1_side1_c','box1_side1_d',...
            'box1_side2_a', 'box1_side2_b', 'box1_side2_c','box1_side2_d','box1_cup1','box1_cup2',...
            'mouse_nose2', 'mouse_body2','mouse_butt2', 'box2_side1_a', 'box2_side1_b', 'box2_side1_c','box2_side1_d',...
            'box2_side2_a', 'box2_side2_b', 'box2_side2_c','box2_side2_d','box2_cup1','box2_cup2'}; % two boxes
    case 'STAAR_SingleMouseOnly'
        dlc_vars = {'Nose','HeadCenter','Tail_base'};
    case 'SIMBA_Aggression'
        dlc_vars = { 'Ear_left_1', 'Ear_right_1','Nose_1', 'Center_1', 'Lat_left_1', 'Lat_right_1', 'Tail_base_1', 'Tail_end_1',... % resident
            'Ear_left_2', 'Ear_right_2','Nose_2', 'Center_2', 'Lat_left_2', 'Lat_right_2', 'Tail_base_2', 'Tail_end_2'}; % intruder

end

% for threechamber
switch behav_type
    case 'threechamber'
        chamber_names = {'A','B'};
    case 'zoom150_threecham'
        chamber_names = {'side1','side2'};
    otherwise
        chamber_names = {''};
end


end