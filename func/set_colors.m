function [color] = set_colors()
% define colors for photometry analysis
color = struct();
color.red = [217 95 95]./255;
color.green = [114 204 123]./255;
color.uv = [.6 .6 .6];
color.DAT = [249, 109, 21]./255; % DA
color.SERT = [80, 182, 187]./255; %5HT

color.Nosepoke_R = [250, 200, 0]./255;
color.Nosepoke_L = [107 188 255]./255;
color.Nosepoke = [191, 77, 255]./255;
color.Magazine = [.7 .7 .7];
color.RewardPortEntry = [.5 .5 .5];
color.Reward = [0 0 0];
color.Shock = [1 0 0];
color.Photostim = [255 140 133]./255;
color.StimsTimesOnset = color.Photostim;
color.StimsTimesOffset = tint(color.Photostim,.5);

color.RewardCue = [.1 .1 .1];
color.RewardNoShock = color.Reward;

color.normal_red = color.red;
color.normal_green = color.green;
color.normal_uv = color.uv;

color.frust_red = color.red.*0.5;
color.frust_green = color.green.*0.5;
color.frust_uv = color.uv.*0.5;

% add drug color here
color.Saline = [.5 .5 .5];
color.MDMA = [116, 63, 176]./255;
color.Vehicle = [.5 .5 .5];
color.Meth = [255 176 82]./255;
color.Saline = [.5 .5 .5];
color.MDMA = [116, 63, 176]./255;
color.Vehicle = [.5 .5 .5];
color.Meth = [255 176 82]./255;
color.CP93129 = [50, 164, 251]./255; % 1b agonist
color.LY344864 = [50, 64, 251]./255; % 1f agonist
color.WAY163909 = [251, 50, 164]./255; % 2c agonist;

color.Attack = shade([1,0,0],.5);
color.AttackOnset = shade([1,0,0],.5);
color.AttackOffset = tint(color.AttackOnset,.5);
color.AggressionOnset = [0 0 0];
color.AggressionOffset = tint(color.AggressionOnset,.5);
color.SocialEntry = [0 0 0];
color.SocialOnset = [.3 .3 .3];
color.SocialOffset = tint(color.SocialOnset,.5);
color.SocialNoAttackOnset = [115, 137, 250]./255;
color.SocialNoAttackOffset = tint(color.SocialNoAttackOnset,.5);

inhib = [55 136 254]./255;
color.InhibConsumption = inhib;
color.InhibDelivery = inhib;
color.InhibNosePoke = inhib;
color.PreInhibConsumption = [.5 .5 .5];
color.PreInhibDelivery = [.5 .5 .5];
color.PreInhibNosePoke = [.5 .5 .5];
color.PostInhibConsumption = [0 0 0];
color.PostInhibDelivery = [0 0 0];
color.PostInhibNosePoke = [0 0 0];

stim = [253 72 72]./255;
color.StimConsumption = stim;
color.StimDelivery = stim;
color.StimNosePoke = stim;
color.PreStimConsumption = [.5 .5 .5];
color.PreStimDelivery = [.5 .5 .5];
color.PreStimNosePoke = [.5 .5 .5];
color.PostStimConsumption = [0 0 0];
color.PostStimDelivery = [0 0 0];
color.PostStimNosePoke = [0 0 0];


end

