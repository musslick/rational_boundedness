%%
%i think location key responses don't match up with actual locations
%%
clear all
close all
remap = NaN(1, 52);
remap(49) = '1';
remap(50) = '2';
remap(51) = '3';
remap(52) = '4';

sub30 = 1:30;
badbehavior = [101 103 117];

badWP =  [1     3     9    14    16    17    18    23    25    26    28];
subject = [102 104:108 110:113 115 119:122 124 127 129 130]

%badbehavior = [101 103 108:111 114 117 119 124];
nansubs = badbehavior- 100;
subjectz = [2     4     5     6     7    12    13    15    16    18    20    21    22    23    25    26    27    28    29    30] + 100;
goodsubjects = subjectz -100;

subjectz = [102 104:116 118:130]
goodsubjects = subjectz -100; %[102 104:108 110:113 115 119:122 124 127 129 130] -100;

for subject = subjectz
    load('subjectparameters.mat')
    keyboardMappings = SubjectParameters.keyboardMappings(subject, :);
    
    load(['SpockData/' num2str(subject) ...
        '/SingleTaskColorNaming_Subject_' num2str(subject) '_VOICE_Block1.mat']);
    disp(['SpockData/' num2str(subject) ...
        '/SingleTaskColorNaming_Subject_' num2str(subject) '_VOICE_Block1.mat']);
    
    RT_data{subject}.CN = updatedTable.BlockResponseVoiceTime(updatedTable.taskAcc == 1);
    RT_data_incorrect{subject}.CN = updatedTable.BlockResponseVoiceTime(updatedTable.taskAcc == 0);
    
    RTCNCongruent(subject-100) = nanmean(updatedTable.BlockResponseVoiceTime(updatedTable.word == updatedTable.color));
    RTCNInCongruent(subject-100) = nanmean(updatedTable.BlockResponseVoiceTime(updatedTable.word ~= updatedTable.color));
    CongVoices = (updatedTable.BlockResponseVoiceTime(updatedTable.word == updatedTable.color));
    CongCNThresh(subject-100) = nanmean(CongVoices) + 2*nanstd(CongVoices);
    CongVoiceAcc(subject-100) = nanmean(updatedTable.taskAcc(updatedTable.word == updatedTable.color)...
        .*(updatedTable.BlockResponseVoiceTime(updatedTable.word == updatedTable.color)<CongCNThresh(subject-100)));
    
    InCongVoices = (updatedTable.BlockResponseVoiceTime(updatedTable.word ~= updatedTable.color));
    InCongCNThresh(subject-100) = nanmean(InCongVoices) + 2*nanstd(InCongVoices);
    InCongVoiceAcc(subject-100) = nanmean(updatedTable.taskAcc(updatedTable.word ~= updatedTable.color)...
        .*(updatedTable.BlockResponseVoiceTime(updatedTable.word ~= updatedTable.color)< InCongCNThresh(subject-100)));
    
    RawCongCNAcc(subject-100) = nanmean(updatedTable.taskAcc(updatedTable.word == updatedTable.color));
    RawInCongCNAcc(subject-100) = nanmean(updatedTable.taskAcc(updatedTable.word ~= updatedTable.color));
    
    
    ColorNamingWordErrors(subject-100) = sum(updatedTable.taskAcc == 0 & ...
        updatedTable.voiceresponse == updatedTable.word & ...
        updatedTable.word ~= updatedTable.color) ...
        /sum(updatedTable.taskAcc == 0 & updatedTable.word ~= updatedTable.color & updatedTable.voiceresponse > 0);
    
    
    matFiles = dir(['SessionData/' num2str(subject) '/SingleTaskWordPointing_*.mat']);
    numfiles = length(matFiles);
    d = matFiles;
    [dx,dx] = sort([d.datenum],'ascend');
    newest = d(dx(1)).name;
    load(['SessionData/' num2str(subject) '/' newest]);
    
    [keyresponse, wordAcc, locAcc ] =   covertkeyboardpresses(subject,'',CSVData);
    
    RT_data{subject}.WP = CSVData.BlockResponseKeyTime(CSVData.wordAcc == 1);
    RT_data_incorrect{subject}.WP = CSVData.BlockResponseKeyTime(CSVData.wordAcc == 0);
    
    RTWPCongruent(subject-100) = nanmean(CSVData.BlockResponseKeyTime(CSVData.word == CSVData.color));
    RTWPInCongruent(subject-100) = nanmean(CSVData.BlockResponseKeyTime(CSVData.word ~= CSVData.color));
    CongWPThresh(subject-100) = nanmean(CSVData.BlockResponseKeyTime(CSVData.word == CSVData.color)) ...
        + 2*nanstd(CSVData.BlockResponseKeyTime(CSVData.word == CSVData.color));
    InCongWPThresh(subject-100) = nanmean(CSVData.BlockResponseKeyTime(CSVData.word ~= CSVData.color)) ...
        + 2*nanstd(CSVData.BlockResponseKeyTime(CSVData.word ~= CSVData.color));
    RawCongWPAcc(subject-100) = nanmean(CSVData.wordAcc(CSVData.word == CSVData.color));
    RawInCongWPAcc(subject-100) = nanmean(CSVData.wordAcc(CSVData.word ~= CSVData.color));
    
    for i = 1:length(keyresponse)
        ColorErrors(i)= wordAcc(i) == 0 & str2num(SubjectParameters.moves{keyboardMappings(CSVData.color(i))}) == (keyresponse(i)) ....
            & CSVData.word(i) ~= CSVData.color(i);
    end
    WordPointingColorErrors(subject-100) = sum(ColorErrors)/sum(CSVData.wordAcc == 0 & CSVData.word ~= CSVData.color & keyresponse > 0);
    
    matFiles = dir(['SessionData/' num2str(subject) '/SingleTaskLocationPointing_*.mat']);
    numfiles = length(matFiles);
    d = matFiles;
    [dx,dx] = sort([d.datenum],'ascend');
    newest = d(dx(1)).name;
    load(['SessionData/' num2str(subject) '/' newest]);
    [keyresponse, wordAcc, locAcc ] =   covertkeyboardpresses(subject,'LP',CSVData);
    
    RT_data{subject}.LP = CSVData.BlockResponseKeyTime(CSVData.taskAcc == 1);
    RT_data_incorrect{subject}.LP = CSVData.BlockResponseKeyTime(CSVData.taskAcc == 0);
    
    RTLPCongruent(subject-100) = nanmean(CSVData.BlockResponseKeyTime(CSVData.word == CSVData.color));
    RTLPInCongruent(subject-100) = nanmean(CSVData.BlockResponseKeyTime(CSVData.word ~= CSVData.color));
    CongLPThresh(subject-100) = nanmean(CSVData.BlockResponseKeyTime(CSVData.word == CSVData.color)) ...
        + 2*nanstd(CSVData.BlockResponseKeyTime(CSVData.word == CSVData.color));
    InCongLPThresh(subject-100) = nanmean(CSVData.BlockResponseKeyTime(CSVData.word ~= CSVData.color)) ...
        + 2*nanstd(CSVData.BlockResponseKeyTime(CSVData.word ~= CSVData.color));
    RawCongLPAcc(subject-100) = nanmean(CSVData.taskAcc(CSVData.word == CSVData.color));
    RawCongLPAccStd(subject-100, :) = (CSVData.taskAcc(CSVData.word == CSVData.color));
    RawInCongLPAcc(subject-100) = nanmean(CSVData.taskAcc(CSVData.word ~= CSVData.color));
    
    
    for i = 1:length(CSVData.wordAcc)
        
        
        WordErrors(i)= locAcc(i) == 0 & str2num(SubjectParameters.moves{keyboardMappings(CSVData.word(i))}) == (keyresponse(i)) ...
            & CSVData.word(i) ~= CSVData.color(i);
    end
    LocationPointingWordErrors(subject-100) = sum(WordErrors)/sum(locAcc == 0 & CSVData.word ~= CSVData.color & keyresponse > 0);
    
    load(['SpockData/' num2str(subject) ...
        '/MultitaskingCNLP_Subject_' num2str(subject) '_VOICE.mat']);
    
    [keyresponse, wordAcc, locAcc ] =   covertkeyboardpresses(subject,'LP',updatedTable);
    
    totalAccuracy = updatedTable.locAcc .* updatedTable.colorAcc;
    voiceRT = updatedTable.BlockResponseVoiceTime(totalAccuracy == 1);
    pointRT = updatedTable.BlockResponseKeyTime(totalAccuracy == 1);
    RT_data{subject}.CNLP = max(voiceRT, pointRT);
    
    voiceRT = updatedTable.BlockResponseVoiceTime(totalAccuracy == 0);
    pointRT = updatedTable.BlockResponseKeyTime(totalAccuracy == 0);
    RT_data_incorrect{subject}.CNLP = max(voiceRT, pointRT);
    
    
    RTCNLPVoiceCongruent(subject-100) = nanmean(updatedTable.BlockResponseVoiceTime(updatedTable.word == updatedTable.color));
    RTCNLPVoiceInCongruent(subject-100) = nanmean(updatedTable.BlockResponseVoiceTime(updatedTable.word ~= updatedTable.color));
    AccCNLPVoiceCongruent(subject-100) = nanmean(updatedTable.colorAcc(updatedTable.word == updatedTable.color)...
        .*(updatedTable.BlockResponseVoiceTime(updatedTable.word == updatedTable.color)<CongCNThresh(subject-100)));
    AccCNLPVoiceInCongruent(subject-100) = nanmean(updatedTable.colorAcc(updatedTable.word ~= updatedTable.color)...
        .*(updatedTable.BlockResponseVoiceTime(updatedTable.word ~= updatedTable.color)<InCongCNThresh(subject-100)));
    AccCNLPPointCongruent(subject-100) = nanmean(updatedTable.locAcc(updatedTable.word == updatedTable.color)...
        .*(updatedTable.BlockResponseKeyTime(updatedTable.word == updatedTable.color)<CongLPThresh(subject-100)));
    AccCNLPPointInCongruent(subject-100) = nanmean(updatedTable.locAcc(updatedTable.word ~= updatedTable.color)...
        .*(updatedTable.BlockResponseKeyTime(updatedTable.word ~= updatedTable.color)<InCongLPThresh(subject-100)));
    
    AccCNLPCong(subject-100, :) = (updatedTable.colorAcc(updatedTable.word == updatedTable.color)...
        .*(updatedTable.BlockResponseVoiceTime(updatedTable.word == updatedTable.color)<CongCNThresh(subject-100))) ...
        .* (updatedTable.locAcc(updatedTable.word == updatedTable.color)...
        .*(updatedTable.BlockResponseKeyTime(updatedTable.word == updatedTable.color)<CongLPThresh(subject-100)));
    
    AccCNLPInCong(subject-100, :) = (updatedTable.colorAcc(updatedTable.word ~= updatedTable.color)...
        .*(updatedTable.BlockResponseVoiceTime(updatedTable.word ~= updatedTable.color)<InCongCNThresh(subject-100))) ...
        .* (updatedTable.locAcc(updatedTable.word ~= updatedTable.color)...
        .*(updatedTable.BlockResponseKeyTime(updatedTable.word ~= updatedTable.color)<InCongLPThresh(subject-100)));
    
    RawAccCNLPVoiceCongruent(subject-100) = nanmean(updatedTable.colorAcc(updatedTable.word == updatedTable.color));
    RawAccCNLPVoiceInCongruent(subject-100) = nanmean(updatedTable.colorAcc(updatedTable.word ~= updatedTable.color));
    RawAccCNLPPointCongruent(subject-100) = nanmean(updatedTable.locAcc(updatedTable.word == updatedTable.color));
    RawAccCNLPPointInCongruent(subject-100) = nanmean(updatedTable.locAcc(updatedTable.word ~= updatedTable.color));
    
    for i = 1:length(keyresponse)
        LPWordErrors(i)= locAcc(i) == 0 & ...
            str2num(SubjectParameters.moves{keyboardMappings(updatedTable.word(i))}) == keyresponse(i) & ...
            updatedTable.location(i) ~= updatedTable.word(i);
        LPColorErrors(i)= locAcc(i) == 0 & ...
            str2num(SubjectParameters.moves{keyboardMappings(updatedTable.color(i))}) == keyresponse(i) & ...
            updatedTable.color(i) ~= updatedTable.word(i);
    end
    
    CNLPWordErrors(subject-100) = sum(LPWordErrors)/sum(locAcc == 0  & updatedTable.color ~= updatedTable.word);
    CNLPColorErrors(subject-100) = sum(LPColorErrors)/sum(locAcc == 0  & updatedTable.color ~= updatedTable.word);
    
    
    
    LPColorNamingWordErrors(subject-100) = sum(updatedTable.colorAcc == 0 & updatedTable.voiceresponse == updatedTable.word) ...
        /sum(updatedTable.colorAcc == 0 & updatedTable.color  ~= updatedTable.word & updatedTable.voiceresponse > 0);
    
    
    load(['SpockData/' num2str(subject) ...
        '/MultitaskingCNWP_Subject_' num2str(subject) '_VOICE.mat']);
    
    totalAccuracy = updatedTable.wordAcc .* updatedTable.colorAcc;
    voiceRT = updatedTable.BlockResponseVoiceTime(totalAccuracy == 1);
    pointRT = updatedTable.BlockResponseKeyTime(totalAccuracy == 1);
    RT_data{subject}.CNWP = max(voiceRT, pointRT);
    
    voiceRT = updatedTable.BlockResponseVoiceTime(totalAccuracy == 0);
    pointRT = updatedTable.BlockResponseKeyTime(totalAccuracy == 0);
    RT_data_incorrect{subject}.CNWP = max(voiceRT, pointRT);
    
    [keyresponse, wordAcc, locAcc ] =   covertkeyboardpresses(subject,'',updatedTable);
    RTCNWPVoiceCongruent(subject-100) = nanmean(updatedTable.BlockResponseVoiceTime(updatedTable.word == updatedTable.color));
    RTCNWPVoiceInCongruent(subject-100) = nanmean(updatedTable.BlockResponseVoiceTime(updatedTable.word ~= updatedTable.color));
    AccCNWPVoiceCongruent(subject-100) = nanmean(updatedTable.colorAcc(updatedTable.word == updatedTable.color)...
        .*(updatedTable.BlockResponseVoiceTime(updatedTable.word == updatedTable.color)<CongCNThresh(subject-100)));
    AccCNWPVoiceInCongruent(subject-100) = nanmean(updatedTable.colorAcc(updatedTable.word ~= updatedTable.color)...
        .*(updatedTable.BlockResponseVoiceTime(updatedTable.word ~= updatedTable.color)<InCongCNThresh(subject-100)));
    AccCNWPPointCongruent(subject-100) = nanmean(updatedTable.wordAcc(updatedTable.word == updatedTable.color)...
        .*(updatedTable.BlockResponseKeyTime(updatedTable.word == updatedTable.color)<CongWPThresh(subject-100)));
    AccCNWPPointInCongruent(subject-100) = nanmean(updatedTable.wordAcc(updatedTable.word ~= updatedTable.color)...
        .*(updatedTable.BlockResponseKeyTime(updatedTable.word ~= updatedTable.color)<InCongWPThresh(subject-100)));
    
    AccCNWPCong(subject-100, :) = (updatedTable.colorAcc(updatedTable.word == updatedTable.color)...
        .*(updatedTable.BlockResponseVoiceTime(updatedTable.word == updatedTable.color)<CongCNThresh(subject-100))) ...
        .* (updatedTable.wordAcc(updatedTable.word == updatedTable.color)...
        .*(updatedTable.BlockResponseKeyTime(updatedTable.word == updatedTable.color)<CongWPThresh(subject-100)));
    
    AccCNWPInCong(subject-100, :) = (updatedTable.colorAcc(updatedTable.word ~= updatedTable.color)...
        .*(updatedTable.BlockResponseVoiceTime(updatedTable.word ~= updatedTable.color)<InCongCNThresh(subject-100))) ...
        .* (updatedTable.wordAcc(updatedTable.word ~= updatedTable.color)...
        .*(updatedTable.BlockResponseKeyTime(updatedTable.word ~= updatedTable.color)<InCongWPThresh(subject-100)));
    
    RawAccCNWPVoiceCongruent(subject-100) = nanmean(updatedTable.wordAcc(updatedTable.word == updatedTable.color));
    RawAccCNWPVoiceInCongruent(subject-100) = nanmean(updatedTable.wordAcc(updatedTable.word ~= updatedTable.color));
    RawAccCNWPPointCongruent(subject-100) = nanmean(updatedTable.wordAcc(updatedTable.word == updatedTable.color));
    RawAccCNWPPointInCongruent(subject-100) = nanmean(updatedTable.wordAcc(updatedTable.word ~= updatedTable.color));
    
    WPColorNamingWordErrors(subject-100) = sum(updatedTable.colorAcc == 0 & updatedTable.voiceresponse == updatedTable.word) ...
        /sum(updatedTable.colorAcc == 0 & updatedTable.color ~= updatedTable.word & updatedTable.voiceresponse > 0);
    for i = 1:length(updatedTable.wordAcc)
        CNWPointColorErrors(i)= wordAcc(i) == 0 & ...
            str2num(SubjectParameters.moves{keyboardMappings(updatedTable.color(i))}) == (keyresponse(i)) ....
            & updatedTable.color(i) ~= updatedTable.word(i);
        
        CNWPointLocationErrors(i)= wordAcc(i) == 0 & ...
            str2num(SubjectParameters.moves{keyboardMappings(updatedTable.location(i))}) == (keyresponse(i)) ....
            & updatedTable.color(i) ~= updatedTable.word(i);
    end
    CNWPColorErrors(subject-100) = sum(CNWPointColorErrors)/ ...
        sum(updatedTable.wordAcc == 0  & updatedTable.color ~= updatedTable.word & keyresponse > 0);
    CNWPLocationErrors(subject-100) = sum(CNWPointLocationErrors)/ ...
        sum(updatedTable.wordAcc == 0  & updatedTable.color ~= updatedTable.word & keyresponse > 0);
    
    matFiles = dir(['SessionData/' num2str(subject) '/MultitaskingCNWP*.mat']);
    numfiles = length(matFiles);
    d = matFiles;
    [dx,dx] = sort([d.datenum],'descend');
    newest = d(dx(1)).name;
    load(['SessionData/' num2str(subject) '/' newest]);
    
    RTCNWPCongruent(subject-100) = nanmean(CSVData.BlockResponseKeyTime(CSVData.word == CSVData.color));
    RTCNWPInCongruent(subject-100) = nanmean(CSVData.BlockResponseKeyTime(CSVData.word ~= CSVData.color));
    
    matFiles = dir(['SessionData/' num2str(subject) '/MultitaskingCNLP*.mat']);
    numfiles = length(matFiles);
    d = matFiles;
    [dx,dx] = sort([d.datenum],'descend');
    newest = d(dx(1)).name;
    load(['SessionData/' num2str(subject) '/' newest]);
    RTCNLPCongruent(subject-100) = nanmean(CSVData.BlockResponseKeyTime(CSVData.word == CSVData.color));
    RTCNLPInCongruent(subject-100) = nanmean(CSVData.BlockResponseKeyTime(CSVData.word ~= CSVData.color));
    
    
end

%%
figure(1)

title('27 subjects')
hold on
% xlim([0 13.5])
subplot(2, 1, 1)
hold on
errorbar(1, nanmean(RTCNCongruent), nanstd(RTCNCongruent)/sqrt(length(RTCNCongruent)))
bar(1, nanmean(RTCNCongruent))
errorbar(2, nanmean(RTLPCongruent), nanstd(RTLPCongruent)/sqrt(length(RTLPCongruent)))
bar(2, nanmean(RTLPCongruent))

errorbar(4, nanmean(RTCNLPVoiceCongruent), nanstd(RTCNLPVoiceCongruent)/sqrt(length(RTCNLPVoiceCongruent)))
bar(4, nanmean(RTCNLPVoiceCongruent))
errorbar(5, nanmean(RTCNLPCongruent), nanstd(RTCNLPCongruent)/sqrt(length(RTCNLPCongruent)))
bar(5, nanmean(RTCNLPCongruent))


errorbar(7, nanmean(RTCNInCongruent), nanstd(RTCNInCongruent)/sqrt(length(RTCNInCongruent)))
bar(7, nanmean(RTCNInCongruent))
errorbar(8, nanmean(RTLPInCongruent), nanstd(RTLPInCongruent)/sqrt(length(RTLPInCongruent)))
bar(8, nanmean(RTLPInCongruent))
errorbar(10, nanmean(RTCNLPVoiceInCongruent), nanstd(RTCNLPVoiceInCongruent)/sqrt(length(RTCNLPVoiceInCongruent)))
bar(10, nanmean(RTCNLPVoiceInCongruent))
errorbar(11, nanmean(RTCNLPInCongruent), nanstd(RTCNLPInCongruent)/sqrt(length(RTCNLPInCongruent)))
bar(11, nanmean(RTCNLPInCongruent))


set(gca, 'fontsize', 14)


ylabel('RT')
xticks([1 2 4 5 7 8 10 11])
xticklabels({'CNCong', 'LPCong', 'CNLPVoice', 'CNLPPoint', 'CNIncong', 'LPIncong', 'CNLPVoice', 'CNLPPoint'})

subplot(2, 1, 2)
hold on
errorbar(1, nanmean(RTCNCongruent), nanstd(RTCNCongruent)/sqrt(length(RTCNCongruent)))
bar(1, nanmean(RTCNCongruent))
errorbar(2, nanmean(RTWPCongruent), nanstd(RTWPCongruent)/sqrt(length(RTWPCongruent)))
bar(2, nanmean(RTWPCongruent))

errorbar(4, nanmean(RTCNWPVoiceCongruent), nanstd(RTCNWPVoiceCongruent)/sqrt(length(RTCNWPVoiceCongruent)))
bar(4, nanmean(RTCNWPVoiceCongruent))
errorbar(5, nanmean(RTCNWPCongruent), nanstd(RTCNWPCongruent)/sqrt(length(RTCNWPCongruent)))
bar(5, nanmean(RTCNWPCongruent))

errorbar(7, nanmean(RTCNInCongruent), nanstd(RTCNInCongruent)/sqrt(length(RTCNInCongruent)))
bar(7, nanmean(RTCNInCongruent))
errorbar(8, nanmean(RTWPInCongruent), nanstd(RTWPInCongruent)/sqrt(length(RTWPInCongruent)))
bar(8, nanmean(RTWPInCongruent))

errorbar(10, nanmean(RTCNWPVoiceInCongruent), nanstd(RTCNWPVoiceInCongruent)/sqrt(length(RTCNWPVoiceInCongruent)))
bar(10, nanmean(RTCNWPVoiceInCongruent))
errorbar(11, nanmean(RTCNWPInCongruent), nanstd(RTCNWPInCongruent)/sqrt(length(RTCNWPInCongruent)))
bar(11, nanmean(RTCNWPInCongruent))


ylabel('RT')
xticks([1 2 4 5 7 8 10 11])
xticklabels({'CNCong', 'WPCong', 'CNWPVoice', 'CNWPPoint', 'CNIncong', 'WPIncong', 'CNWPVoice', 'CNWPPoint'})


set(gca, 'fontsize', 14)
%%
figure(2)

hold on
% xlim([0 13.5])
subplot(2, 1, 1)
hold on
title('multitasking tasks thresholded within task')
errorbar(1, nanmean(AccCNLPVoiceCongruent), nanstd(AccCNLPVoiceCongruent)/sqrt(length(AccCNLPVoiceCongruent)))
bar(1, nanmean(AccCNLPVoiceCongruent))
errorbar(2, nanmean(AccCNLPPointCongruent), nanstd(AccCNLPPointCongruent)/sqrt(length(AccCNLPPointCongruent)))
bar(2, nanmean(AccCNLPPointCongruent))

errorbar(4, nanmean(AccCNLPVoiceInCongruent), nanstd(AccCNLPVoiceInCongruent)/sqrt(length(AccCNLPVoiceInCongruent)))
bar(4, nanmean(AccCNLPVoiceInCongruent))
errorbar(5, nanmean(AccCNLPPointInCongruent), nanstd(AccCNLPPointInCongruent)/sqrt(length(AccCNLPPointInCongruent)))
bar(5, nanmean(AccCNLPPointInCongruent))


set(gca, 'fontsize', 14)


ylabel('accuracy')
xticks([1 2 4 5])
xticklabels({ 'CNLPVoiceCong', 'CNLPPointCong',  'CNLPVoiceInCong', 'CNLPPointInCong'})
ylim([0 1])
grid on

subplot(2, 1, 2)
hold on
errorbar(1, nanmean(AccCNWPVoiceCongruent), nanstd(AccCNWPVoiceCongruent)/sqrt(length(AccCNWPVoiceCongruent)))
bar(1, nanmean(AccCNWPVoiceCongruent))
errorbar(2, nanmean(AccCNWPPointCongruent), nanstd(AccCNWPPointCongruent)/sqrt(length(AccCNWPPointCongruent)))
bar(2, nanmean(AccCNWPPointCongruent))

errorbar(4, nanmean(AccCNWPVoiceInCongruent), nanstd(AccCNWPVoiceInCongruent)/sqrt(length(AccCNWPVoiceInCongruent)))
bar(4, nanmean(AccCNWPVoiceInCongruent))
errorbar(5, nanmean(AccCNWPPointInCongruent), nanstd(AccCNWPPointInCongruent)/sqrt(length(AccCNWPPointInCongruent)))
bar(5, nanmean(AccCNWPPointInCongruent))


set(gca, 'fontsize', 14)
grid on


ylabel('accuracy')
ylim([0 1])
xticks([1 2 4 5])
xticklabels({ 'CNWPVoiceCong', 'CNWPPointCong',  'CNWPVoiceInCong', 'CNWPPointInCong'})


set(gca, 'fontsize', 14)
%%
figure(3)
hold on
% nansubs = [1 3 7 17];
AccCNLPCong(nansubs, :) = NaN;
AccCNLPInCong(nansubs, :) = NaN;
AccCNWPCong(nansubs, :) = NaN;
AccCNWPInCong(nansubs, :) = NaN;

ACCCNLPMinusWPCong = nanmean(AccCNLPCong')-nanmean(AccCNWPCong');
ACCCNLPMinusWPInCong = nanmean(AccCNLPInCong')-nanmean(AccCNWPInCong');


errorbar([1 2 4 5 7 8], nanmean([nanmean(AccCNLPCong'); nanmean(AccCNLPInCong'); ...
    nanmean(AccCNWPCong'); nanmean(AccCNWPInCong'); ...
    (ACCCNLPMinusWPCong); (ACCCNLPMinusWPInCong)]'),...
    ...
    nanstd([nanmean(AccCNLPCong'); nanmean(AccCNLPInCong'); ...
    nanmean(AccCNWPCong'); nanmean(AccCNWPInCong'); ...
    (ACCCNLPMinusWPCong); (ACCCNLPMinusWPInCong)]'), '.k');

bar([1 2 4 5 7 8], nanmean([nanmean(AccCNLPCong'); nanmean(AccCNLPInCong'); ...
    nanmean(AccCNWPCong'); nanmean(AccCNWPInCong'); ...
    ACCCNLPMinusWPCong; ACCCNLPMinusWPInCong; ]'));
ylabel('accuracy')
ylim([0 1])
xticks([1 2 4 5 7 8])
xticklabels({ 'CNLPCong', 'CNLPIncong',  'CNWPCong', 'CNWPInCong', 'LP-WP Cong', 'LP-WP InCong'})


set(gca, 'fontsize', 14)

%%
colors = distinguishable_colors(30);
figure(4)
hold on
subplot(2, 1, 1)
hold on
title('individual and overall differences in multitasking accuracy (each subject thresholded on 2std single task RT)')
violin([nanmean(AccCNLPCong'); nanmean(AccCNLPInCong'); ...
    nanmean(AccCNWPCong'); nanmean(AccCNWPInCong')]', 'bw',.08, 'facecolor',[.9 .9 .9])
for i = 1:30
    % scatter([ones(1, 30)],nanmean(AccCNLPCong'));
    % scatter([2*ones(1, 30)],nanmean(AccCNLPInCong'));
    % scatter([3*ones(1, 30)],nanmean(AccCNWPCong'));
    % scatter([4*ones(1, 30)],nanmean(AccCNWPInCong'));
    if i <=15
        scatter(1+.004*i,nanmean(AccCNLPCong(i, :)), 'MarkerEdgeColor', colors(i, :));
        scatter(2+.004*i,nanmean(AccCNLPInCong(i, :)), 'MarkerEdgeColor', colors(i, :));
        scatter(3+.004*i,nanmean(AccCNWPCong(i, :)), 'MarkerEdgeColor', colors(i, :));
        scatter(4+.004*i,nanmean(AccCNWPInCong(i, :)), 'MarkerEdgeColor', colors(i, :));
    else
        scatter(1-.004*i,nanmean(AccCNLPCong(i, :)), 'MarkerEdgeColor', colors(i, :));
        scatter(2-.004*i,nanmean(AccCNLPInCong(i, :)), 'MarkerEdgeColor', colors(i, :));
        scatter(3-.004*i,nanmean(AccCNWPCong(i, :)), 'MarkerEdgeColor', colors(i, :));
        scatter(4-.004*i,nanmean(AccCNWPInCong(i, :)), 'MarkerEdgeColor', colors(i, :));
    end
end
ylabel('accuracy')
ylim([0 1])
xticks([1 2 3 4 5 6])

xticklabels({ 'CNLPCong', 'CNLPIncong',  'CNWPCong', 'CNWPInCong'})
legend({'Violin', 'Mean', 'Median'})
grid on
set(gca, 'fontsize', 14)

subplot(2, 1, 2)
hold on
title('within subject differences in accuracy between multitasking conditions (each subject thresholded on 2std single task RT)')
violin([(ACCCNLPMinusWPCong); (ACCCNLPMinusWPInCong)]', 'bw',.2, 'facecolor',[.9 .9 .9])
% scatter([1*ones(1, 30)],ACCCNLPMinusWPCong);
% scatter([2*ones(1, 30)],(ACCCNLPMinusWPInCong));
for i = 1:30
    if i <=15
        scatter(1+.004*i,(ACCCNLPMinusWPCong(i)), 'MarkerEdgeColor', colors(i, :));
        scatter(2+.004*i,(ACCCNLPMinusWPInCong(i)), 'MarkerEdgeColor', colors(i, :));
    else
        scatter(1-.004*i,(ACCCNLPMinusWPCong(i)), 'MarkerEdgeColor', colors(i, :));
        scatter(2-.004*i,(ACCCNLPMinusWPInCong(i)), 'MarkerEdgeColor', colors(i, :));
    end
end

ylabel('difference in accuracy')
% ylim([0 1])
xticks([1 2 ])
plot([0 3], [0 0], 'k', 'linewidth', 4)
xticklabels({  'CNLP-CNWP Cong', 'CNLP-CNWP InCong'})
legend({'Violin', 'Mean', 'Median'})
set(gca, 'fontsize', 14)
grid on

%%
figure(5)

% nansubs = [1 3 9 14 16 17 18 23 25 26 28]
%nansubs = [1 3 17];
RawCongCNAcc(nansubs) = NaN;
RawCongLPAcc(nansubs) = NaN;
RawCongWPAcc(nansubs) = NaN;
RawInCongCNAcc(nansubs) = NaN;
RawInCongLPAcc(nansubs) = NaN;
RawInCongWPAcc(nansubs) = NaN;

RawAccCNLPVoiceCongruent(nansubs) = NaN;
RawAccCNLPPointCongruent(nansubs) = NaN;
RawAccCNWPVoiceCongruent(nansubs) = NaN;
RawAccCNWPPointCongruent(nansubs) = NaN;

RawAccCNLPVoiceInCongruent(nansubs) = NaN;
RawAccCNLPPointInCongruent(nansubs) = NaN;
RawAccCNWPVoiceInCongruent(nansubs) = NaN;
RawAccCNWPPointInCongruent(nansubs) = NaN;

% subplot(2, 1, 1)
% title('raw accuracy')
% 
% hold on
% errorbar(1, nanmean(RawCongCNAcc), nanstd(RawCongCNAcc)/sqrt(length(RawCongCNAcc)))
% bar(1, nanmean(RawCongCNAcc))
% errorbar(2, nanmean(RawCongLPAcc), nanstd(RawCongLPAcc)/sqrt(length(RawCongLPAcc)))
% bar(2, nanmean(RawCongLPAcc))
% 
% errorbar(4, nanmean(RawAccCNLPVoiceCongruent), nanstd(RawAccCNLPVoiceCongruent)/sqrt(length(RawAccCNLPVoiceCongruent)))
% bar(4, nanmean(RawAccCNLPVoiceCongruent))
% errorbar(5, nanmean(RawAccCNLPPointCongruent), nanstd(RawAccCNLPPointCongruent)/sqrt(length(RawAccCNLPPointCongruent)))
% bar(5, nanmean(RawAccCNLPPointCongruent))
% 
% 
% 
% errorbar(7, nanmean(RawInCongCNAcc), nanstd(RawInCongCNAcc)/sqrt(length(RawInCongCNAcc)))
% bar(7, nanmean(RawInCongCNAcc))
% errorbar(8, nanmean(RawInCongLPAcc), nanstd(RawInCongLPAcc)/sqrt(length(RawInCongLPAcc)))
% bar(8, nanmean(RawInCongLPAcc))
% 
% errorbar(10, nanmean(RawAccCNLPVoiceInCongruent), nanstd(RawAccCNLPVoiceInCongruent)/sqrt(length(RawAccCNLPVoiceInCongruent)))
% bar(10, nanmean(RawAccCNLPVoiceInCongruent))
% errorbar(11, nanmean(RawAccCNLPPointInCongruent), nanstd(RawAccCNLPPointInCongruent)/sqrt(length(RawAccCNLPPointInCongruent)))
% bar(11, nanmean(RawAccCNLPPointInCongruent))
% grid on
% set(gca, 'fontsize', 14)
% 
% 
% ylabel('Accuracy')
% xticks([1 2 4 5 7 8 10 11])
% xticklabels({'CNCong', 'LPCong', 'CNLPVoice', 'CNLPPoint', 'CNIncong', 'LPIncong', 'CNLPVoice', 'CNLPPoint'})
% 
% subplot(2, 1, 2)
% hold on
% errorbar(1, nanmean(RawCongCNAcc), nanstd(RawCongCNAcc)/sqrt(length(RawCongCNAcc)))
% bar(1, nanmean(RawCongCNAcc))
% errorbar(2, nanmean(RawCongWPAcc), nanstd(RawCongWPAcc)/sqrt(length(RawCongWPAcc)))
% bar(2, nanmean(RawCongWPAcc))
% 
% errorbar(4, nanmean(RawAccCNWPVoiceCongruent), nanstd(RawAccCNWPVoiceCongruent)/sqrt(length(RawAccCNWPVoiceCongruent)))
% bar(4, nanmean(RawAccCNWPVoiceCongruent))
% errorbar(5, nanmean(RawAccCNWPPointCongruent), nanstd(RawAccCNWPPointCongruent)/sqrt(length(RawAccCNWPPointCongruent)))
% bar(5, nanmean(RawAccCNWPPointCongruent))
% 
% 
% errorbar(7, nanmean(RawInCongCNAcc), nanstd(RawInCongCNAcc)/sqrt(length(RawInCongCNAcc)))
% bar(7, nanmean(RawInCongCNAcc))
% errorbar(8, nanmean(RawInCongWPAcc), nanstd(RawInCongWPAcc)/sqrt(length(RawInCongWPAcc)))
% bar(8, nanmean(RawInCongWPAcc))
% 
% errorbar(10, nanmean(RawAccCNWPVoiceInCongruent), nanstd(RawAccCNWPVoiceInCongruent)/sqrt(length(RawAccCNWPVoiceInCongruent)))
% bar(10, nanmean(RawAccCNWPVoiceInCongruent))
% errorbar(11, nanmean(RawAccCNWPPointInCongruent), nanstd(RawAccCNWPPointInCongruent)/sqrt(length(RawAccCNWPPointInCongruent)))
% bar(11, nanmean(RawAccCNWPPointInCongruent))
% grid on
% set(gca, 'fontsize', 14)
% 
% 
% ylabel('Accuracy')
% xticks([1 2 4 5 7 8 10 11])
% xticklabels({'CNCong', 'WPCong', 'CNWPVoice', 'CNWPPoint', 'CNIncong', 'WPIncong', 'CNWPVoice', 'CNWPPoint'})
% 
% 
% set(gca, 'fontsize', 14)

% %%
% figure(6)
% hold on
% % nansubs = [1     3     9    14    16    17    18    23    25    26    28];
% ColorNamingWordErrors(nansubs) = NaN;
% LocationPointingWordErrors(nansubs) = NaN;
% WordPointingColorErrors(nansubs) = NaN;
% 
% LPColorNamingWordErrors(nansubs) = NaN;
% WPColorNamingWordErrors(nansubs) = NaN;
% 
% CNLPWordErrors(nansubs) = NaN;
% CNWPColorErrors(nansubs) = NaN;
% 
% colors = summer(30);
% errorbar(1, nanmean(ColorNamingWordErrors), nanstd(ColorNamingWordErrors)/sqrt(length(ColorNamingWordErrors)))
% bar(1, nanmean(ColorNamingWordErrors))
% errorbar(2, nanmean(LocationPointingWordErrors), nanstd(LocationPointingWordErrors)/sqrt(length(LocationPointingWordErrors)))
% bar(2, nanmean(LocationPointingWordErrors))
% errorbar(3, nanmean(WordPointingColorErrors), nanstd(WordPointingColorErrors)/sqrt(length(WordPointingColorErrors)))
% bar(3, nanmean(WordPointingColorErrors))
% 
% errorbar(5, nanmean(LPColorNamingWordErrors), nanstd(LPColorNamingWordErrors)/sqrt(length(LPColorNamingWordErrors)))
% bar(5, nanmean(LPColorNamingWordErrors))
% errorbar(6, nanmean(CNLPWordErrors), nanstd(CNLPWordErrors)/sqrt(length(CNLPWordErrors)))
% bar(6, nanmean(CNLPWordErrors))
% errorbar(7, nanmean(CNLPColorErrors), nanstd(CNLPColorErrors)/sqrt(length(CNLPColorErrors)))
% bar(7, nanmean(CNLPColorErrors))
% 
% errorbar(9, nanmean(WPColorNamingWordErrors), nanstd(WPColorNamingWordErrors)/sqrt(length(WPColorNamingWordErrors)))
% bar(9, nanmean(WPColorNamingWordErrors))
% errorbar(10, nanmean(CNWPColorErrors), nanstd(CNWPColorErrors)/sqrt(length(CNWPColorErrors)))
% bar(10, nanmean(CNWPColorErrors))
% errorbar(11, nanmean(CNWPLocationErrors), nanstd(CNWPLocationErrors)/sqrt(length(CNWPLocationErrors)))
% bar(11, nanmean(CNWPLocationErrors))
% 
% plot([0 12], [.33 .33], '-k', 'linewidth', 2)
% 
% set(gca, 'fontsize', 14)
% 
% ax = gca;
% set(gca,'XTick',[1 2 3 5 6 7 9 10 11])
% set(gca,'XTickLabel',{'CN:WordSubs', 'LP:WordSubs', 'WP:ColorSubs', ...
%     'CNLP:CN WordSubs', 'CNLP:LP WordSubs', 'CNLP:LP ColorSubs',...
%     'CNWP:CN WordSubs', 'CNWP:WP ColorSubs', 'CNWP:WP LocationSubs'})
% 
% ylabel('% errors from category of interest')
% 
% set(gca, 'fontsize', 14)
% xtickangle(45)
% grid on
% %%
% figure(7)
% hold on
% title('Subjects can only multitask with non-overlapping stimuli')
% 
% colors = summer(30);
% 
% 
% bar([1 2 4 5], [nanmean(nanmean(AccCNLPCong(goodsubjects, :)')), nanmean(nanmean(AccCNWPCong(goodsubjects, :)')), ...
%     nanmean(nanmean(AccCNLPInCong(goodsubjects, :)')), ...
%     nanmean(nanmean(AccCNWPInCong(goodsubjects, :)'))], 'FaceColor', [.9 .9 .9], 'EdgeColor', [.95 .95 .95])
% 
% % nanmean(RawCongLPAcc) nanmean(RawInCongCNAcc)
% for i = 1:30
%     if sum(goodsubjects == i) == 1
%         plot([1 2], [nanmean(AccCNLPCong(i, :)) nanmean(AccCNWPCong(i, :))], 'Color', [.6 .6 .6])
%         plot([4 5], [nanmean(AccCNLPInCong(i, :)), nanmean(AccCNWPInCong(i, :))], 'Color', [.6 .6 .6])
%         scatter(1,nanmean(AccCNLPCong(i, :)), 'MarkerEdgeColor', colors(i, :), 'MarkerFaceColor', colors(i, :));
%         scatter(4,nanmean(AccCNLPInCong(i, :)), 'MarkerEdgeColor', colors(i, :), 'MarkerFaceColor', colors(i, :));
%         scatter(2,nanmean(AccCNWPCong(i, :)), 'MarkerEdgeColor', colors(i, :), 'MarkerFaceColor', colors(i, :));
%         scatter(5,nanmean(AccCNWPInCong(i, :)),'MarkerEdgeColor', colors(i, :), 'MarkerFaceColor', colors(i, :));
%     end
% end
% ylabel('accuracy (RT must be within 2 stds of single task)')
% ylim([0 1])
% xticks([1 2 4 5])
% xlabel('congruent trials                         incongruent trials')
% 
% xticklabels({ 'CNLP',   'CNWP', 'CNLP', 'CNWP'})
% % legend({'Violin', 'Mean', 'Median'})
% grid on
% set(gca, 'fontsize', 14)
% 
% 
% % bar([1 2 3], [nanmean(RawInCongCNAcc(goodsubjects)), nanmean(RawInCongLPAcc(goodsubjects)), nanmean(RawInCongWPAcc(goodsubjects))], 'FaceColor', [.9 .9 .9], 'EdgeColor', [.95 .95 .95])
% % 
% % for i = 1:30
% %     if RawCongWPAcc(i) > 0
% %         scatter(1,RawCongCNAcc(i)-i*.001, 'MarkerEdgeColor', colors(i, :), 'MarkerFaceColor', colors(i, :));
% %         scatter(2,RawCongLPAcc(i)-i*.001, 'MarkerEdgeColor', colors(i, :), 'MarkerFaceColor', colors(i, :));
% %         scatter(3,RawCongWPAcc(i)-i*.001, 'MarkerEdgeColor', colors(i, :), 'MarkerFaceColor', colors(i, :));
% %          plot([ 1 2 3], [RawCongCNAcc(i)-i*.001, RawCongLPAcc(i)-i*.001, RawCongWPAcc(i)-i*.001], 'Color', [.6 .6 .6])
% %     end
% % end
% %%
% %%
% figure(8)
% hold on
% % nansubs = [1 3 17];
% ColorNamingWordErrors(nansubs) = NaN;
% LocationPointingWordErrors(nansubs) = NaN;
% WordPointingColorErrors(nansubs) = NaN;
% 
% LPColorNamingWordErrors(nansubs) = NaN;
% WPColorNamingWordErrors(nansubs) = NaN;
% 
% CNLPWordErrors(nansubs) = NaN;
% CNWPColorErrors(nansubs) = NaN;
% 
% colors = summer(30);
% 
% bar(1, nanmean(LPColorNamingWordErrors), 'FaceColor', [.9 .9 .9], 'EdgeColor', [.95 .95 .95])
% bar(2, nanmean(CNLPWordErrors), 'FaceColor', [.9 .9 .9], 'EdgeColor', [.95 .95 .95])
% bar(3, nanmean(CNLPColorErrors), 'FaceColor', [.9 .9 .9], 'EdgeColor', [.95 .95 .95])
% 
% bar(5, nanmean(CNWPColorErrors), 'FaceColor', [.9 .9 .9], 'EdgeColor', [.95 .95 .95])
% bar(6, nanmean(CNWPLocationErrors), 'FaceColor', [.9 .9 .9], 'EdgeColor', [.95 .95 .95])
% bar(7, nanmean(WPColorNamingWordErrors), 'FaceColor', [.9 .9 .9], 'EdgeColor', [.95 .95 .95])
% 
% plot([0 8], [.33 .33], '-k', 'linewidth', 2)
% 
% for i = 1:30
%     scatter(1,nanmean(LPColorNamingWordErrors(i)) + i *.001, 'MarkerEdgeColor', colors(i, :), 'MarkerFaceColor', colors(i, :));
%     scatter(2, nanmean(CNLPWordErrors(i)) + i *.001, 'MarkerEdgeColor', colors(i, :), 'MarkerFaceColor', colors(i, :));
%     scatter(3,nanmean(CNLPColorErrors(i)) + i *.001, 'MarkerEdgeColor', colors(i, :), 'MarkerFaceColor', colors(i, :));
%     
%     scatter(5, nanmean(CNLPWordErrors(i)) + i *.001, 'MarkerEdgeColor', colors(i, :), 'MarkerFaceColor', colors(i, :));
%     
%     scatter(6,nanmean(CNLPColorErrors(i)) + i *.001, 'MarkerEdgeColor', colors(i, :), 'MarkerFaceColor', colors(i, :));
%     scatter(7,nanmean(LPColorNamingWordErrors(i))+ i *.001, 'MarkerEdgeColor', colors(i, :), 'MarkerFaceColor', colors(i, :));
%     plot([3 5], [(nanmean(CNLPColorErrors(i))+ i *.001), ( i *.001+ nanmean(CNLPWordErrors(i)))], 'Color', [.6 .6 .6])
% end
% 
% 
% 
% set(gca, 'fontsize', 14)
% 
% ax = gca;
% set(gca,'XTick',[1 2 3 5 6 7])
% set(gca,'XTickLabel',{ ...
%     'CNLP:CN WordSubs',  'CNLP:LP WordSubs', 'CNLP:LP ColorSubs',...
%     'CNWP:WP ColorSubs', 'CNWP:WP LocationSubs', 'CNWP:CN WordSubs',})
% 
% ylabel('% errors from category of interest')
% 
% set(gca, 'fontsize', 14)
% xtickangle(45)
% grid on
% %%
% 
% %%
% figure(9)
% hold on
% title('Accuracy on single tasks')
% 
% colors = summer(30);
% 
% % goodsubjects = goodsubjects -100;
% bar([1 2 3], [nanmean(RawInCongCNAcc(goodsubjects)), nanmean(RawInCongLPAcc(goodsubjects)), nanmean(RawInCongWPAcc(goodsubjects))], 'FaceColor', [.9 .9 .9], 'EdgeColor', [.95 .95 .95])
% 
% for i = 1:30
%     if RawCongWPAcc(i) > 0
%         scatter(1,RawCongCNAcc(i)-i*.001, 'MarkerEdgeColor', colors(i, :), 'MarkerFaceColor', colors(i, :));
%         scatter(2,RawCongLPAcc(i)-i*.001, 'MarkerEdgeColor', colors(i, :), 'MarkerFaceColor', colors(i, :));
%         scatter(3,RawCongWPAcc(i)-i*.001, 'MarkerEdgeColor', colors(i, :), 'MarkerFaceColor', colors(i, :));
%          plot([ 1 2 3], [RawCongCNAcc(i)-i*.001, RawCongLPAcc(i)-i*.001, RawCongWPAcc(i)-i*.001], 'Color', [.6 .6 .6])
%     end
% end
% ylabel('accuracy')
% ylim([0 1])
% xticks([1 2 3])
% 
% xticklabels({ 'ColorNaming',   'LocationPointing', 'WordPointing'})
% % legend({'Violin', 'Mean', 'Median'})
% grid on
% ylim([0 1.02])
% set(gca, 'fontsize', 14)
%%
%%
figure(10)

% nansubs = [1 3 9 14 16 17 18 23 25 26 28]
%nansubs = [1 3 17];
RawCongCNAcc(nansubs) = NaN;
RawCongLPAcc(nansubs) = NaN;
RawCongWPAcc(nansubs) = NaN;
RawInCongCNAcc(nansubs) = NaN;
RawInCongLPAcc(nansubs) = NaN;
RawInCongWPAcc(nansubs) = NaN;

RawAccCNLPVoiceCongruent(nansubs) = NaN;
RawAccCNLPPointCongruent(nansubs) = NaN;
RawAccCNWPVoiceCongruent(nansubs) = NaN;
RawAccCNWPPointCongruent(nansubs) = NaN;

RawAccCNLPVoiceInCongruent(nansubs) = NaN;
RawAccCNLPPointInCongruent(nansubs) = NaN;
RawAccCNWPVoiceInCongruent(nansubs) = NaN;
RawAccCNWPPointInCongruent(nansubs) = NaN;

subplot(2, 1, 1)
title('raw accuracy')

hold on
bar(1, nanmean(RawCongCNAcc),'FaceColor', [.9 .9 .9], 'EdgeColor', [.95 .95 .95])
bar(2, nanmean(RawCongLPAcc), 'FaceColor', [.9 .9 .9], 'EdgeColor', [.95 .95 .95])
bar(4, nanmean(RawAccCNLPVoiceCongruent), 'FaceColor', [.9 .9 .9], 'EdgeColor', [.95 .95 .95])
bar(5, nanmean(RawAccCNLPPointCongruent), 'FaceColor', [.9 .9 .9], 'EdgeColor', [.95 .95 .95])

bar(7, nanmean(RawInCongCNAcc), 'FaceColor', [.9 .9 .9], 'EdgeColor', [.95 .95 .95])
bar(8, nanmean(RawInCongLPAcc), 'FaceColor', [.9 .9 .9], 'EdgeColor', [.95 .95 .95])

bar(10, nanmean(RawAccCNLPVoiceInCongruent), 'FaceColor', [.9 .9 .9], 'EdgeColor', [.95 .95 .95])
bar(11, nanmean(RawAccCNLPPointInCongruent), 'FaceColor', [.9 .9 .9], 'EdgeColor', [.95 .95 .95])

for i = 1:30
    if RawCongWPAcc(i) > 0
        scatter(1,RawCongCNAcc(i)-i*.001, 'MarkerEdgeColor', colors(i, :), 'MarkerFaceColor', colors(i, :));
        scatter(2,RawCongLPAcc(i)-i*.001, 'MarkerEdgeColor', colors(i, :), 'MarkerFaceColor', colors(i, :));
        scatter(4,RawAccCNLPVoiceCongruent(i)-i*.001, 'MarkerEdgeColor', colors(i, :), 'MarkerFaceColor', colors(i, :));
        scatter(5,RawAccCNLPPointCongruent(i)-i*.001, 'MarkerEdgeColor', colors(i, :), 'MarkerFaceColor', colors(i, :));
        scatter(7,RawInCongCNAcc(i)-i*.001, 'MarkerEdgeColor', colors(i, :), 'MarkerFaceColor', colors(i, :));
        scatter(8,RawInCongLPAcc(i)-i*.001, 'MarkerEdgeColor', colors(i, :), 'MarkerFaceColor', colors(i, :));
        scatter(10,RawAccCNLPVoiceInCongruent(i)-i*.001, 'MarkerEdgeColor', colors(i, :), 'MarkerFaceColor', colors(i, :));
                scatter(11,RawAccCNLPPointInCongruent(i)-i*.001, 'MarkerEdgeColor', colors(i, :), 'MarkerFaceColor', colors(i, :));
    end
end

grid on
set(gca, 'fontsize', 14)


ylabel('Accuracy')
xticks([1 2 4 5 7 8 10 11])
xticklabels({'CNCong', 'LPCong', 'CNLPVoice', 'CNLPPoint', 'CNIncong', 'LPIncong', 'CNLPVoice', 'CNLPPoint'})
xtickangle(30)
subplot(2, 1, 2)
hold on
bar(1, nanmean(RawCongCNAcc), 'FaceColor', [.9 .9 .9], 'EdgeColor', [.95 .95 .95])
bar(2, nanmean(RawCongWPAcc), 'FaceColor', [.9 .9 .9], 'EdgeColor', [.95 .95 .95])

bar(4, nanmean(RawAccCNWPVoiceCongruent), 'FaceColor', [.9 .9 .9], 'EdgeColor', [.95 .95 .95])
bar(5, nanmean(RawAccCNWPPointCongruent) ,'FaceColor', [.9 .9 .9], 'EdgeColor', [.95 .95 .95])

bar(7, nanmean(RawInCongCNAcc),'FaceColor', [.9 .9 .9], 'EdgeColor', [.95 .95 .95])
bar(8, nanmean(RawInCongWPAcc), 'FaceColor', [.9 .9 .9], 'EdgeColor', [.95 .95 .95])

bar(10, nanmean(RawAccCNWPVoiceInCongruent), 'FaceColor', [.9 .9 .9], 'EdgeColor', [.95 .95 .95])
bar(11, nanmean(RawAccCNWPPointInCongruent), 'FaceColor', [.9 .9 .9], 'EdgeColor', [.95 .95 .95])
grid on
set(gca, 'fontsize', 14)

for i = 1:30
    if RawCongWPAcc(i) > 0
        scatter(1,RawCongCNAcc(i)-i*.001, 'MarkerEdgeColor', colors(i, :), 'MarkerFaceColor', colors(i, :));
        scatter(2,RawCongWPAcc(i)-i*.001, 'MarkerEdgeColor', colors(i, :), 'MarkerFaceColor', colors(i, :));
        scatter(4,RawAccCNWPVoiceCongruent(i)-i*.001, 'MarkerEdgeColor', colors(i, :), 'MarkerFaceColor', colors(i, :));
        scatter(5,RawAccCNWPPointCongruent(i)-i*.001, 'MarkerEdgeColor', colors(i, :), 'MarkerFaceColor', colors(i, :));
        scatter(7,RawInCongCNAcc(i)-i*.001, 'MarkerEdgeColor', colors(i, :), 'MarkerFaceColor', colors(i, :));
        scatter(8,RawInCongWPAcc(i)-i*.001, 'MarkerEdgeColor', colors(i, :), 'MarkerFaceColor', colors(i, :));
        scatter(10,RawAccCNWPVoiceInCongruent(i)-i*.001, 'MarkerEdgeColor', colors(i, :), 'MarkerFaceColor', colors(i, :));
                scatter(11,RawAccCNWPPointInCongruent(i)-i*.001, 'MarkerEdgeColor', colors(i, :), 'MarkerFaceColor', colors(i, :));
    end
end


ylabel('Accuracy')
xticks([1 2 4 5 7 8 10 11])
xticklabels({'CNCong', 'WPCong', 'CNWPVoice', 'CNWPPoint', 'CNIncong', 'WPIncong', 'CNWPVoice', 'CNWPPoint'})
xtickangle(30)

set(gca, 'fontsize', 14)
%%
figure(11)

% nansubs = [1 3 9 14 16 17 18 23 25 26 28]
%nansubs = [1 3 17];
%single tasks
RTCNCongruent(nansubs) = NaN;
RTLPCongruent(nansubs) = NaN;
RTWPCongruent(nansubs) = NaN;

RTCNInCongruent(nansubs) = NaN;
RTLPInCongruent(nansubs) = NaN;
RTWPInCongruent(nansubs) = NaN;

%multitasks voice
RTCNLPVoiceCongruent(nansubs) = NaN;
RTCNWPVoiceCongruent(nansubs) = NaN;
RTCNLPVoiceInCongruent(nansubs) = NaN;
RTCNWPVoiceInCongruent(nansubs) = NaN;

%multitaskspointing
RTCNWPCongruent(nansubs) = NaN;
RTCNWPInCongruent(nansubs) = NaN;
RTCNLPCongruent(nansubs) = NaN;
RTCNLPInCongruent(nansubs) = NaN;

subplot(2, 1, 1)
title('RT')

hold on
bar(1, nanmean(RTCNCongruent),'FaceColor', [.9 .9 .9], 'EdgeColor', [.95 .95 .95])
bar(2, nanmean(RTLPCongruent), 'FaceColor', [.9 .9 .9], 'EdgeColor', [.95 .95 .95])
bar(4, nanmean(RTCNLPVoiceCongruent), 'FaceColor', [.9 .9 .9], 'EdgeColor', [.95 .95 .95])
bar(5, nanmean(RTCNLPCongruent), 'FaceColor', [.9 .9 .9], 'EdgeColor', [.95 .95 .95])

bar(7, nanmean(RTCNInCongruent), 'FaceColor', [.9 .9 .9], 'EdgeColor', [.95 .95 .95])
bar(8, nanmean(RTCNLPInCongruent), 'FaceColor', [.9 .9 .9], 'EdgeColor', [.95 .95 .95])
bar(10, nanmean(RTCNLPVoiceInCongruent), 'FaceColor', [.9 .9 .9], 'EdgeColor', [.95 .95 .95])
bar(11, nanmean(RTCNInCongruent), 'FaceColor', [.9 .9 .9], 'EdgeColor', [.95 .95 .95])

for i = goodsubjects
        scatter(1,RTCNCongruent(i)-i*.001, 'MarkerEdgeColor', colors(i, :), 'MarkerFaceColor', colors(i, :));
        scatter(2,RTLPCongruent(i)-i*.001, 'MarkerEdgeColor', colors(i, :), 'MarkerFaceColor', colors(i, :));
        scatter(4,RTCNLPVoiceCongruent(i)-i*.001, 'MarkerEdgeColor', colors(i, :), 'MarkerFaceColor', colors(i, :));
        scatter(5,RTCNLPCongruent(i)-i*.001, 'MarkerEdgeColor', colors(i, :), 'MarkerFaceColor', colors(i, :));
        scatter(7,RTCNInCongruent(i)-i*.001, 'MarkerEdgeColor', colors(i, :), 'MarkerFaceColor', colors(i, :));
        scatter(8,RTCNLPInCongruent(i)-i*.001, 'MarkerEdgeColor', colors(i, :), 'MarkerFaceColor', colors(i, :));
        scatter(10,RTCNLPVoiceInCongruent(i)-i*.001, 'MarkerEdgeColor', colors(i, :), 'MarkerFaceColor', colors(i, :));
                scatter(11,RTCNLPInCongruent(i)-i*.001, 'MarkerEdgeColor', colors(i, :), 'MarkerFaceColor', colors(i, :));
end

grid on
set(gca, 'fontsize', 14)


ylabel('RT')
xticks([1 2 4 5 7 8 10 11])
xticklabels({'CNCong', 'LPCong', 'CNLPVoice', 'CNLPPoint', 'CNIncong', 'LPIncong', 'CNLPVoice', 'CNLPPoint'})
xtickangle(30)
ylim([.2 1.2])
subplot(2, 1, 2)
hold on

bar(1, nanmean(RTCNCongruent),'FaceColor', [.9 .9 .9], 'EdgeColor', [.95 .95 .95])
bar(2, nanmean(RTWPCongruent), 'FaceColor', [.9 .9 .9], 'EdgeColor', [.95 .95 .95])
bar(4, nanmean(RTCNWPVoiceCongruent), 'FaceColor', [.9 .9 .9], 'EdgeColor', [.95 .95 .95])
bar(5, nanmean(RTCNWPCongruent), 'FaceColor', [.9 .9 .9], 'EdgeColor', [.95 .95 .95])

bar(7, nanmean(RTCNInCongruent), 'FaceColor', [.9 .9 .9], 'EdgeColor', [.95 .95 .95])
bar(8, nanmean(RTCNWPInCongruent), 'FaceColor', [.9 .9 .9], 'EdgeColor', [.95 .95 .95])
bar(10, nanmean(RTCNWPVoiceInCongruent), 'FaceColor', [.9 .9 .9], 'EdgeColor', [.95 .95 .95])
bar(11, nanmean(RTCNInCongruent), 'FaceColor', [.9 .9 .9], 'EdgeColor', [.95 .95 .95])

for i = goodsubjects
        scatter(1,RTCNCongruent(i)-i*.001, 'MarkerEdgeColor', colors(i, :), 'MarkerFaceColor', colors(i, :));
        scatter(2,RTWPCongruent(i)-i*.001, 'MarkerEdgeColor', colors(i, :), 'MarkerFaceColor', colors(i, :));
        scatter(4,RTCNWPVoiceCongruent(i)-i*.001, 'MarkerEdgeColor', colors(i, :), 'MarkerFaceColor', colors(i, :));
        scatter(5,RTCNWPCongruent(i)-i*.001, 'MarkerEdgeColor', colors(i, :), 'MarkerFaceColor', colors(i, :));
        scatter(7,RTCNInCongruent(i)-i*.001, 'MarkerEdgeColor', colors(i, :), 'MarkerFaceColor', colors(i, :));
        scatter(8,RTCNWPInCongruent(i)-i*.001, 'MarkerEdgeColor', colors(i, :), 'MarkerFaceColor', colors(i, :));
        scatter(10,RTCNWPVoiceInCongruent(i)-i*.001, 'MarkerEdgeColor', colors(i, :), 'MarkerFaceColor', colors(i, :));
                scatter(11,RTCNWPInCongruent(i)-i*.001, 'MarkerEdgeColor', colors(i, :), 'MarkerFaceColor', colors(i, :));
end

grid on

ylabel('RT')
xticks([1 2 4 5 7 8 10 11])
xticklabels({'CNCong', 'WPCong', 'CNWPVoice', 'CNWPPoint', 'CNIncong', 'WPIncong', 'CNWPVoice', 'CNWPPoint'})
xtickangle(30)
ylim([.2 1.2])
set(gca, 'fontsize', 14)
%%

%%
figure(20)
hold on
title('Subjects can only multitask with non-overlapping stimuli')

colors = summer(30);


bar([1 2 3 7 6 5 9 10 11 15 14 13], [nanmean(RawCongCNAcc(goodsubjects)), ...
     nanmean(RawCongLPAcc(goodsubjects)), ...
    nanmean(nanmean(AccCNLPCong(goodsubjects, :)')), ...
    ...
    nanmean(RawCongCNAcc(goodsubjects)), ...
        nanmean(RawCongWPAcc(goodsubjects)), ...
    nanmean(nanmean(AccCNWPCong(goodsubjects, :)')), ...
    ...
    nanmean(RawInCongCNAcc(goodsubjects)), ... %9
     nanmean(RawInCongLPAcc(goodsubjects)), ...
    nanmean(nanmean(AccCNLPInCong(goodsubjects, :)')), ...
    ....
        nanmean(RawInCongCNAcc(goodsubjects)), ... 
    nanmean(RawInCongWPAcc(goodsubjects)), ...
    nanmean(nanmean(AccCNWPInCong(goodsubjects, :)'))], 'FaceColor', [.9 .9 .9], 'EdgeColor', [.95 .95 .95])
xticks([1 2 3 5 6 7 9 10 11 13 14 15]);
xticklabels({ 'CN', 'LP', 'CNLP',   'CNWP','WP','CN',  'CN', 'LP', 'CNLP', 'CNWP' ,'WP','CN' })

% nanmean(RawCongLPAcc) nanmean(RawInCongCNAcc)
 %plot([1 2 3 5 6 7 9 10 11 13 14 15],
for i = 1:30
    if sum(goodsubjects == i) == 1
        plot([1 2 3 5], ...
    [nanmean(RawCongCNAcc(i))-i*.0005, ...
     nanmean(RawCongLPAcc(i))-i*.0005, ...
    nanmean(AccCNLPCong(i, :))+i*.0005, ...
    nanmean(AccCNWPCong(i, :))+i*.0005],'Color', [.6 .6 .6])

        plot([7 6 5], ...
    [nanmean(RawCongCNAcc(i))-i*.0005, ...
     nanmean(RawCongWPAcc(i))-i*.0005, ...
    nanmean(AccCNWPCong(i, :))+i*.0005],'Color', [.6 .6 .6])


        plot([9 10 11 13], ...
    [nanmean(RawInCongCNAcc(i))-i*.0005, ...
     nanmean(RawInCongLPAcc(i))-i*.0005, ...
    nanmean(AccCNLPInCong(i, :))+i*.0005, ...
    nanmean(AccCNWPInCong(i, :))+i*.0005],'Color', [.6 .6 .6])

        plot([15 14 13], ...
    [nanmean(RawInCongCNAcc(i))-i*.0005, ...
     nanmean(RawInCongWPAcc(i))-i*.0005, ...
    nanmean(AccCNWPInCong(i, :))+i*.0005],'Color', [.6 .6 .6])

        scatter(1, nanmean(RawCongCNAcc(i))-i*.0005, 'MarkerEdgeColor', colors(i, :), 'MarkerFaceColor', colors(i, :));
        scatter(2,nanmean(RawCongLPAcc(i))-i*.0005, 'MarkerEdgeColor', colors(i, :), 'MarkerFaceColor', colors(i, :));
        scatter(3,nanmean(AccCNLPCong(i, :))+i*.0005, 'MarkerEdgeColor', colors(i, :), 'MarkerFaceColor', colors(i, :));
        scatter(7, nanmean(RawCongCNAcc(i))-i*.0005, 'MarkerEdgeColor', colors(i, :), 'MarkerFaceColor', colors(i, :));
        scatter(6,nanmean(RawCongWPAcc(i))-i*.0005, 'MarkerEdgeColor', colors(i, :), 'MarkerFaceColor', colors(i, :));
        scatter(5,nanmean(AccCNWPCong(i, :))+i*.0005,'MarkerEdgeColor', colors(i, :), 'MarkerFaceColor', colors(i, :));
        scatter(9, nanmean(RawInCongCNAcc(i))-i*.0005, 'MarkerEdgeColor', colors(i, :), 'MarkerFaceColor', colors(i, :));
        scatter(10, nanmean(RawInCongLPAcc(i))-i*.0005, 'MarkerEdgeColor', colors(i, :), 'MarkerFaceColor', colors(i, :));
        scatter(11, nanmean(AccCNLPInCong(i, :))+i*.0005, 'MarkerEdgeColor', colors(i, :), 'MarkerFaceColor', colors(i, :));
       scatter(15, nanmean(RawInCongCNAcc(i))-i*.0005, 'MarkerEdgeColor', colors(i, :), 'MarkerFaceColor', colors(i, :));
        scatter(14, nanmean(RawInCongWPAcc(i))-i*.0005, 'MarkerEdgeColor', colors(i, :), 'MarkerFaceColor', colors(i, :));
        scatter(13,nanmean(AccCNWPInCong(i, :))+i*.0005, 'MarkerEdgeColor', colors(i, :), 'MarkerFaceColor', colors(i, :));
    end
end
ylabel('accuracy (RT must be within 2 stds of single task)')
ylim([0 1.01])

xlabel('congruent trials                              incongruent trials')


% legend({'Violin', 'Mean', 'Median'})
grid on
set(gca, 'fontsize', 24)

%% FINAL ACCURACY PLOT (SEBASTIAN)

plotAccuracy = 1;

nSubj = length(goodsubjects);

fig1 = figure(20);
set(fig1, 'Position', [100 100 900 250]);

% colors
plotSettings;
barColor = colors(cWeak, :);
scatterColor = colors(cContrast1, :);

% labels
xLabels = {'CN', 'LP', 'WP', 'CN+LP', 'CN+WP'};

% plot
for plotID = 1:2
    
    subplot(1,2,plotID);

    x = 1:5;

    scatterdata_x = repmat(x, nSubj, 1);
    scatterdata_y = nan(size(scatterdata_x));
    
    if(plotID == 1)
        bardata_mean = [nanmean(RawCongCNAcc(goodsubjects)), ...
                         nanmean(RawCongLPAcc(goodsubjects)), ...
                         nanmean(RawCongWPAcc(goodsubjects)), ...
                         nanmean(nanmean(AccCNLPCong(goodsubjects, :)')), ...
                         nanmean(nanmean(AccCNWPCong(goodsubjects, :)'))] * 100; 

        bardata_sem =   [nanstd(RawCongCNAcc(goodsubjects)), ...
                         nanstd(RawCongLPAcc(goodsubjects)), ...
                         nanstd(RawCongWPAcc(goodsubjects)), ...
                         nanstd(nanmean(AccCNLPCong(goodsubjects, :)')), ...
                         nanstd(nanmean(AccCNWPCong(goodsubjects, :)'))] * 100 / sqrt(nSubj);  
                     
        scatterdata_y(:,1) = RawCongCNAcc(goodsubjects) * 100;
        scatterdata_y(:,2) = RawCongLPAcc(goodsubjects) * 100;
        scatterdata_y(:,3) = RawCongWPAcc(goodsubjects) * 100;
        scatterdata_y(:,4) = nanmean(AccCNLPCong(goodsubjects, :)') * 100;
        scatterdata_y(:,5) = nanmean(AccCNWPCong(goodsubjects, :)') * 100;
    
    elseif(plotID == 2)            
        bardata_mean = [nanmean(RawInCongCNAcc(goodsubjects)), ...
                     nanmean(RawInCongLPAcc(goodsubjects)), ...
                     nanmean(RawInCongWPAcc(goodsubjects)), ...
                     nanmean(nanmean(AccCNLPInCong(goodsubjects, :)')), ...
                     nanmean(nanmean(AccCNWPInCong(goodsubjects, :)'))] * 100; 

        bardata_sem =   [nanstd(RawInCongCNAcc(goodsubjects)), ...
                         nanstd(RawInCongLPAcc(goodsubjects)), ...
                         nanstd(RawInCongWPAcc(goodsubjects)), ...
                         nanstd(nanmean(AccCNLPInCong(goodsubjects, :)')), ...
                         nanstd(nanmean(AccCNWPInCong(goodsubjects, :)'))] * 100 / sqrt(nSubj);  
                     
        scatterdata_y(:,1) = RawInCongCNAcc(goodsubjects) * 100;
        scatterdata_y(:,2) = RawInCongLPAcc(goodsubjects) * 100;
        scatterdata_y(:,3) = RawInCongWPAcc(goodsubjects) * 100;
        scatterdata_y(:,4) = nanmean(AccCNLPInCong(goodsubjects, :)') * 100;
        scatterdata_y(:,5) = nanmean(AccCNWPInCong(goodsubjects, :)') * 100;             

    end

    hold on;
    bar(x, bardata_mean, 'FaceColor', barColor);
    errorbar(x, bardata_mean, bardata_sem, '.k');
    scatter(scatterdata_x(:), scatterdata_y(:), markerSize, scatterColor, 'LineWidth', 1);
    plot(scatterdata_x', scatterdata_y', 'Color', [scatterColor 0.2]);
    hold off;

    ylim([0 100]);
    xlim([min(x)-0.5 max(x)+0.5]);
    set(gca, 'XTick', x);
    set(gca, 'XTickLabels', xLabels);
    set(gca, 'FontSize',  fontSize_ylabel);
    ylabel('Accuracy (%)','FontSize', fontSize_ylabel, 'FontName', fontName, 'FontSize', fontSize_ylabel);
    if(plotID == 1)
        title('Congruent Stroop Stimuli');
    elseif(plotID == 2)
        title('Incongruent Stroop Stimuli');
    end

end

%% TOWNSEND & WENGER ANALYSIS
plotSEM = 1;

% PEFFORM ANALYSIS

maxRT = 1.5;
RT_x = 0.00:0.01:maxRT;
nSubj = length(goodsubjects);

CNLP_C_all = nan(nSubj, length(RT_x));
CNWP_C_all = nan(nSubj, length(RT_x));

for subj_idx = 1:nSubj
    subject = goodsubjects(subj_idx)+100;
    
    % CNLP
    taskRT_all_AB = RT_data{subject}.CNLP;
    taskRT_all_A = RT_data{subject}.CN;
    taskRT_all_B = RT_data{subject}.LP;
    [CN_LP, CN_LP_1, min_CN_LP, CNLP_C]  = computeTownsendWengerCDF(taskRT_all_AB, taskRT_all_A, taskRT_all_B, RT_x);
    CNLP_C_all(subj_idx, :) = CNLP_C;
    
    % CNWP
    taskRT_all_AB = RT_data{subject}.CNWP;
    taskRT_all_A = RT_data{subject}.CN;
    taskRT_all_B = RT_data{subject}.WP;
    [CN_WP, CN_WP_1, min_CN_WP, CNWP_C]  = computeTownsendWengerCDF(taskRT_all_AB, taskRT_all_A, taskRT_all_B, RT_x);
    CNWP_C_all(subj_idx, :) = CNWP_C;
    
end

% remove invalid values
CNLP_C_all(CNLP_C_all == -Inf) = nan;
CNWP_C_all(CNWP_C_all == -Inf) = nan;
CNLP_C_all(CNLP_C_all == Inf) = nan;
CNWP_C_all(CNWP_C_all == Inf) = nan;

% plot
close all;
fig2 = figure(2);
set(fig2, 'Position', [100 100 400 200]);
plotSettings;

CNLP_mean = nanmean(CNLP_C_all);
CNLP_sem = nanstd(CNLP_C_all)./sqrt(sum(~isnan(CNLP_C_all)));
CNLP_sem(CNLP_sem == 0) = nan;

CNWP_mean = nanmean(CNWP_C_all);
CNWP_sem = nanstd(CNWP_C_all)./sqrt(sum(~isnan(CNWP_C_all)));
CNWP_sem(CNWP_sem == 0) = nan;

hold on;
if(plotSEM)
%     errorbar(RT_x, CNLP_mean, CNLP_sem, '-', 'LineWidth', 1, 'Color', [0 0 0 0.1]);
%     errorbar(RT_x, CNWP_mean, CNWP_sem, '-', 'LineWidth', 1, 'Color', [0 0 0 0.1]);
        s = shadedErrorBar(RT_x,CNLP_C_all,{@nanmean,@nansem},'lineprops', {'-', 'Color',  colors(cContrast2,:)});
        s.patch.FaceColor = colors(cContrast2,:);
        set(s.edge,'LineWidth',1,'LineStyle','-')
        
        s = shadedErrorBar(RT_x,CNWP_C_all,{@nanmean,@nansem},'lineprops', {'-', 'Color',  colors(cContrast3,:)});
        s.patch.FaceColor = colors(cContrast3,:);
        set(s.patch, 'EdgeColor', colors(cContrast3,:));
        set(s.edge,'LineWidth',1,'LineStyle','-')
end
p1 = plot(RT_x, CNLP_mean, '-k', 'LineWidth', 3, 'Color', colors(cContrast2,:));
p2 = plot(RT_x, CNWP_mean, '-k', 'LineWidth', 3, 'Color', colors(cContrast3,:));
hold off;

xlim([0, max(RT_x)]);
ylabel('Capacity C','FontSize', fontSize_ylabel, 'FontName', fontName, 'FontSize', fontSize_ylabel);
xlabel('Time t in Seconds','FontSize', fontSize_ylabel, 'FontName', fontName, 'FontSize', fontSize_ylabel);
leg = legend([p1 p2], 'CN+LP', 'CN+WP', 'location', 'northwest');
set(leg, 'FontSize', fontSize_ylabel); 
set(gca, 'FontSize', fontSize_gca)


%% TOWNSEND & ALTIERI ANALYSIS
plotSEM = 0;
plotIncorrect = 1;

% PEFFORM ANALYSIS

maxRT = 1.5;
RT_x = 0.01:0.01:maxRT;
nSubj = length(goodsubjects);

CNLP_A_AND_CF = nan(nSubj, length(RT_x));
CNWP_A_AND_CF = nan(nSubj, length(RT_x));
CNLP_A_AND_IF = nan(nSubj, length(RT_x));
CNWP_A_AND_IF = nan(nSubj, length(RT_x));

for subj_idx = 1:nSubj
    subject = goodsubjects(subj_idx)+100;
    
    % CNLP CORRECT AND FAST
    taskRT_correct_AB = RT_data{subject}.CNLP;
    taskRT_correct_A = RT_data{subject}.CN;
    taskRT_correct_B = RT_data{subject}.LP;
    taskRT_incorrect_AB = RT_data_incorrect{subject}.CNLP;
    taskRT_incorrect_A = RT_data_incorrect{subject}.CN;
    taskRT_incorrect_B =  RT_data_incorrect{subject}.LP;
    
    taskRT_correct_AB(isnan(taskRT_correct_AB)) = [];
    taskRT_correct_A(isnan(taskRT_correct_A)) = [];
    taskRT_correct_B(isnan(taskRT_correct_B)) = [];
    taskRT_incorrect_AB(isnan(taskRT_incorrect_AB)) = [];
    taskRT_incorrect_A(isnan(taskRT_incorrect_A)) = [];
    taskRT_incorrect_B(isnan(taskRT_incorrect_B)) = [];
    
    [A_AND_CF, A_AND_IF]  = computeTownsendAltieriCapacity(taskRT_correct_AB, taskRT_correct_A, taskRT_correct_B, taskRT_incorrect_AB, taskRT_incorrect_A, taskRT_incorrect_B, RT_x);
    CNLP_A_AND_CF(subj_idx, :) = A_AND_CF;
    CNLP_A_AND_IF(subj_idx, :) = A_AND_IF;
    
    % CNWP CORRECT AND FAST
    taskRT_correct_AB = RT_data{subject}.CNWP;
    taskRT_correct_A = RT_data{subject}.CN;
    taskRT_correct_B = RT_data{subject}.WP;
    taskRT_incorrect_AB = RT_data_incorrect{subject}.CNWP;
    taskRT_incorrect_A = RT_data_incorrect{subject}.CN;
    taskRT_incorrect_B =  RT_data_incorrect{subject}.WP;
    
    % remove nans
    taskRT_correct_AB(isnan(taskRT_correct_AB)) = [];
    taskRT_correct_A(isnan(taskRT_correct_A)) = [];
    taskRT_correct_B(isnan(taskRT_correct_B)) = [];
    taskRT_incorrect_AB(isnan(taskRT_incorrect_AB)) = [];
    taskRT_incorrect_A(isnan(taskRT_incorrect_A)) = [];
    taskRT_incorrect_B(isnan(taskRT_incorrect_B)) = [];
    
    [A_AND_CF, A_AND_IF]  = computeTownsendAltieriCapacity(taskRT_correct_AB, taskRT_correct_A, taskRT_correct_B, taskRT_incorrect_AB, taskRT_incorrect_A, taskRT_incorrect_B, RT_x);
    CNWP_A_AND_CF(subj_idx, :) = A_AND_CF;
    CNWP_A_AND_IF(subj_idx, :) = A_AND_IF;
    
end

% remove invalid values
CNLP_A_AND_CF(CNLP_A_AND_CF == -Inf) = nan;
CNWP_A_AND_CF(CNWP_A_AND_CF == -Inf) = nan;
CNLP_A_AND_IF(CNLP_A_AND_IF == -Inf) = nan;
CNWP_A_AND_IF(CNWP_A_AND_IF == -Inf) = nan;

CNLP_A_AND_CF(CNLP_A_AND_CF == Inf) = nan;
CNWP_A_AND_CF(CNWP_A_AND_CF == Inf) = nan;
CNLP_A_AND_IF(CNLP_A_AND_IF == Inf) = nan;
CNWP_A_AND_IF(CNWP_A_AND_IF == Inf) = nan;

% PLOT CORRECT AND FAST
fig2 = figure(2);
if(plotIncorrect)
    set(fig2, 'Position', [100 100 700 200]);
    subplot(1,2,1);
else
    set(fig2, 'Position', [100 100 300 200]);
end
plotSettings;

CNLP_A_AND_CF_mean = nanmean(CNLP_A_AND_CF);
CNLP_A_AND_CF_sem = nanstd(CNLP_A_AND_CF)./sqrt(sum(~isnan(CNLP_A_AND_CF)));
CNLP_A_AND_CF_sem(CNLP_A_AND_CF_sem == 0) = nan;

CNWP_A_AND_CF_mean = nanmean(CNWP_A_AND_CF);
CNWP_A_AND_CF_sem = nanstd(CNWP_A_AND_CF)./sqrt(sum(~isnan(CNWP_A_AND_CF)));
CNWP_A_AND_CF_sem(CNWP_A_AND_CF_sem == 0) = nan;

hold on;
plot(RT_x, CNLP_A_AND_CF_mean, '--k', 'LineWidth', 3);
plot(RT_x, CNWP_A_AND_CF_mean, '-k', 'LineWidth', 3);

if(plotSEM)
    errorbar(RT_x, CNLP_A_AND_CF_mean, CNLP_A_AND_CF_sem, '-k', 'LineWidth', 1);
    errorbar(RT_x, CNWP_A_AND_CF_mean, CNWP_A_AND_CF_sem, '-k', 'LineWidth', 1);
end
hold off;

title('Fast and Correct','FontSize', fontSize_title, 'FontName', fontName, 'FontSize', fontSize_title);
ylabel('A(t)','FontSize', fontSize_ylabel, 'FontName', fontName, 'FontSize', fontSize_ylabel);
xlabel('Time t in Seconds','FontSize', fontSize_ylabel, 'FontName', fontName, 'FontSize', fontSize_ylabel);
leg = legend('CN+LP', 'CN+WP');
set(leg, 'FontSize', fontSize_ylabel); 

if(plotIncorrect)
% PLOT INCORRECT AND FAST
subplot(1,2,2);

CNLP_A_AND_IF_mean = nanmean(CNLP_A_AND_IF);
CNLP_A_AND_IF_sem = nanstd(CNLP_A_AND_IF)./sqrt(sum(~isnan(CNLP_A_AND_IF)));
CNLP_A_AND_IF_sem(CNLP_A_AND_IF_sem == 0) = nan;

CNWP_A_AND_IF_mean = nanmean(CNWP_A_AND_IF);
CNWP_A_AND_IF_sem = nanstd(CNWP_A_AND_IF)./sqrt(sum(~isnan(CNWP_A_AND_IF)));
CNWP_A_AND_IF_sem(CNWP_A_AND_IF_sem == 0) = nan;

hold on;
plot(RT_x, CNLP_A_AND_IF_mean, '--k', 'LineWidth', 3);
plot(RT_x, CNWP_A_AND_IF_mean, '-k', 'LineWidth', 3);

if(plotSEM)
    errorbar(RT_x, CNLP_A_AND_IF_mean, CNLP_A_AND_IF_sem, '-k', 'LineWidth', 1);
    errorbar(RT_x, CNWP_A_AND_IF_mean, CNWP_A_AND_IF_sem, '-k', 'LineWidth', 1);
end
hold off;

title('Fast and Incorrect','FontSize', fontSize_title, 'FontName', fontName, 'FontSize', fontSize_title);
ylabel('A(t)','FontSize', fontSize_ylabel, 'FontName', fontName, 'FontSize', fontSize_ylabel);
xlabel('Time t in Seconds','FontSize', fontSize_ylabel, 'FontName', fontName, 'FontSize', fontSize_ylabel);
leg = legend('CN+LP', 'CN+WP');
set(leg, 'FontSize', fontSize_ylabel, 'location', 'southwest'); 

end