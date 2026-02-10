function [outliers] = medAbsDeviation(data)

% ===============================================
% Automatically remove artefacts occurring as sharp peaks by computing
% a robust z-scores (based on the median and the median absolute
% devition) and replacing and interpolating all values that fall beyond
% 4
med_abs_dev = 1.4826 * median(abs(data - median(data, 'omitnan')), 'omitnan');  %  the constant 1.4826 assumes normality of the data
med_abs_dev_scores = (data - median(data, 'omitnan')) / med_abs_dev;    
OUTL_THRESHOLD = 3;
% numel(data(med_abs_dev_scores > OUTL_THRESHOLD))
% figure; 
% plot(data); hold on
% plot(find(abs(med_abs_dev_scores) > OUTL_THRESHOLD), data(abs(med_abs_dev_scores) > OUTL_THRESHOLD), 'x', 'Color', 'r')
outliers = abs(med_abs_dev_scores) > OUTL_THRESHOLD;

