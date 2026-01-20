function logInfo = util_getLogInfo(neuralFilename,ephysLog)
    log_i = find(strcmp(ephysLog.Session,neuralFilename));
    
    ap_stereo = str2num(ephysLog.AP_Stereotaxic{log_i});
    ml_stereo = str2num(ephysLog.ML_Stereotaxic{log_i});
    electrode_depth = str2num(ephysLog.ElectrodeSettleDepth{log_i});
    acs_ch = str2num(ephysLog.accPranBank{log_i});
    electrode_spc = str2num(ephysLog.ElectrodeSpacing{log_i});
    ctx_top_ch = str2num(ephysLog.ctxPran{log_i});
    probeName = {ephysLog.probeLabel{log_i}};
    
    if isempty(acs_ch)
        acs_ch = NaN;
    end
    if isempty(ctx_top_ch)
        ctx_top_ch = NaN;
    end
    
    neuralFilename = {neuralFilename};
    logInfo = table(neuralFilename,log_i,ap_stereo,ml_stereo,...
        electrode_depth,acs_ch,electrode_spc,ctx_top_ch,probeName);
end
