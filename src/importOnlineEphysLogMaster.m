function ephysLog = importOnlineEphysLogMaster(dirs)
try
sessionSheet = GetGoogleSpreadsheet('1ivwAb6mOP0DCNM_x0owi1eP5aGump0zCu-dGorI4XjM');
ephysLog = cell2table(sessionSheet(2:end,:),'VariableNames',sessionSheet(1,:));

catch
    ephysLog = readtable([dirs.results '\ephysLog_w_PransData.xlsx']);
end
ephysLog = sortrows(ephysLog,"Session","ascend");
ephysLog = sortrows(ephysLog,"Date","ascend");
assignProbeLabels
end
