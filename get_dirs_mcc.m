function dirs = get_dirs_mcc(user)
user = user(1:end-1);
switch user
    case 'rri-martinez3'
        dirs.root = 'E:\conflict\analyses';
        dirs.toolbox = 'E:\conflict\analyses\toolbox';
        dirs.dajo_toolbox = 'E:\conflict\analyses\2022-dajo-toolbox-main';
        dirs.data = 'E:\conflict\data';
        dirs.results = 'E:\conflict\results';
        
    case 'DataCruncher'
        dirs.root = 'L:\conflict\analyses';
        dirs.toolbox = 'L:\conflict\analyses\toolbox';
        dirs.dajo_toolbox = 'L:\conflict\analyses\2022-dajo-toolbox-main';
        dirs.data = 'L:\conflict\data';
        dirs.results = 'L:\conflict\results';

end

addpath(genpath(dirs.root));
addpath(genpath(dirs.toolbox));
addpath(genpath(dirs.dajo_toolbox));
addpath(genpath(dirs.data));

end

