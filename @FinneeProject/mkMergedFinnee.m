function obj = mkMergedFinnee(obj, dataset2beMerged, partitions, Tag4files2merged, Tag4newfiles, isAligned)

IdF = find(strcmp(obj.Summary.FinneeType, Tag4files2merged));
FolderOut = fullfile(obj.OutputFolder, 'MergedFinnee');
ParallelMe = true;


c = cvpartition(numel(IdF),'KFold', partitions);
for ii = 1:partitions
    IdX = find(c.test(ii));
    
    for jj = 1:numel(IdX)
        FinnesIn{ii, jj} = fullfile(obj.Summary.FolderID{IdX(jj)}, [obj.Summary.FileID{IdX(jj)}, '.fin']);
    
    end

    if isAligned
        plignMe{ii} = {obj.Summary.AlignMe{IdX}};
    end
end

FileID = {}; StartTime = {}; EndTime = {};
FolderID = {}; FinneeType = {}; Dataset4Quali = []; Dataset4Quant = [];
dts4Ql = 1; dts4Qt = 1; AlignMe = {};


for ii = 1:partitions
    dirs = {FinnesIn{ii, :}}';
    if isempty(dirs{end})
        dirs(end) = [];

    end
    StartTime{ii, 1} = datetime;
    FinneeType{ii, 1} = Tag4newfiles;
    fileID = [Tag4newfiles, '_', num2str(ii)];
    disp([fileID, ' ... STARTING ...']),
    FileID{ii, 1} = fileID; AlignMe{ii, 1} = {};



    if isAligned
        myFinnee = Finnee('Multiple', 'overwrite', ...
            'TgtDatasets', dataset2beMerged, ...
            'Alignment', plignMe{ii}, ...
            'FinneesIn', dirs, ...
            'folderOut', FolderOut, ...
            'fileID', fileID, ...
            'spikes', 0);

    else
        myFinnee = Finnee('Multiple', 'overwrite', ...
            'TgtDatasets', dataset2beMerged, ...
            'Alignment', [], ...
            'FinneesIn', dirs, ...
            'folderOut', FolderOut, ...
            'fileID', fileID, ...
            'spikes', 0);

    end
    
    dts4Ql = 3;
    disp([fileID, ' ... DONE!']),
    FolderID{ii, 1} = FolderOut;
    EndTime{ii, 1} = datetime;
    Dataset4Quali(ii, 1) = dts4Ql;
    Dataset4Quant(ii, 1) = dts4Qt;

end

obj.Summary = [obj.Summary; table(FileID, FolderID, FinneeType, StartTime, EndTime, Dataset4Quali, Dataset4Quant, AlignMe)];
myProject = obj; %#ok<*NASGU>
save(fullfile(obj.Path2Project, 'myProject.mat'), 'myProject')

end