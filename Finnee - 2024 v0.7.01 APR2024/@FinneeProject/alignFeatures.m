function [obj, pSmu] = alignFeatures(obj, Id2PeakList)
minRpl = 50;

myMPL = obj.FeaturesLists{Id2PeakList}.Original;
orderFit = 2;
load(myMPL.Name)

MyTmNor = table();
MyTmNor.ID = myMrgPeakList.Summary.idMF;

% initialisation
for ii = 1:numel(myMrgPeakList.AdditionalInformation.ListFinnee.Name)
    MyTmNor.(myMrgPeakList.AdditionalInformation.ListFinnee.Name{ii}) = nan(numel(MyTmNor.ID), 1);

end

% Fill the table
for ii = 1:numel(MyTmNor.ID)
    myMF = myMrgPeakList.MergedFeatures{MyTmNor.ID(ii)};

    StackMe.chrom_apex = myMF.chrom_Time_IA;
    StackMe.chrom_Iapex = myMF.chrom_intensity_apex;
    StackMe.Id2Dts = myMF.Id2Dts;
    vector = StackMe.Id2Dts;
    ListVec = unique(StackMe.Id2Dts);

    for jj = 1:numel(ListVec)
        idx = find(StackMe.Id2Dts == ListVec(jj));
        if numel(idx) == 1
            MyTmNor(ii, ListVec(jj)+1) = {StackMe.chrom_apex(idx)};

        else
            [~, MX] = max(StackMe.chrom_Iapex(idx));
            MyTmNor(ii, ListVec(jj)+1) = {StackMe.chrom_apex(idx(MX))};

        end
    end


end

MyTmNor(sum(~isnan(table2array(MyTmNor(:, 2:end))), 2) < minRpl, :) = [];
MyTmNor.AvgTm = mean(table2array(MyTmNor(:, 2:end)), 2, 'omitnan');

myMrgPeakList.AdditionalInformation.myAlignment.Data4 = MyTmNor;
for ii = 1:numel(myMrgPeakList.AdditionalInformation.ListFinnee.Name)
    clear XY_
    XY = [MyTmNor.(ii+1), MyTmNor.AvgTm-MyTmNor.(ii+1)];
    XY(any(isnan(XY), 2), :) = [];
    [p,S,mu] = polyfit(XY(:, 1), XY(:, 2), orderFit);
    [XY_(:,1), XY_(:,2)] = polyval(p, XY(:, 1), S, mu);
    XY_(:, 3) = abs(XY(:, 2) - XY_(:, 1)) > XY_(:, 2);
    XY(XY_(:, 3) == 1, :) = [];
    [p,S,mu] = polyfit(XY(:, 1), XY(:, 2), orderFit);

    myMrgPeakList.AdditionalInformation.myAlignment.pSmu{ii} = {p, S, mu};


end
pSmu = myMrgPeakList.AdditionalInformation.myAlignment.pSmu;
save(myMrgPeakList.Name, 'myMrgPeakList')

for ii = 1:numel(myMrgPeakList.AdditionalInformation.ListFinnee.Name)
    name = myMrgPeakList.AdditionalInformation.ListFinnee.Name{ii};
    IdF = find(strcmp(obj.Summary.FileID, name));
    obj.Summary.AlignMe{IdF} = pSmu{ii};

end

myProject = obj; %#ok<*NASGU>
save(fullfile(obj.Path2Project, 'myProject.mat'), 'myProject')

end