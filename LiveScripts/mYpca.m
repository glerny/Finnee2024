Summary = myMrgPeakList.TargetedAnalysis{2}.Summary;
[~, IdX] = intersect(Area.ID, Summary.ID);
Ar = Area.Values(IdX, :);
 [wcoeff,score,latent,tsquared,explained] = pca(Ar','VariableWeights','variance');

 figure
 hold on
 for ii = 1:16
      scatter3(score(myMap.Code == ii, 1), score(myMap.Code == ii, 2), score(myMap.Code == ii, 3));
 end
