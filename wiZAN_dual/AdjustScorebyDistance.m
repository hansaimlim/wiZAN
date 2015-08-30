function Dist = AdjustScorebyDistance(P, center) %P: raw score matrix; center: hypothesized center value
Dist=(center-P).^2;
end
