function plotOptimization(x,y,labels,labeltext,xlab,ylab)

figure;
xylab=[x,y,labels];
FigureLabels=cell(numel(unique(labels)),1);
labelList=unique(labels);
for nl=1:numel(labelList)
    lb=labelList(nl);
    xy=xylab(xylab(:,3)==lb,[1 2]);
    FigureLabels{nl,1}=[labeltext '=' num2str(lb)];
    %scatter(xy(:,1),xy(:,2),'filled');
    plot(xy(:,1),xy(:,2),'-o');
    hold on
end
FigureLabels=string(FigureLabels);
legend(FigureLabels);
xlabel(xlab);
ylabel(ylab);
maxx=x(y(:)==max(y(:)));
maxy=y(y(:)==max(y(:)));
maxl=labels(y(:)==max(y(:)));
xposition=find(x==maxx);
if xposition >= numel(x)/2
    txt = ['(' num2str(maxx) ',' num2str(maxy) ') at ' labeltext '=' num2str(maxl) ' \rightarrow '];
    text(maxx,maxy,txt,'HorizontalAlignment','right');
else
    txt = [' \leftarrow(' num2str(maxx) ',' num2str(maxy) ') at ' labeltext '=' num2str(maxl)];
    text(maxx,maxy,txt,'HorizontalAlignment','left');
end
end