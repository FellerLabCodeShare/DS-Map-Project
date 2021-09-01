function silValue = silTest(prefDirections)


%Remove NaN
prefDirections(isnan(prefDirections)) = [];

%Specify how many clusters to analyze, from 2 to numClust
numClust = 20;

silValue = NaN(1,numClust);
if length(prefDirections) > numClust

    %Turn prefDir into cartesian coordinates, assume vector length is 1 for all
    [x,y] = pol2cart(prefDirections,ones(length(prefDirections),1));

    %Initialize silValue
    silValue = NaN(length(x),numClust);



    %Run k means, for now we will do clusters = 2
    for j = 2:numClust
        [idx,C] = kmeans([x y],j);


        for i = 1:length(x)
            %What cluster are we in
            clustNum = idx(i);
            
            b = NaN(j,1);
            
            for m = 1:j
                b(m) = mean(sqrt((x(idx == m)-x(i)).^2 + (y(idx == m)-y(i)).^2));
            end
            
            a = b(clustNum);
            b(clustNum) = [];
            bmin = min(b);
            
            silValue(i,j) = (bmin-a)/max([a,bmin]);

        end

    end

end

end



%Below is the old script which does the calculation in an incorrect manner
%
% %Remove NaN
% prefDirections(isnan(prefDirections)) = [];
% 
% %Specify how many clusters to analyze, from 2 to numClust
% numClust = 8;
% 
% silValue = NaN(1,numClust);
% if length(prefDirections) > numClust
% 
%     %Turn prefDir into cartesian coordinates, assume vector length is 1 for all
%     [x,y] = pol2cart(prefDirections,ones(length(prefDirections),1));
% 
%     %Initialize silValue
%     silValue = NaN(length(x),numClust);
% 
% 
% 
%     %Run k means, for now we will do clusters = 2
%     for j = 2:numClust
%         [idx,C] = kmeans([x y],j);
% 
% 
%         for i = 1:length(x)
%             %What cluster are we in
%             clustNum = idx(i);
% 
%             %What's the nearest cluster to current point
%             tempX = x(idx ~= clustNum);
%             tempY = y(idx ~= clustNum);
%             tempInd = idx(idx ~= clustNum);
%             nearestInd = find(sqrt((tempX-x(i)).^2 + (tempY-y(i)).^2)-min(sqrt((tempX-x(i)).^2 + (tempY-y(i)).^2)) == 0,1,'first');
%             nearestClust = tempInd(nearestInd);
% 
%             a = mean(sqrt((x(idx == clustNum)-x(i)).^2 + (y(idx == clustNum)-y(i)).^2));
%             b = mean(sqrt((x(idx == nearestClust)-x(i)).^2 + (y(idx == nearestClust)-y(i)).^2));
%             silValue(i,j) = (b-a)/max([a,b]);
% 
%         end
% 
%     end
% 
% end
% 
% end