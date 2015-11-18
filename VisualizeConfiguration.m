function y = VisualizeConfiguration (X, Y, Lines)

% Stand alone function used to debug and visualize as images the Lines corresponding
% to the nanomagnets at (X, Y). The user can click on a nanomagnet and the
% corresponding Line is shown as an image.
%

close all; figure;
selected = [];
while (1 == 1)
    subplot(1,2,1);
    plot(X, Y, 'o'); hold on;
    if (size(selected) > 0) plot(X(selected), Y(selected), 'ro');  end;
    
    for i=1:length(Lines)
        for (j=1:Lines(i).length)
            if (ismember(i, selected))
                OutputImage(Lines(i).y(j), Lines(i).x(j)) = 255;
            else
                OutputImage(Lines(i).y(j), Lines(i).x(j)) = 120;
            end;
        end;
    end;
    subplot(1,2,2); imagesc(OutputImage);

    [x y] = ginput(1);
    d = (X-x).^2 + (Y-y).^2;
    [mind i] = min(d);
    selected = [selected i];
end;