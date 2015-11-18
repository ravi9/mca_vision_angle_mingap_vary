function OutputImage = DrawSelectedLines (Lines, Layout, selected)

% Service function to draw the selected Lines in a Layout as an image.

N = length(Layout.indices);

for i=1:N
    L = Lines(Layout.indices(i)).length;
    for j=1:L
        OutputImage(Lines(Layout.indices(i)).y(j), Lines(Layout.indices(i)).x(j)) = 255*selected(i);
    end;
    %% highlight the end points
    for (ii=-1:1)
        for (jj=-1:1)
            x = max(1, Lines(Layout.indices(i)).y(1)+ii);
            y = max(1, Lines(Layout.indices(i)).x(1)+jj);
            OutputImage(x, y) = 255*selected(i);
            x = max(1, Lines(Layout.indices(i)).y(L)+ii);
            y = max(1, Lines(Layout.indices(i)).x(L)+jj);
            OutputImage(x, y) = 255*selected(i);
        end
    end

end;
OutputImage = 255 - OutputImage;