function toggleVisibility(checkbox, plotHandle)
    if checkbox.Value
        plotHandle.Visible = 'on';
    else
        plotHandle.Visible = 'off';
    end
end