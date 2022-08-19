function plot8traj_3Dgaussian_AH(data, smooth_window, d_legend, d_marker_loc, d_marker_size, d_marker_legend, l_color, l_opacity, l_width, p_color, p_size, p_freq, a_xlim, a_ylim, a_xlabel, a_ylabel, a_title, a_fsize)
% Kanha Batra's neural trajectory script, adaptation to other smoothing
% Aniek
%      
    % Purpose: plot 8 trajectories in the same figure using gaussian
    %          smoothening
    %
    % Inputs:
    % data : 1x6 cell array
    %     each element is a double array with dimensions 2xtimexcases 
    %     (denoted by T and C respectively)
    % smooth_window : Array of 2 numbers, frames used for smoothing
    %     Always put 0 as second to not use future smoothing
    % d_legend : 1x6 cell array
    %     each element is a string for the legend entry of the
    %     corresponding trajectory
    % d_marker_loc : 1x3 double array
    %     each element represents the onset of an event in the time series,
    %     eg baseline, cue onset and trial end; each value must lie between
    %     0 and T
    % d_marker_size : double greater than 0
    %     size of the event marker
    % d_marker_legend : 1x3 cell array
    %     each element is a string for the name of the event
    % l_color : 1x6 cell array
    %     each element is a 1x3 double array for the color of the
    %     corresponding trajectory; each element of the array is between 0
    %     and 255, i.e., the RGB values
    % l_opacity : double between 0 and 1
    %     opacity for the line plot
    % l_width : double greater than 0
    %     width of the line plot
    % p_color : 1x6 cell array
    %     each element is a 1x3 double array for the color of the
    %     corresponding time bin markers; each element of the array is
    %     between 0 and 255, i.e., the RGB values
    % p_size : double greater than 0
    %     size of time bin markers
    % p_freq : double between 0 and T
    %     frequency of markers, eg if p_freq = 5, every 5th bin will be
    %     marked
    %
    % Optional Input:
    % a_xlim : 1x2 double array
    %     x-axis limits
    % a_ylim : 1x2 double array
    %     y-axis limits
    % a_xlabel : string
    %     title of the x-axis units
    % a_ylabel : string
    %     title of the y-axis units
    % a_title : string
    %     title of the overall plot
    % a_fsize : double
    %     font size for the axes titles, axes labels, figure title and
    %     legend
    %
    % Output:
    % <none> : a new figure is constructed
    
    % Code:
    d1 = smoothdata(mean(data{1}, 3), 2, 'gaussian', smooth_window);
    d2 = smoothdata(mean(data{2}, 3), 2, 'gaussian', smooth_window);
    d3 = smoothdata(mean(data{3}, 3), 2, 'gaussian', smooth_window);
    d4 = smoothdata(mean(data{4}, 3), 2, 'gaussian', smooth_window);
    d5 = smoothdata(mean(data{5}, 3), 2, 'gaussian', smooth_window);
    d6 = smoothdata(mean(data{6}, 3), 2, 'gaussian', smooth_window);
    d7 = smoothdata(mean(data{7}, 3), 2, 'gaussian', smooth_window);
    d8 = smoothdata(mean(data{8}, 3), 2, 'gaussian', smooth_window);

    hold on;
    p1 = plot3(d1(1, :), d1(2, :), d1(3,:), 'DisplayName', d_legend{1});
    p1.Color(1: 3) = l_color{1}; p1.Color(4) = l_opacity; p1.LineWidth = l_width;
    p1.Marker = '.'; p1.MarkerFaceColor = p_color{1}; p1.MarkerEdgeColor = p_color{1}; p1.MarkerIndices = [1: p_freq: size(d1, 2)]; p1.MarkerSize = p_size;
    e11 = scatter3(d1(1, d_marker_loc(1)), d1(2, d_marker_loc(1)), d1(3, d_marker_loc(1)), d_marker_size, [0, 0, 0]/255, 'o', 'filled', 'HandleVisibility', 'off');
    e12 = scatter3(d1(1, d_marker_loc(2)), d1(2, d_marker_loc(2)), d1(3, d_marker_loc(2)), d_marker_size, [0, 0, 0]/255, '>', 'filled', 'HandleVisibility', 'off');
    e13 = scatter3(d1(1, d_marker_loc(3)), d1(2, d_marker_loc(3)), d1(3, d_marker_loc(3)), d_marker_size, [0, 0, 0]/255, 'square', 'filled', 'HandleVisibility', 'off');
    
    p2 = plot3(d2(1, :), d2(2, :),d2(3, :), 'DisplayName', d_legend{2});
    p2.Color(1: 3) = l_color{2}; p2.Color(4) = l_opacity; p2.LineWidth = l_width;
    p2.Marker = '.'; p2.MarkerFaceColor = p_color{2}; p2.MarkerEdgeColor = p_color{2}; p2.MarkerIndices = [1: p_freq: size(d2, 2)]; p2.MarkerSize = p_size;
    e21 = scatter3(d2(1, d_marker_loc(1)), d2(2, d_marker_loc(1)), d2(3, d_marker_loc(1)), d_marker_size, [0, 0, 0]/255, 'o', 'filled', 'HandleVisibility', 'off');
    e22 = scatter3(d2(1, d_marker_loc(2)), d2(2, d_marker_loc(2)), d2(3, d_marker_loc(2)), d_marker_size, [0, 0, 0]/255, '>', 'filled', 'HandleVisibility', 'off');
    e23 = scatter3(d2(1, d_marker_loc(3)), d2(2, d_marker_loc(3)), d2(3, d_marker_loc(3)), d_marker_size, [0, 0, 0]/255, 'square', 'filled', 'HandleVisibility', 'off');
   
    p3 = plot3(d3(1, :), d3(2, :),d3(3, :), 'DisplayName', d_legend{3});
    p3.Color(1: 3) = l_color{3}; p3.Color(4) = l_opacity; p3.LineWidth = l_width;
    p3.Marker = '.'; p3.MarkerFaceColor = p_color{3}; p3.MarkerEdgeColor = p_color{3}; p3.MarkerIndices = [1: p_freq: size(d3, 2)]; p3.MarkerSize = p_size;
    e31 = scatter3(d3(1, d_marker_loc(1)), d3(2, d_marker_loc(1)), d3(3, d_marker_loc(1)), d_marker_size, [0, 0, 0]/255, 'o', 'filled', 'HandleVisibility', 'off');
    e32 = scatter3(d3(1, d_marker_loc(2)), d3(2, d_marker_loc(2)), d3(3, d_marker_loc(2)), d_marker_size, [0, 0, 0]/255, '>', 'filled', 'HandleVisibility', 'off');
    e33 = scatter3(d3(1, d_marker_loc(3)), d3(2, d_marker_loc(3)), d3(3, d_marker_loc(3)), d_marker_size, [0, 0, 0]/255, 'square', 'filled', 'HandleVisibility', 'off');
   
    p4 = plot3(d4(1, :), d4(2, :),d4(3, :), 'DisplayName', d_legend{4});
    p4.Color(1: 3) = l_color{4}; p4.Color(4) = l_opacity; p4.LineWidth = l_width;
    p4.Marker = '.'; p4.MarkerFaceColor = p_color{4}; p4.MarkerEdgeColor = p_color{4}; p4.MarkerIndices = [1: p_freq: size(d4, 2)]; p4.MarkerSize = p_size;
    e41 = scatter3(d4(1, d_marker_loc(1)), d4(2, d_marker_loc(1)), d4(3, d_marker_loc(1)), d_marker_size, [0, 0, 0]/255, 'o', 'filled', 'HandleVisibility', 'off');
    e42 = scatter3(d4(1, d_marker_loc(2)), d4(2, d_marker_loc(2)), d4(3, d_marker_loc(2)), d_marker_size, [0, 0, 0]/255, '>', 'filled', 'HandleVisibility', 'off');
    e43 = scatter3(d4(1, d_marker_loc(3)), d4(2, d_marker_loc(3)), d4(3, d_marker_loc(3)), d_marker_size, [0, 0, 0]/255, 'square', 'filled', 'HandleVisibility', 'off');
    
    p5 = plot3(d5(1, :), d5(2, :), d5(3, :), 'DisplayName', d_legend{5});
    p5.Color(1: 3) = l_color{5}; p5.Color(4) = l_opacity; p5.LineWidth = l_width;
    p5.Marker = '.'; p5.MarkerFaceColor = p_color{5}; p5.MarkerEdgeColor = p_color{5}; p5.MarkerIndices = [1: p_freq: size(d5, 2)]; p5.MarkerSize = p_size;
    e51 = scatter3(d5(1, d_marker_loc(1)), d5(2, d_marker_loc(1)), d5(3, d_marker_loc(1)), d_marker_size, [0, 0, 0]/255, 'o', 'filled', 'HandleVisibility', 'off');
    e52 = scatter3(d5(1, d_marker_loc(2)), d5(2, d_marker_loc(2)), d5(3, d_marker_loc(2)), d_marker_size, [0, 0, 0]/255, '>', 'filled', 'HandleVisibility', 'off');
    e53 = scatter3(d5(1, d_marker_loc(3)), d5(2, d_marker_loc(3)), d5(3, d_marker_loc(3)), d_marker_size, [0, 0, 0]/255, 'square', 'filled', 'HandleVisibility', 'off');
    
    p6 = plot3(d6(1, :), d6(2, :), d6(3, :), 'DisplayName', d_legend{6});
    p6.Color(1: 3) = l_color{6}; p6.Color(4) = l_opacity; p6.LineWidth = l_width;
    p6.Marker = '.'; p6.MarkerFaceColor = p_color{6}; p6.MarkerEdgeColor = p_color{6}; p6.MarkerIndices = [1: p_freq: size(d6, 2)]; p6.MarkerSize = p_size;
    e61 = scatter3(d6(1, d_marker_loc(1)), d6(2, d_marker_loc(1)), d6(3, d_marker_loc(1)), d_marker_size, [0, 0, 0]/255, 'o', 'filled', 'HandleVisibility', 'off');
    e62 = scatter3(d6(1, d_marker_loc(2)), d6(2, d_marker_loc(2)), d6(3, d_marker_loc(2)), d_marker_size, [0, 0, 0]/255, '>', 'filled', 'HandleVisibility', 'off');
    e63 = scatter3(d6(1, d_marker_loc(3)), d6(2, d_marker_loc(3)), d6(3, d_marker_loc(3)), d_marker_size, [0, 0, 0]/255, 'square', 'filled', 'HandleVisibility', 'off');
    
    
    p7 = plot3(d7(1, :), d7(2, :), d7(3, :), 'DisplayName', d_legend{7});
    p7.Color(1: 3) = l_color{7}; p7.Color(4) = l_opacity; p7.LineWidth = l_width;
    p7.Marker = '.'; p7.MarkerFaceColor = p_color{7}; p7.MarkerEdgeColor = p_color{7}; p7.MarkerIndices = [1: p_freq: size(d7, 2)]; p7.MarkerSize = p_size;
    e71 = scatter3(d7(1, d_marker_loc(1)), d7(2, d_marker_loc(1)), d7(3, d_marker_loc(1)), d_marker_size, [0, 0, 0]/255, 'o', 'filled', 'HandleVisibility', 'off');
    e72 = scatter3(d7(1, d_marker_loc(2)), d7(2, d_marker_loc(2)), d7(3, d_marker_loc(2)), d_marker_size, [0, 0, 0]/255, '>', 'filled', 'HandleVisibility', 'off');
    e73 = scatter3(d7(1, d_marker_loc(3)), d7(2, d_marker_loc(3)), d7(3, d_marker_loc(3)), d_marker_size, [0, 0, 0]/255, 'square', 'filled', 'HandleVisibility', 'off');
    
    
    p8 = plot3(d8(1, :), d8(2, :), d8(3, :), 'DisplayName', d_legend{8});
    p8.Color(1: 3) = l_color{8}; p8.Color(4) = l_opacity; p8.LineWidth = l_width;
    p8.Marker = '.'; p8.MarkerFaceColor = p_color{8}; p8.MarkerEdgeColor = p_color{8}; p8.MarkerIndices = [1: p_freq: size(d8, 2)]; p8.MarkerSize = p_size;
    e81 = scatter3(d8(1, d_marker_loc(1)), d8(2, d_marker_loc(1)), d8(3, d_marker_loc(1)), d_marker_size, [0, 0, 0]/255, 'o', 'filled', 'HandleVisibility', 'off');
    e82 = scatter3(d8(1, d_marker_loc(2)), d8(2, d_marker_loc(2)), d8(3, d_marker_loc(2)), d_marker_size, [0, 0, 0]/255, '>', 'filled', 'HandleVisibility', 'off');
    e83 = scatter3(d8(1, d_marker_loc(3)), d8(2, d_marker_loc(3)), d8(3, d_marker_loc(3)), d_marker_size, [0, 0, 0]/255, 'square', 'filled', 'HandleVisibility', 'off');
    
       
    e1 = scatter3(NaN, NaN, NaN, d_marker_size, [0, 0, 0]/255, 'o', 'filled', 'Displayname', d_marker_legend{1});
    e2 = scatter3(NaN, NaN, NaN, d_marker_size, [0, 0, 0]/255, '>', 'filled', 'Displayname', d_marker_legend{2});
    e3 = scatter3(NaN, NaN, NaN, d_marker_size, [0, 0, 0]/255, 'square', 'filled', 'Displayname', d_marker_legend{3});
    
    hold off;
    
    if nargin > 12
        xlim(a_xlim);
        ylim(a_ylim);
        xlabel(a_xlabel);
        ylabel(a_ylabel);
        title(a_title);
        set(findall(gcf,'-property','FontSize'),'FontSize',a_fsize)
    end
end
