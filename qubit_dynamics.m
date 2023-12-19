% MATLAB Final
% Julian Barossi


function qubitSummationModel()
    % Figure for the sliders and the plot
    f = figure('Name', 'Qubit States on Bloch Sphere', 'NumberTitle', 'off', 'Position', [100, 100, 800, 600]);

    % First  sliders for the first spin vector
    % Slider for phi1 
    sliderPhi1 = uicontrol('Parent', f, 'Style', 'slider', 'Position', [100, 120, 300, 20],'value', 0, 'min', 0, 'max', 2*pi, 'SliderStep', [1/300 1/30]);
    txtPhi1 = uicontrol('Style', 'text', 'Position', [100, 140, 300, 20],'String', 'First Spin - Phi: 0 rad');

    % Slider for theta1 (polar angle)
    sliderTheta1 = uicontrol('Parent', f, 'Style', 'slider', 'Position', [100, 90, 300, 20],'value', 0, 'min', 0, 'max', pi, 'SliderStep', [1/300 1/30]);
    txtTheta1 = uicontrol('Style', 'text', 'Position', [100, 110, 300, 20],'String', 'First Spin - Theta: 0 rad');

    % Second sliders for the second spin vector
    % Slider for phi2 (azimuthal angle)
    sliderPhi2 = uicontrol('Parent', f, 'Style', 'slider', 'Position', [400, 120, 300, 20], ...
                           'value', 0, 'min', 0, 'max', 2*pi, 'SliderStep', [1/300 1/30]);
    txtPhi2 = uicontrol('Style', 'text', 'Position', [400, 140, 300, 20], ...
                        'String', 'Second Spin - Phi: 0 rad');

    % Slider for theta2 (polar angle)
    sliderTheta2 = uicontrol('Parent', f, 'Style', 'slider', 'Position', [400, 90, 300, 20], ...
                             'value', 0, 'min', 0, 'max', pi, 'SliderStep', [1/300 1/30]);
    txtTheta2 = uicontrol('Style', 'text', 'Position', [400, 110, 300, 20], ...
                          'String', 'Second Spin - Theta: 0 rad');

    % Text for displaying the sum vector's phi and theta
    txtSumPhiTheta = uicontrol('Style', 'text', 'Position', [550, 30, 250, 40], ...
                               'String', 'Sum Vector - Phi: 0 rad, Theta: 0°');

 

    % Initial plot setup
    update_plot();

    % Callbacks for the sliders
    sliderPhi1.Callback = @(es, ed) update_plot();
    sliderTheta1.Callback = @(es, ed) update_plot();
    sliderPhi2.Callback = @(es, ed) update_plot();
    sliderTheta2.Callback = @(es, ed) update_plot();

    function update_plot()
        % Clear 
        cla;

        % slider value
        phi1 = sliderPhi1.Value;
        theta1 = sliderTheta1.Value;
        phi2 = sliderPhi2.Value;
        theta2 = sliderTheta2.Value;

        % Update slider text 
        txtPhi1.String = sprintf('First Spin - Phi: %.2f rad', phi1);
        txtTheta1.String = sprintf('First Spin - Theta: %.2f rad', theta1);
        txtPhi2.String = sprintf('Second Spin - Phi: %.2f rad', phi2);
        txtTheta2.String = sprintf('Second Spin - Theta: %.2f rad', theta2);

        % Cartesian coordinates 
        [x1, y1, z1] = sph2cart(phi1, pi/2 - theta1, 1);
        [x2, y2, z2] = sph2cart(phi2, pi/2 - theta2, 1);

        % Sum of the spins
        xs = x1 + x2;
        ys = y1 + y2;
        zs = z1 + z2;
        norm = sqrt(xs^2 + ys^2 + zs^2);
        xs = xs / norm;
        ys = ys / norm;
        zs = zs / norm;

        % spherical coordinates
        [phiSum, thetaSum, ~] = cart2sph(xs, ys, zs);
        thetaSum = pi/2 - thetaSum; % Adjust for Bloch sphere representation
        txtSumPhiTheta.String = sprintf('Sum Vector - Phi: %.2f rad, Theta: %.2f°', phiSum, rad2deg(thetaSum));

        % Plot both spins and sum
        spinVector1 = plot3([0, x1], [0, y1], [0, z1], 'b', 'LineWidth', 2); hold on;
        spinVector2 = plot3([0, x2], [0, y2], [0, z2], 'r', 'LineWidth', 2);
        sumVector = plot3([0, xs], [0, ys], [0, zs], 'g', 'LineWidth', 2);

        % Bloch sphere with a transparent surface and grid
        [X, Y, Z] = sphere(50);
        sphereSurface = surface(X, Y, Z, 'FaceColor', 'none', 'EdgeColor', [0.5 0.5 0.5], 'EdgeAlpha', 0.5);

        % Axsis
        axis equal; xlabel('x'); ylabel('y'); zlabel('z');
        title('Spin Vectors Summation on Bloch Sphere');
        legend([spinVector1, spinVector2, sumVector], {'First Spin', 'Second Spin', 'Sum of Spins'}, 'Location', 'bestoutside');

        % |0> and |1> state labels
        text(0, 0, 1.1, '|0>', 'HorizontalAlignment', 'center');
        text(0, 0, -1.1, '|1>', 'HorizontalAlignment', 'center');

        view(43, 24); % view angle

        hold off;
    end
end
