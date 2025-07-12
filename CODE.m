%% --- PUSH-DOWN BIOMECCANICA DEL BRACCIO ---
% Questo script analizza e confronta due esecuzioni (errata/corretta) di un movimento biomeccanico di push-down
%% --- PARAMETRI ANTROPOMETRICI ---
% Definire le lunghezze dei segmenti e la forza applicata
l1 = 0.305;          % Lunghezza spalla-gomito in metri
l2 = 0.27;           % Lunghezza gomito-polso in metri
F_ext = 40*9.81;     % Forza verticale applicata al polso [N] (massa*gravità)

%% --- LETTURA DATI DA FILE CSV ---
% Legge i dati dei marker da file per due esecuzioni: errata e corretta
dati_errato = readtable('dati_marker_2.csv');   % esecuzione errata
dati_corretto = readtable('dati_marker_1.csv'); % esecuzione corretta

%% --- ESTRAZIONE DELLE COORDINATE DEI MARKER ---
% Estrae i vettori tempo e le coordinate (x,y) dei tre marker (spalla, gomito, polso) per entrambe le esecuzioni
t_err = dati_errato.time(:);                        % tempo per l'esecuzione errata
P1_err = [dati_errato.x1 dati_errato.y1];           % marker spalla (errato)
P2_err = [dati_errato.x2 dati_errato.y2];           % marker gomito (errato)
P3_err = [dati_errato.x3 dati_errato.y3];           % marker polso (errato)
t_corr = dati_corretto.time(:);                     % tempo per l'esecuzione corretta
P1_corr = [dati_corretto.x1 dati_corretto.y1];      % marker spalla (corretto)
P2_corr = [dati_corretto.x2 dati_corretto.y2];      % marker gomito (corretto)
P3_corr = [dati_corretto.x3 dati_corretto.y3];      % marker polso (corretto)

%% --- CALCOLO ANGOLI E JACOBIANO ---
% Calcola per ogni frame:
% - theta1: angolo tra spalla e gomito
% - theta2: angolo tra gomito e polso
% - J: matrice Jacobiana per la cinematica differenziale del braccio
[theta1_err, theta2_err, J_err] = calcola_cinematica(P1_err, P2_err, P3_err, l1, l2);
[theta1_corr, theta2_corr, J_corr] = calcola_cinematica(P1_corr, P2_corr, P3_corr, l1, l2);

%% --- COSTRUZIONE VETTORE FORZA APPLICATA ---
% Crea una matrice forza costante (verticale verso il basso) applicata al polso per ogni frame
F_err = repmat([0; -F_ext], 1, length(t_err));
F_corr = repmat([0; -F_ext], 1, length(t_corr));

%% --- CALCOLO DELLE COPPIE ARTICOLARI (TAU) ---
% Calcola le coppie articolari richieste in spalla e gomito (tau1 e tau2) tramite il Jacobiano e la forza
tau_err = calcola_tau(J_err, F_err);
tau_corr = calcola_tau(J_corr, F_corr);

%% --- VISUALIZZAZIONE E ANIMAZIONE ---
% Lancia la funzione di animazione che mostra:
% - il movimento dei marker e delle leve (errato/corretto)
% - le traiettorie articolari
% - i grafici delle coppie articolari
% - la tabella dettagliata dei valori
animazione_pushdown_separata(t_err, P1_err, P2_err, P3_err, theta1_err, theta2_err, J_err, tau_err, ...
                             t_corr, P1_corr, P2_corr, P3_corr, theta1_corr, theta2_corr, J_corr, tau_corr, F_ext);

%% --- FUNZIONI DI SUPPORTO ---
function [theta1, theta2, J] = calcola_cinematica(P1, P2, P3, l1, l2)
    N = size(P1,1);
    theta1 = zeros(N,1);
    theta2 = zeros(N,1);
    J = zeros(2,2,N);
    for i = 1:N
        v1 = P2(i,:) - P1(i,:);
        v2 = P3(i,:) - P2(i,:);
        theta1(i) = atan2(v1(2), v1(1));
        theta2(i) = atan2(v2(2), v2(1)) - theta1(i);
        th1 = theta1(i); th2 = theta2(i);
        J(:,:,i) = [-l1*sin(th1)-l2*sin(th1+th2), -l2*sin(th1+th2);
                     l1*cos(th1)+l2*cos(th1+th2),  l2*cos(th1+th2)];
    end
end

function tau = calcola_tau(J, F)
    N = size(J,3);
    tau = zeros(2,N);
    for i = 1:N
        tau(:,i) = J(:,:,i)' * F(:,i);
    end
end

function [xa,ya,ang_c] = arco(C, v_ref, v_mov, r)
    a1 = atan2(v_ref(2), v_ref(1));
    a2 = atan2(v_mov(2), v_mov(1));
    da = wrapToPi(a2-a1);
    if da >= 0
        theta = linspace(a1, a2, 30);
    else
        theta = linspace(a2, a1, 30);
        theta = fliplr(theta);
    end
    xa = C(1) + r*cos(theta);
    ya = C(2) + r*sin(theta);
    ang_c = (theta(1) + theta(end))/2;
end

function [Fx, Fy] = ellisse_di_forza(J)
    theta = linspace(0, 2*pi, 50);
    tau = [cos(theta); sin(theta)];
    if rank(J) < 2 || any(isnan(J(:))) || any(isinf(J(:)))
        Fx = nan(1,50); Fy = nan(1,50);
        return
    end
    F = (J')\tau;
    Fx = F(1,:);
    Fy = F(2,:);
end

function [dX, dY] = diagonale_maggiore_ellisse(Fx, Fy)
    pts = [Fx(:) Fy(:)]';
    [~,~,V] = svd(cov(pts'));
    a = max(sqrt(sum(pts.^2,1)));
    dir = V(:,1);
    dX = [-a*dir(1), a*dir(1)];
    dY = [-a*dir(2), a*dir(2)];
end

function area = area_ellisse_di_forza(Fx, Fy)
    pts = [Fx(:) Fy(:)]';
    if all(~isnan(pts(:)))
        c = mean(pts,2);
        pts0 = pts - c;
        Sigma = cov(pts0');
        [U,S,~] = svd(Sigma);
        a = sqrt(S(1,1))*2; % semiasse maggiore
        b = sqrt(S(2,2))*2; % semiasse minore
        area = pi*a*b;
    else
        area = NaN;
    end
end

function lung = lunghezza_diagonale_ellisse(Fx, Fy)
    % Calcola la lunghezza della diagonale maggiore dell'ellisse di forza (in N)
    if all(~isnan(Fx)) && all(~isnan(Fy))
        [dX, dY] = diagonale_maggiore_ellisse(Fx, Fy);
        lung = sqrt((dX(2)-dX(1))^2 + (dY(2)-dY(1))^2);
    else
        lung = NaN;
    end
end

function dist = distanza_spalla_diagonale(P_spalla, P_polso, Fx, Fy)
    % Calcola la distanza minima tra la spalla e la diagonale maggiore dell'ellisse (braccio del momento)
    if all(~isnan(Fx)) && all(~isnan(Fy))
        [dX, dY] = diagonale_maggiore_ellisse(Fx, Fy);
        % Estremi della diagonale maggiore (nel sistema del polso)
        D1 = P_polso + [dX(1), dY(1)];
        D2 = P_polso + [dX(2), dY(2)];
        % Distanza punto-retta (segmento D1-D2)
        v = D2 - D1;
        w = P_spalla - D1;
        c1 = dot(w,v);
        c2 = dot(v,v);
        b = c1 / c2;
        Pb = D1 + b*v; % Proiezione ortogonale su segmento (può uscire dai limiti)
        % Clamp Pb al segmento
        b = max(0, min(1, b));
        Pb = D1 + b*v;
        dist = norm(P_spalla - Pb);
    else
        dist = NaN;
    end
end

function animazione_pushdown_separata(t1, P1_1, P2_1, P3_1, theta1_1, theta2_1, J1, tau1, ...
                                      t2, P1_2, P2_2, P3_2, theta1_2, theta2_2, J2, tau2, ~)
    fig = figure('Name','Animazione Pushdown + Analisi','NumberTitle','off', ...
                 'Position',[50 50 2200 950]);

    ax1 = subplot(2,2,1); hold(ax1,'on');
    ax2 = subplot(2,2,2); hold(ax2,'on');
    ax3 = subplot(2,2,3); hold(ax3,'on');
    ax4 = subplot(2,2,4); hold(ax4,'on');

    pos1 = get(ax1,'Position'); pos2 = get(ax2,'Position');
    set(ax1, 'Position', [pos1(1)-0.07, pos1(2)-0.09, pos1(3)*1.5, pos1(4)*1.15]);
    set(ax2, 'Position', [pos2(1)-0.12, pos2(2)-0.09, pos2(3)*1.5, pos2(4)*1.15]);
    pos3 = get(ax3,'Position'); pos4 = get(ax4,'Position'); delta = 0.05;
    set(ax3,'Position',[pos3(1)+delta, pos3(2)-0.05, pos3(3)-delta, pos3(4)]);
    set(ax4,'Position',[pos4(1),       pos4(2)-0.05, pos4(3)-delta, pos4(4)]);

    N1 = size(P1_1,1); N2 = size(P1_2,1); N = min(N1,N2);
    t_min = min([t1(1), t2(1)]);
    t_max = max([t1(end), t2(end)]);
    allx = [P1_1(:,1); P2_1(:,1); P3_1(:,1); P1_2(:,1); P2_2(:,1); P3_2(:,1)];
    ally = [P1_1(:,2); P2_1(:,2); P3_1(:,2); P1_2(:,2); P2_2(:,2); P3_2(:,2)];
    buffer = 0.2;
    minx = min(allx)-buffer; maxx = max(allx)+buffer;
    miny = min(ally)-buffer; maxy = max(ally)+buffer;

    h1 = plot(ax1, [P1_1(1,1) P2_1(1,1) P3_1(1,1)], [P1_1(1,2) P2_1(1,2) P3_1(1,2)], 'ro-','LineWidth',5,'MarkerFaceColor','r');
    hF1 = quiver(ax1, P3_1(1,1), P3_1(1,2), 0, -0.15, 'k','LineWidth',2,'MaxHeadSize',4,'AutoScale','off');
    hF1Label = text(ax1, P3_1(1,1)+0.03, P3_1(1,2)-0.09, 'F', 'FontSize',14, 'FontWeight','bold', 'Color','k');
    title(ax1, 'Esecuzione ERRATA');
    xlabel(ax1, 'x [m]'); ylabel(ax1, 'y [m]');
    grid(ax1,'on'); axis(ax1,'equal');
    set(ax1,'XLim',[minx maxx],'YLim',[miny maxy]);
    hArc1 = plot(ax1, nan, nan, 'g-', 'LineWidth',2.5);
    hArc2 = plot(ax1, nan, nan, 'c-', 'LineWidth',2.5);
    hLineRef1 = plot(ax1, nan, nan, 'g--','LineWidth',2.5);
    hLineMov1 = plot(ax1, nan, nan, 'g-','LineWidth',2.5);
    hLineRef2 = plot(ax1, nan, nan, 'c--','LineWidth',2.5);
    hLineMov2 = plot(ax1, nan, nan, 'c-','LineWidth',2.5);
    hSpallaLabel = text(ax1,0,0,'','FontSize',11,'FontWeight','bold','Color','k','HorizontalAlignment','left');
    hGomitoLabel = text(ax1,0,0,'','FontSize',11,'FontWeight','bold','Color','k','HorizontalAlignment','left');
    hPolsoLabel  = text(ax1,0,0,'','FontSize',11,'FontWeight','bold','Color','k','HorizontalAlignment','left');
    hEll1 = plot(ax1, nan, nan, 'm-', 'LineWidth', 2);
    hDiag1 = plot(ax1, nan, nan, 'y-', 'LineWidth', 2);

    h2 = plot(ax2, [P1_2(1,1) P2_2(1,1) P3_2(1,1)], [P1_2(1,2) P2_2(1,2) P3_2(1,2)], 'bo-','LineWidth',5,'MarkerFaceColor','b');
    hF2 = quiver(ax2, P3_2(1,1), P3_2(1,2), 0, -0.15, 'k','LineWidth',2,'MaxHeadSize',4,'AutoScale','off');
    hF2Label = text(ax2, P3_2(1,1)+0.03, P3_2(1,2)-0.09, 'F', 'FontSize',14, 'FontWeight','bold', 'Color','k');
    title(ax2, 'Esecuzione CORRETTA');
    xlabel(ax2, 'x [m]'); ylabel(ax2, 'y [m]');
    grid(ax2,'on'); axis(ax2,'equal');
    set(ax2,'XLim',[minx maxx],'YLim',[miny maxy]);
    hArc3 = plot(ax2, nan, nan, 'g-', 'LineWidth',2.5);
    hArc4 = plot(ax2, nan, nan, 'c-', 'LineWidth',2.5);
    hLineRef3 = plot(ax2, nan, nan, 'g--','LineWidth',2.5);
    hLineMov3 = plot(ax2, nan, nan, 'g-','LineWidth',2.5);
    hLineRef4 = plot(ax2, nan, nan, 'c--','LineWidth',2.5);
    hLineMov4 = plot(ax2, nan, nan, 'c-','LineWidth',2.5);
    hSpallaLabel2 = text(ax2,0,0,'','FontSize',11,'FontWeight','bold','Color','k','HorizontalAlignment','left');
    hGomitoLabel2 = text(ax2,0,0,'','FontSize',11,'FontWeight','bold','Color','k','HorizontalAlignment','left');
    hPolsoLabel2  = text(ax2,0,0,'','FontSize',11,'FontWeight','bold','Color','k','HorizontalAlignment','left');
    hEll2 = plot(ax2, nan, nan, 'm-', 'LineWidth', 2);
    hDiag2 = plot(ax2, nan, nan, 'y-', 'LineWidth', 2);

    varNames = {'x1','y1','x2','y2','x3','y3','θ₁ [deg]','θ₂ [deg]',...
        'J11','J12','J21','J22','τ₁ [Nm]','τ₂ [Nm]'};
    tblDataErr = cell(1,numel(varNames));
    tblDataCorr = cell(1,numel(varNames));
    uitblErr = uitable('Data',tblDataErr,'ColumnName',varNames,...
        'Position',[60 870 870 80],'FontSize',12,'BackgroundColor',[1 0.95 0.95]);
    uitblCorr = uitable('Data',tblDataCorr,'ColumnName',varNames,...
        'Position',[1000 870 870 80],'FontSize',12,'BackgroundColor',[0.95 0.95 1]);

    uitblAlpha = uitable('Data',{'--','--'},'ColumnName',{'errata','corretta'},...
        'Position',[20 540 200 48],'FontSize',12,'BackgroundColor',[1 1 0.9]);
    uicontrol('Style','text','String','Pendenza diagonale maggiore ellisse [gradi]: (risp. asse x)','Position',[20 610 350 40],'FontSize',12,'BackgroundColor',[0.9 0.9 0.9]);

    uitblAreaEll = uitable('Data',{'--','--'},'ColumnName',{'errata','corretta'},...
        'Position',[20 480 200 48],'FontSize',12,'BackgroundColor',[0.9 1 0.9]);
    uicontrol('Style','text','String','Area ellisse di forza [N^2]','Position',[20 530 200 40],'FontSize',12,'BackgroundColor',[0.9 1 0.9]);

    % --- Nuova tabella lunghezza diagonale maggiore ellisse (senza riduzioni) ---
    uitblLungDiag = uitable('Data',{'--','--'},'ColumnName',{'errata','corretta'}, ...
        'Position',[20 420 200 48],'FontSize',12,'BackgroundColor',[0.9 1 1]);
    uicontrol('Style','text','String','Lunghezza diagonale maggiore ellisse [N]','Position',[20 470 200 40],'FontSize',12,'BackgroundColor',[0.9 1 1]);

    % --- Nuova tabella distanza minima tra spalla e diagonale maggiore ellisse, in mm ---
    uitblDistSpallaDiag = uitable('Data',{'--','--'},'ColumnName',{'errata','corretta'}, ...
        'Position',[20 360 200 48],'FontSize',12,'BackgroundColor',[1 0.9 1]);
    uicontrol('Style','text','String','Distanza min. spalla-diagonale maggiore ellisse [mm]','Position',[20 410 200 40],'FontSize',12,'BackgroundColor',[1 0.9 1]);

    uicontrol('Style','text','String','Frame:','Position',[17 180 80 40],'FontSize',12,'BackgroundColor',[0.9 0.9 0.9]);
    editFrame = uicontrol('Style','edit','String','1','Position',[100 190 70 40],'FontSize',12,'BackgroundColor',[0.8 0.8 1]);
    uicontrol('Style','text','String','Tempo [s]:','Position',[20 140 70 40],'FontSize',12,'BackgroundColor',[0.9 0.9 0.9]);
    editTime = uicontrol('Style','edit','String',sprintf('%.3f',t_min),'Position',[100 140 80 40],'FontSize',12);
    btnGoTime = uicontrol('Style','pushbutton','String','Tempo','Position',[180 140 80 40],'FontSize',12);
    btnStart = uicontrol('Style','pushbutton','String','Start','Position',[20 250 70 50],'FontSize',12,'BackgroundColor',[0.6 1 0.6]);
    btnStop  = uicontrol('Style','pushbutton','String','Stop','Position',[90 250 70 50],'FontSize',12,'BackgroundColor',[1 0.6 0.6]);
    btnRestart = uicontrol('Style','pushbutton','String','Restart','Position',[160 250 70 50],'FontSize',12,'BackgroundColor',[1 1 0.6]);
    rangeStr = sprintf('Frame: 1-%d | Tempo: %.3f-%.3f s',N,t_min,t_max);
    uicontrol('Style','text','String',rangeStr,'Position',[10 90 300 25],'FontSize',11,'ForegroundColor','k');

    cla(ax3);
    cla(ax4);
    hTau1Err = plot(ax3, t1(1), tau1(1,1), 'r', 'LineWidth', 1.5); hold(ax3,'on');
    hTau1Corr = plot(ax3, t2(1), tau2(1,1), 'b', 'LineWidth', 1.5);
    xlabel(ax3, 'Tempo [s]');
    ylabel(ax3, '\tau_1 [Nm]');
    title(ax3,'Andamento \tau_1 (spalla) nel tempo');
    legend(ax3, {'Scorretto','Corretto'});
    grid(ax3,'on');
    set(ax3, 'XLim', [t_min t_max]);
    padding = 0.05;
    y1min = min([tau1(1,:), tau2(1,:)]);
    y1max = max([tau1(1,:), tau2(1,:)]);
    yrange1 = y1max - y1min;
    set(ax3, 'YLim', [y1min - padding*yrange1, y1max + padding*yrange1]);

    hTau2Err = plot(ax4, t1(1), tau1(2,1), 'r', 'LineWidth', 1.5); hold(ax4,'on');
    hTau2Corr = plot(ax4, t2(1), tau2(2,1), 'b', 'LineWidth', 1.5);
    xlabel(ax4, 'Tempo [s]');
    ylabel(ax4, '\tau_2 [Nm]');
    title(ax4,'Andamento \tau_2 (gomito) nel tempo');
    legend(ax4, {'Scorretto','Corretto'});
    grid(ax4,'on');
    set(ax4, 'XLim', [t_min t_max]);
    y2min = min([tau1(2,:), tau2(2,:)]);
    y2max = max([tau1(2,:), tau2(2,:)]);
    yrange2 = y2max - y2min;
    set(ax4, 'YLim', [y2min - padding*yrange2, y2max + padding*yrange2]);

    isRunning = false;
    set(btnStart,'Callback',@(src,evt) assignin('base','isRunning',true));
    set(btnStop, 'Callback',@(src,evt) assignin('base','isRunning',false));
    set(btnRestart,'Callback',@(src,evt) restartAnim());
    set(editFrame,'Callback',@(src,evt) goToFrame());
    set(btnGoTime,'Callback',@(src,evt) goToTime());
    set(editTime,'Callback',@(src,evt) goToTime());

    label_off = 0.015; arc_text_off = 0.03;
    frame = 1; lastFrame = -1; lastTime = -1;

    max_ell_size = 1.5; % << LIMITA LA VISUALIZZAZIONE
    riduzione = 0.3;    % << ELLISSE RIDOTTA PER ANIMAZIONE

    function goToFrame()
        frameStr = get(editFrame,'String');
        f = str2double(frameStr);
        if isnan(f), f = 1; end
        frame = max(1, min(N, round(f)));
        set(editFrame,'String',num2str(frame));
        set(editTime,'String',sprintf('%.3f', t1(frame)));
        lastFrame = frame; lastTime = t1(frame);
        updateTables();
        aggiornaCurveTau(frame);
        aggiornaEllisse(frame);
    end

    function goToTime()
        tval = str2double(get(editTime,'String'));
        if isnan(tval), tval = t_min; end
        tval = max(t_min, min(t_max, tval));
        if any(diff(t1)<=0)
            frame = find(t1>=tval,1,'first');
            if isempty(frame), frame = N; end
        else
            [~,frame] = min(abs(t1-tval));
        end
        set(editFrame,'String',num2str(frame));
        set(editTime,'String',sprintf('%.3f', t1(frame)));
        lastFrame = frame; lastTime = t1(frame);
        updateTables();
        aggiornaCurveTau(frame);
        aggiornaEllisse(frame);
    end

    function restartAnim()
        frame = 1;
        set(editFrame,'String','1');
        set(editTime,'String',sprintf('%.3f', t1(1)));
        assignin('base','isRunning',false);
        lastFrame = 1; lastTime = t1(1);
        updateTables();
        aggiornaCurveTau(frame);
        aggiornaEllisse(frame);
    end

    function updateTables()
        tblDataErr = {P1_1(frame,1), P1_1(frame,2), P2_1(frame,1), P2_1(frame,2), ...
                      P3_1(frame,1), P3_1(frame,2), ...
                      rad2deg(theta1_1(frame)), rad2deg(theta2_1(frame)), ...
                      J1(1,1,frame), J1(1,2,frame), J1(2,1,frame), J1(2,2,frame), ...
                      tau1(1,frame), tau1(2,frame)};
        tblDataCorr = {P1_2(frame,1), P1_2(frame,2), P2_2(frame,1), P2_2(frame,2), ...
                       P3_2(frame,1), P3_2(frame,2), ...
                       rad2deg(theta1_2(frame)), rad2deg(theta2_2(frame)), ...
                       J2(1,1,frame), J2(1,2,frame), J2(2,1,frame), J2(2,2,frame), ...
                       tau2(1,frame), tau2(2,frame)};
        set(uitblErr,'Data',tblDataErr);
        set(uitblCorr,'Data',tblDataCorr);
    end

    function aggiornaCurveTau(frame)
        set(hTau1Err, 'XData', t1(1:frame), 'YData', tau1(1,1:frame));
        set(hTau1Corr, 'XData', t2(1:frame), 'YData', tau2(1,1:frame));
        set(hTau2Err, 'XData', t1(1:frame), 'YData', tau1(2,1:frame));
        set(hTau2Corr, 'XData', t2(1:frame), 'YData', tau2(2,1:frame));
        drawnow limitrate
    end

    function aggiornaEllisse(frame)
        % --- ERRATO ---
        Jnow = J1(:,:,frame);
        [Fx1_full, Fy1_full] = ellisse_di_forza(Jnow);
        area1 = area_ellisse_di_forza(Fx1_full, Fy1_full); % calcolo area REALE
        lung_diag1 = lunghezza_diagonale_ellisse(Fx1_full, Fy1_full); % lunghezza diagonale REALE
        dist_spalla_diag1 = distanza_spalla_diagonale(P1_1(frame,:), P3_1(frame,:), Fx1_full, Fy1_full);

        Fx1 = Fx1_full;
        Fy1 = Fy1_full;
        if max(abs(Fx1)) > max_ell_size || max(abs(Fy1)) > max_ell_size
            scale = max_ell_size / max([max(abs(Fx1)), max(abs(Fy1))]);
            Fx1 = Fx1 * scale;
            Fy1 = Fy1 * scale;
        end
        Fx1 = Fx1 * riduzione;
        Fy1 = Fy1 * riduzione;
        if all(~isnan(Fx1))
            [dX1, dY1] = diagonale_maggiore_ellisse(Fx1, Fy1);
            set(hDiag1, 'XData', P3_1(frame,1) + dX1, 'YData', P3_1(frame,2) + dY1);
            alfa1 = atan2(dY1(2)-dY1(1), dX1(2)-dX1(1));
            alfa1deg = rad2deg(alfa1);
            alphaData = get(uitblAlpha, 'Data');
            alphaData{1,1} = sprintf('%.1f',alfa1deg);
            set(uitblAlpha,'Data',alphaData);
        else
            set(hDiag1, 'XData', nan, 'YData', nan);
            alphaData = get(uitblAlpha, 'Data');
            alphaData{1,1} = '--';
            set(uitblAlpha,'Data',alphaData);
        end
        set(hEll1, 'XData', P3_1(frame,1) + Fx1, 'YData', P3_1(frame,2) + Fy1);

        % --- CORRETTO ---
        Jnow2 = J2(:,:,frame);
        [Fx2_full, Fy2_full] = ellisse_di_forza(Jnow2);
        area2 = area_ellisse_di_forza(Fx2_full, Fy2_full); % calcolo area REALE
        lung_diag2 = lunghezza_diagonale_ellisse(Fx2_full, Fy2_full); % lunghezza diagonale REALE
        dist_spalla_diag2 = distanza_spalla_diagonale(P1_2(frame,:), P3_2(frame,:), Fx2_full, Fy2_full);

        Fx2 = Fx2_full;
        Fy2 = Fy2_full;
        if max(abs(Fx2)) > max_ell_size || max(abs(Fy2)) > max_ell_size
            scale2 = max_ell_size / max([max(abs(Fx2)), max(abs(Fy2))]);
            Fx2 = Fx2 * scale2;
            Fy2 = Fy2 * scale2;
        end
        Fx2 = Fx2 * riduzione;
        Fy2 = Fy2 * riduzione;
        if all(~isnan(Fx2))
            [dX2, dY2] = diagonale_maggiore_ellisse(Fx2, Fy2);
            set(hDiag2, 'XData', P3_2(frame,1) + dX2, 'YData', P3_2(frame,2) + dY2);
            alfa2 = atan2(dY2(2)-dY2(1), dX2(2)-dX2(1));
            alfa2deg = rad2deg(alfa2);
            alphaData = get(uitblAlpha, 'Data');
            alphaData{1,2} = sprintf('%.1f',alfa2deg);
            set(uitblAlpha,'Data',alphaData);
        else
            set(hDiag2, 'XData', nan, 'YData', nan);
            alphaData = get(uitblAlpha, 'Data');
            alphaData{1,2} = '--';
            set(uitblAlpha,'Data',alphaData);
        end
        set(hEll2, 'XData', P3_2(frame,1) + Fx2, 'YData', P3_2(frame,2) + Fy2);

        % --- Aggiorna tabella area ellisse ---
        areaData = get(uitblAreaEll, 'Data');
        if ~isnan(area1)
            areaData{1,1} = sprintf('%.2f', area1);
        else
            areaData{1,1} = '--';
        end
        if ~isnan(area2)
            areaData{1,2} = sprintf('%.2f', area2);
        else
            areaData{1,2} = '--';
        end
        set(uitblAreaEll, 'Data', areaData);

        % --- Aggiorna tabella lunghezza diagonale maggiore (senza riduzioni) ---
        lungData = get(uitblLungDiag, 'Data');
        if ~isnan(lung_diag1)
            lungData{1,1} = sprintf('%.2f', lung_diag1);
        else
            lungData{1,1} = '--';
        end
        if ~isnan(lung_diag2)
            lungData{1,2} = sprintf('%.2f', lung_diag2);
        else
            lungData{1,2} = '--';
        end
        set(uitblLungDiag, 'Data', lungData);

        % --- Aggiorna tabella distanza minima spalla-diagonale maggiore ellisse (in mm) ---
        dist_mm1 = dist_spalla_diag1 * 1000;
        dist_mm2 = dist_spalla_diag2 * 1000;
        distData = get(uitblDistSpallaDiag, 'Data');
        if ~isnan(dist_mm1)
            distData{1,1} = sprintf('%.1f', dist_mm1);
        else
            distData{1,1} = '--';
        end
        if ~isnan(dist_mm2)
            distData{1,2} = sprintf('%.1f', dist_mm2);
        else
            distData{1,2} = '--';
        end
        set(uitblDistSpallaDiag, 'Data', distData);
    end

    assignin('base','isRunning',false);
    updateTables();
    aggiornaCurveTau(frame);
    aggiornaEllisse(frame);
    while ishandle(fig)
        isRunning = evalin('base','isRunning');
        frameStr = get(editFrame,'String');
        tStr = get(editTime,'String');
        f = str2double(frameStr);
        tval = str2double(tStr);
        if isnan(f), f = 1; end
        if isnan(tval), tval = t_min; end

        if abs(tval-lastTime) > 1e-8
            tval = max(t_min, min(t_max, tval));
            if any(diff(t1)<=0)
                frame = find(t1>=tval,1,'first');
                if isempty(frame), frame = N; end
            else
                [~,frame] = min(abs(t1-tval));
            end
            set(editFrame,'String',num2str(frame));
            set(editTime,'String',sprintf('%.3f', t1(frame)));
            lastFrame = frame; lastTime = t1(frame);
            updateTables();
            aggiornaCurveTau(frame);
            aggiornaEllisse(frame);
        elseif f ~= lastFrame
            frame = max(1, min(N, round(f)));
            set(editFrame,'String',num2str(frame));
            set(editTime,'String',sprintf('%.3f', t1(frame)));
            lastFrame = frame; lastTime = t1(frame);
            updateTables();
            aggiornaCurveTau(frame);
            aggiornaEllisse(frame);
        end

        set(h1, 'XData',[P1_1(frame,1) P2_1(frame,1) P3_1(frame,1)], ...
                'YData',[P1_1(frame,2) P2_1(frame,2) P3_1(frame,2)]);
        set(hF1, 'XData', P3_1(frame,1), 'YData', P3_1(frame,2), ...
            'UData', 0, 'VData', -0.15);
        set(hF1Label, 'Position', [P3_1(frame,1)+0.03, P3_1(frame,2)-0.09]);
        vRef1 = [1 0]; vMov1 = P2_1(frame,:) - P1_1(frame,:);
        [xa,ya,~] = arco(P1_1(frame,:), vRef1, vMov1, 0.07);
        set(hArc1,'XData',xa,'YData',ya);
        set(hLineRef1,'XData',[P1_1(frame,1), P1_1(frame,1)+0.07*vRef1(1)], ...
            'YData',[P1_1(frame,2), P1_1(frame,2)+0.07*vRef1(2)]);
        vMov1u = vMov1/norm(vMov1);
        set(hLineMov1,'XData',[P1_1(frame,1), P1_1(frame,1)+0.07*vMov1u(1)], ...
            'YData',[P1_1(frame,2), P1_1(frame,2)+0.07*vMov1u(2)]);
        vRef2 = P2_1(frame,:) - P1_1(frame,:); vMov2 = P3_1(frame,:) - P2_1(frame,:);
        [xa2,ya2,~] = arco(P2_1(frame,:), vRef2, vMov2, 0.07);
        set(hArc2,'XData',xa2,'YData',ya2);
        vRef2u = vRef2/norm(vRef2); vMov2u = vMov2/norm(vMov2);
        set(hLineRef2,'XData',[P2_1(frame,1), P2_1(frame,1)+0.07*vRef2u(1)], ...
            'YData',[P2_1(frame,2), P2_1(frame,2)+0.07*vRef2u(2)]);
        set(hLineMov2,'XData',[P2_1(frame,1), P2_1(frame,1)+0.07*vMov2u(1)], ...
            'YData',[P2_1(frame,2), P2_1(frame,2)+0.07*vMov2u(2)]);
        set(hSpallaLabel, 'Position', P1_1(frame,:) + [+0.03 +0.07], 'String', 'spalla');
        set(hGomitoLabel, 'Position', P2_1(frame,:) + [-0.18 -0.015], 'String', 'gomito');
        set(hPolsoLabel,  'Position', P3_1(frame,:) + [+0.05 +0.02], 'String', 'polso');

        set(h2, 'XData',[P1_2(frame,1) P2_2(frame,1) P3_2(frame,1)], ...
                'YData',[P1_2(frame,2) P2_2(frame,2) P3_2(frame,2)]);
        set(hF2, 'XData', P3_2(frame,1), 'YData', P3_2(frame,2), ...
            'UData', 0, 'VData', -0.15);
        set(hF2Label, 'Position', [P3_2(frame,1)+0.03, P3_2(frame,2)-0.09]);
        vRef3 = [1 0]; vMov3 = P2_2(frame,:) - P1_2(frame,:);
        [xa3,ya3,~] = arco(P1_2(frame,:), vRef3, vMov3, 0.07);
        set(hArc3,'XData',xa3,'YData',ya3);
        set(hLineRef3,'XData',[P1_2(frame,1), P1_2(frame,1)+0.07*vRef3(1)], ...
            'YData',[P1_2(frame,2), P1_2(frame,2)+0.07*vRef3(2)]);
        vMov3u = vMov3/norm(vMov3);
        set(hLineMov3,'XData',[P1_2(frame,1), P1_2(frame,1)+0.07*vMov3u(1)], ...
            'YData',[P1_2(frame,2), P1_2(frame,2)+0.07*vMov3u(2)]);
        vRef4 = P2_2(frame,:) - P1_2(frame,:); vMov4 = P3_2(frame,:) - P2_2(frame,:);
        [xa4,ya4,~] = arco(P2_2(frame,:), vRef4, vMov4, 0.07);
        set(hArc4,'XData',xa4,'YData',ya4);
        vRef4u = vRef4/norm(vRef4); vMov4u = vMov4/norm(vMov4);
        set(hLineRef4,'XData',[P2_2(frame,1), P2_2(frame,1)+0.07*vRef4u(1)], ...
            'YData',[P2_2(frame,2), P2_2(frame,2)+0.07*vRef4u(2)]);
        set(hLineMov4,'XData',[P2_2(frame,1), P2_2(frame,1)+0.07*vMov4u(1)], ...
            'YData',[P2_2(frame,2), P2_2(frame,2)+0.07*vMov4u(2)]);
        set(hSpallaLabel2, 'Position', P1_2(frame,:) + [+0.03 +0.07], 'String', 'spalla');
        set(hGomitoLabel2, 'Position', P2_2(frame,:) + [-0.18 -0.015], 'String', 'gomito');
        set(hPolsoLabel2,  'Position', P3_2(frame,:) + [+0.05 +0.02], 'String', 'polso');

        aggiornaCurveTau(frame);
        aggiornaEllisse(frame);

        drawnow;
        if isRunning
            pause(0.03);
            frame = frame + 1;
            if frame > N, frame = 1; end
            set(editFrame,'String',num2str(frame));
            set(editTime,'String',sprintf('%.3f', t1(frame)));
            lastFrame = frame; lastTime = t1(frame);
            updateTables();
            aggiornaCurveTau(frame);
            aggiornaEllisse(frame);
        else
            pause(0.05);
        end
    end
end
