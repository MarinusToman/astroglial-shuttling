%% Environment clean-up
clc         % Clear command window
clear       % Clear workspace
close all   % Close all figures
%% Simulation parameters
simTime     = 100;              % Duration of the simulation; unit: s;
iStart      = 10;               % Time to start presynaptic stimulation; unit: s;
iStop       = 60;               % Time to stop presynaptic stimulation; unit: s;
I_mag       = 10;               % Magnitude of external stimulus; unit: uA/cm^2
freq        = 10;               % Presynaptic stimulation frequency; unit: Hz;
dt          = 1e-5;             % Time step size; unit: s; Typically 10us(1e-5s);
t           = 0:dt:simTime;     % Time points
steps       = simTime / dt;     % Number of time steps in simulation
savePeriod  = 100;              % Save only every 100 time steps
resCount    = 1;                % Save result counter
saveCount   = savePeriod;       % Save counter
morph       = Model_Morphology; % Get cellular morphology
% Synapse coverage
synCov = 'c'; % Synapse coverage configuration (A,B,C)
switch synCov
    case {'A','a'}
        morph.ecsVol = morph.cleftVol;
    case {'C','c'}
        morph.ecsVol = morph.ecsVol * 7;
    otherwise
        morph.ecsVol = morph.ecsVol;
end
% Plot and save info
freqChar    = num2str(freq);    % Get freq as a string
curDate     = Get_Date_String;  % Get current date in string form
figureFile  = ['Astro_Shuttling_Model_SynCov',synCov,'_',freqChar,'Hz_',curDate]; % Figure file location
saveFile    = [figureFile,'.mat']; % Save file location
% Other
pulseWidth  = 3;                % External stimulus pulse width; unit: ms
Start_App   = iStart/dt;
Stop_App    = iStop/dt;
spikeFlag   = 0;
poissonFlag = 0; % Poisson stimulation?
bAP         = 1; % Back proagating AP
preb4post   = 1; % Pre fires before post
spDiff      = 2; % Difference between pre-post spikes; unit: ms
V_ast_var   = 0; % Variable astrocyte membrane potential?
ISI         = round((1/freq)/dt); % Interspike interval
pulseWidth  = pulseWidth * 100; % Convert to same unit as dt; unit: 10us
% Constants
zNa         = 1;
zK          = 1;
zCa         = 2;
F           = 96485.33;         % Faradays constant; unit: C/mol;
Cm_post     = 0.01; % Postsynaptic membrane capacitance; unit: F/m^2;
Cm_ast      = 0.01; % Astrocytic membrane capacitance; unit: F/m^2;
gNMDA=1;
gAMPA=1;
%% Variable generation
store       = Model_Storage(steps,savePeriod); % Variable generation
%% Initial conditions
% Steady state
V0_ast      = -0.080742766;     % Astrocytic membrane potential; unit: V
V0_pre      = -0.07;            % Presynaptic resting potential; unit: V
V0_post     = -0.07;            % Postsynaptic resting potential; unit: V
% Concentrations - Ref: Wade et al. (2019)
Ca_neu      = 50e-9;            % Ca concentration in neuron; unit: M
Ca_e        = 1.5e-3;           % Ca concentration in ECS; unit: M
Ca_i        = 100e-9;           % Ca concentration in astrocyte; unit: M
Na_e        = 135e-3;           % Na concentration in ECS; unit: M
Na_i        = 15e-3;            % Na concentration in Ast; unit: M
K_e         = 4e-3;             % K concentration in GECS; unit: M
K_i         = 100e-3;           % K concentration in Ast; unit: M
% Initial conditions
% Voltage
V_pre       = V0_pre;           % Presynaptic membrane potential; unit: V
V_post      = V0_post;          % Postsynaptic membrane potential; unit: V
V_ast       = V0_ast;           % Astrocytic membrane potential; unit: V
% Concentrations
Na_psc      = Na_i;             % Na concentration in PsC; unit: M (mol/L)
Na_ecs      = Na_e;             % Na concentration in ECS; unit: M
Na_pre      = Na_i;             % Na concentration in Pre; unit: M
Na_post     = Na_i;             % Na concentration in Post; unit: M
K_psc       = K_i;              % K concentration in PsC; unit: M
K_ecs       = K_e;              % K concentration in ECS; unit: M
K_pre       = K_i;              % K concentration in Pre; unit: M
K_post      = K_i;              % K concentration in Post; unit: M
Ca_pre      = Ca_neu;           % Ca concentration in PreSyn; unit: M
Ca_post     = Ca_neu;           % Ca concentration in PostSyn; unit: M
Ca_ecs      = Ca_e;             % Ca concentration in ECS; unit: M
Ca_psc      = Ca_i;             % Ca concentration in PsC; unit: M
Ca_ast      = Ca_i;             % Ca concentration in Ast; unit: M
Glu_ecs     = 0;                % Glutamate concentration
JNaEAAT     = 0;
% Neuron params
GluRel      = 0;                % Glutamate release
I_pre       = 0;
nmdar       = 0;
ampar       = 0;
I_stim      = 0;                % External stimulus
I_stim_post = 0;
Q = 0.01;
R = 0;
spDiff = spDiff * 100; % Convert to 10us (1 dt)
[aM, bM]    = m_equations(V_pre, V0_pre);
[aN, bN]    = n_equations(V_pre, V0_pre);
[aH, bH]    = h_equations(V_pre, V0_pre);
[m_inf, ~]  = mca_equations(V_pre*1000);
xHH(1)      = (aM / (aM + bM));
xHH(2)      = (aN / (aN + bN));
xHH(3)      = (aH / (aH + bH));
xHH(4)      = m_inf;
E_Ca        = Nernst_Potential(Ca_ecs, Ca_pre, zCa) * 1000;

% if poissonFlag
    rng(1); % seed random number generator
    vt1 = rand(size(t));
    spikes1 = (freq*dt) > vt1;%poisson
    vt2 = rand(size(t));
    spikes2 = (freq*dt) > vt2;%poisson
% end
%% Presimulation setup
tic; % Start stopwatch timer

for i=1:10000
    %% Presynaptic neuron
    % Active currents
    [V_pre,xHH,INa_vg,IK_vg,ICa_vg,~] = ... 
        Neuron_Model(V_pre, V0_pre, xHH, dt, I_stim,I_pre,E_Ca); % AP
    [INa_nka_pre,IK_nka_pre] = NKA(K_ecs,Na_pre,0); % NKA
    [ICa_pm_pre] = PMCA(Ca_pre); % PMCA 
    % Total currents
    INa_pre     = INa_nka_pre;
    IK_pre      = IK_nka_pre;
    ICa_pre     = ICa_vg + ICa_pm_pre;
    % Get leak conductances
    [gNa_l_pre] = Leak_Conduct(INa_pre,Na_ecs,Na_pre,V_pre,zNa);
    [gK_l_pre]  = Leak_Conduct(IK_pre,K_ecs,K_pre,V_pre,zK);
    [gCa_l_pre] = Leak_Conduct(ICa_pre,Ca_ecs,Ca_pre,V_pre,zCa);
    % Get passive leak currents
    [INa_l_pre] = Leak(gNa_l_pre,Na_ecs,Na_pre,V_pre,zNa);
    [IK_l_pre]  = Leak(gK_l_pre,K_ecs,K_pre,V_pre,zK);
    [ICa_l_pre] = Leak(gCa_l_pre,Ca_ecs,Ca_pre,V_pre,zCa);
    % Total currents
    INa_pre     = INa_nka_pre + INa_l_pre;
    IK_pre      = IK_nka_pre + IK_l_pre;
    ICa_pre     = ICa_vg + ICa_pm_pre + ICa_l_pre;
    I_pre       = INa_pre + IK_pre + ICa_pre;
    %% Postsynapse
    % Active currents
    [INa_nka_post,IK_nka_post] = NKA(K_ecs,Na_post,0); % NKA
    [INa_nmda_post,ICa_nmda_post,IK_nmda_post,nmda_dr] = NMDA(Glu_ecs, nmdar, V_post,gNMDA); % NMDA
    [INa_ampa_post,ampa_dr] = AMPA(Glu_ecs, ampar, V_post,gAMPA); % AMPA
    [ICa_pm_post] = PMCA(Ca_post);
    % Total currents
    INa_post    = INa_nmda_post + INa_nka_post + INa_ampa_post;
    IK_post     = IK_nka_post + IK_nmda_post;
    ICa_post    = ICa_nmda_post + ICa_pm_post;
    % Get leak conductances
    [gNa_l_post]= Leak_Conduct(INa_post,Na_ecs,Na_post,V_post,zNa);
    [gK_l_post] = Leak_Conduct(IK_post,K_ecs,K_post,V_post,zK);
    [gCa_l_post]= Leak_Conduct(ICa_post,Ca_ecs,Ca_post,V_post,zCa);
    % Get passive leak currents
    [INa_l_post]= Leak(gNa_l_post,Na_ecs,Na_post,V_post,zNa);
    [IK_l_post] = Leak(gK_l_post,K_ecs,K_post,V_post,zK);
    [ICa_l_post]= Leak(gCa_l_post,Ca_ecs,Ca_post,V_post,zCa);
    % Total currents
    INa_post    = INa_nmda_post + INa_ampa_post+ INa_nka_post + INa_l_post;
    IK_post     = IK_nka_post + IK_l_post + IK_nmda_post;
    ICa_post    = ICa_nmda_post + ICa_l_post + ICa_pm_post;
    I_post      = INa_post+IK_post+ICa_post;
    dv_post     = -I_post * (dt/Cm_post);%
    V_post      = V_post + dv_post; % New membrane potential
    % Update NMDA activation
    nmdar = nmdar + (dt * nmda_dr);
    % Update AMPA activation
    ampar = ampar + (dt * ampa_dr);
    %% Astrocyte
    % Active currents
    [INa_ncx, ICa_ncx] = NCX(Na_ecs,Na_psc,Ca_ecs,Ca_psc,V_ast);% NCX
    [INa_nka_psc, IK_nka_psc] = NKA(K_ecs,Na_psc,1); % NKA
    [IK_kir] = Kir(K_ecs,K_psc,V_ast); % Kir
    % Process
%     [INa_pf] = Process(Na_i,Na_psc,zNa,V_ast,V0_ast,morph.processLength);
%     [IK_pf] = Process(K_i,K_psc,zK,V_ast,V0_ast,morph.processLength);
%     [ICa_pf] = Process(Ca_ast,Ca_psc,zCa,V_ast,V0_ast,morph.processLength);
    % Total currents
    INa_psc_tot = INa_ncx + INa_nka_psc;
    IK_psc_tot = IK_nka_psc + IK_kir;
    ICa_psc_tot = ICa_ncx;
    % Get leak conductances
    [gNa_l_psc] = Leak_Conduct(INa_psc_tot,Na_ecs,Na_psc,V_ast,zNa);
    [gK_l_psc] = Leak_Conduct(IK_psc_tot,K_ecs,K_psc,V_ast,zK);
    [gCa_l_psc] = Leak_Conduct(ICa_psc_tot,Ca_ecs,Ca_psc,V_ast,zCa);
    % Get passive leak currents
    [INa_l_psc] = Leak(gNa_l_psc,Na_ecs,Na_psc,V_ast,zNa);
    [IK_l_psc] = Leak(gK_l_psc,K_ecs,K_psc,V_ast,zK);
    [ICa_l_psc] = Leak(gCa_l_psc,Ca_ecs,Ca_psc,V_ast,zCa);
    % Total Astrocyte currents
    INa_psc_tot = INa_ncx + INa_nka_psc + INa_l_psc;
    IK_psc_tot = IK_nka_psc + IK_kir + IK_l_psc;
    ICa_psc_tot = ICa_ncx + ICa_l_psc;

    [Glu_ecs,Q,R] = Glu_Rel_Ca(Ca_pre,Glu_ecs,Q,R,dt);
end
% Display elapsed time
clc;
elapsed = toc; % get elapsed time
fprintf(['Presimulation took ',num2str(elapsed), ' seconds']);

%% Start simulation
tic; % Start stopwatch timer

for i = 1 : steps
    %% External stimulus for presynaptic neuron
    if (i >= Start_App && i <= Stop_App) % Single stimulus
        if poissonFlag
            I_stim = I_mag*100*spikes1(i);
%             I_stim_post = I_mag*spikes2(i-spDiff)*5;
            if preb4post
                I_stim_post = I_mag*spikes2(i-spDiff)*5;
            else
                I_stim_post = I_mag*spikes2(i+spDiff)*5;
            end
        else
            I_stim = I_mag*(mod(i-1,ISI)<=pulseWidth); % External stimulus magnitude
%             I_stim_post = I_mag*spikes2(i)*5;
            if preb4post
                I_stim_post = I_mag*(mod(i-spDiff,ISI)<=pulseWidth)/30;
            else
                I_stim_post = I_mag*(mod(i+spDiff,ISI)<=pulseWidth)/30;
            end
        end
    else % Otherwise no stimulation
        I_stim = 0; % External stimulus magnitude
    end
    %% Presynaptic neuron
    % Active currents
    E_Ca = Nernst_Potential(Ca_ecs, Ca_pre, zCa) * 1000;
    [V_pre,xHH,INa_vg,IK_vg,ICa_vg,~] = ... 
        Neuron_Model(V_pre, V0_pre, xHH, dt, I_stim,I_pre,E_Ca);
    [INa_nka_pre,IK_nka_pre] = NKA(K_ecs,Na_pre,0); % NKA
    [ICa_pm_pre] = PMCA(Ca_pre); % PMCA
    % Leak currents
    INa_l_pre = Leak(gNa_l_pre,Na_ecs,Na_pre,V_pre,zNa);
    IK_l_pre = Leak(gK_l_pre,K_ecs,K_pre,V_pre,zK);
    ICa_l_pre = Leak(gCa_l_pre,Ca_ecs,Ca_pre,V_pre,zCa); 
    % Glutamate release
    [Glu_ecs,Q,R] = Glu_Rel_Ca(Ca_pre,Glu_ecs,Q,R,dt);
    %% Postsynaptic neuron
    % Active currents
    [INa_nka_post,IK_nka_post] = NKA(K_ecs,Na_post,0); % NKA
    [INa_nmda_post,ICa_nmda_post,IK_nmda_post,nmda_dr] = NMDA(Glu_ecs, nmdar, V_post,gNMDA); % NMDA
    [INa_ampa_post,ampa_dr] = AMPA(Glu_ecs, ampar, V_post,gAMPA); % AMPA
    [ICa_pm_post] = PMCA(Ca_post);
    % Leak currents
    [INa_l_post] = Leak(gNa_l_post,Na_ecs,Na_post,V_post,zNa);
    [IK_l_post] = Leak(gK_l_post,K_ecs,K_post,V_post,zK);
    [ICa_l_post] = Leak(gCa_l_post,Ca_ecs,Ca_post,V_post,zCa);
    %% Astrocyte currents/fluxes
    % Active currents
    [INaEAAT,IKEAAT] = EAAT2(Glu_ecs);
    [INa_ncx, ICa_ncx] = NCX(Na_ecs,Na_psc,Ca_ecs,Ca_psc,V_ast); % NCX
    [INa_nka_psc, IK_nka_psc] = NKA(K_ecs,Na_psc,1); % NKA
    IK_kir = Kir(K_ecs,K_psc,V_ast); % Kir
    % Leak currents
    INa_l_psc = Leak(gNa_l_psc,Na_ecs,Na_psc,V_ast,zNa);
    IK_l_psc = Leak(gK_l_psc,K_ecs,K_psc,V_ast,zK);
    ICa_l_psc = Leak(gCa_l_psc,Ca_ecs,Ca_psc,V_ast,zCa);
    % Process
    [INa_pf] = Process(Na_i,Na_psc,zNa,V_ast,V0_ast,morph.processLength);
    [IK_pf] = Process(K_i,K_psc,zK,V_ast,V0_ast,morph.processLength);
    [ICa_pf] = Process(Ca_ast,Ca_psc,zCa,V_ast,V0_ast,morph.processLength);
    %% ECS currents
    INa_dif_ecs = ECS_Diff(Na_ecs,Na_e,zNa);
    IK_dif_ecs = ECS_Diff(K_ecs,K_e,zK);
    ICa_dif_ecs = ECS_Diff(Ca_ecs,Ca_e,zCa);
    %% Update voltage
    % postsynaptic
    INa_post        = INa_nmda_post + INa_nka_post + INa_ampa_post+INa_l_post;
    IK_post         = IK_nka_post + IK_nmda_post +IK_l_post;
    ICa_post        = ICa_nmda_post + ICa_pm_post + ICa_l_post;
    if bAP
        I_post          = INa_post+IK_post+ICa_post-I_stim_post;
    else
        I_post          = INa_post+IK_post+ICa_post;
    end
    dv_post         = -I_post * (dt/Cm_post);%
    V_post          = V_post + dv_post; % New membrane potential
    % astrocytic
    INa_psc_tot     = INa_ncx + INa_nka_psc + INa_l_psc + INaEAAT;
    IK_psc_tot      = IK_nka_psc + IK_kir + IK_l_psc + IKEAAT;
    ICa_psc_tot     = ICa_ncx + ICa_l_psc;
    I_ast           = INa_psc_tot + IK_psc_tot + ICa_psc_tot;
    if V_ast_var
        dv_ast          = -I_ast * (dt/Cm_ast); % Membrane potential change
        V_ast           = V_ast + dv_ast; % New membrane potential
    end
    %% Convert current densities to currents
    % Presynapse
    INa_vg          = INa_vg * morph.synSA;
    IK_vg           = IK_vg * morph.synSA;
    ICa_vg          = ICa_vg * morph.synSA;
    INa_nka_pre     = INa_nka_pre * morph.synSA;
    IK_nka_pre      = IK_nka_pre * morph.synSA;
    ICa_pm_pre      = ICa_pm_pre * morph.synSA;
    INa_l_pre       = INa_l_pre * morph.synSA;
    IK_l_pre        = IK_l_pre * morph.synSA;
    ICa_l_pre       = ICa_l_pre * morph.synSA;
    % Postsynapse
    INa_nmda_post   = INa_nmda_post * morph.synSA;
    ICa_nmda_post   = ICa_nmda_post * morph.synSA;
    IK_nmda_post    = IK_nmda_post * morph.synSA;
    INa_nka_post    = INa_nka_post * morph.synSA;
    IK_nka_post     = IK_nka_post * morph.synSA;
    INa_l_post      = INa_l_post * morph.synSA;
    IK_l_post       = IK_l_post * morph.synSA;
    ICa_l_post      = ICa_l_post * morph.synSA;
    INa_ampa_post   = INa_ampa_post * morph.synSA;
    ICa_pm_post     = ICa_pm_post * morph.synSA;
    % Astrocyte
    INa_ncx         = INa_ncx * morph.pscSA; 
    INa_nka_psc     = INa_nka_psc * morph.pscSA;     
    INaEAAT         = INaEAAT * morph.pscSA;
    IK_nka_psc      = IK_nka_psc * morph.pscSA;
    IK_kir          = IK_kir * morph.pscSA;     
    IKEAAT          = IKEAAT * morph.pscSA;
    ICa_ncx         = ICa_ncx * morph.pscSA; 
    ICa_l_psc       = ICa_l_psc * morph.pscSA;
    IK_l_psc        = IK_l_psc * morph.pscSA;
    INa_l_psc       = INa_l_psc * morph.pscSA;
    % Process
    INa_pf          = INa_pf * morph.processCSA;
    IK_pf           = IK_pf * morph.processCSA;
    ICa_pf          = ICa_pf * morph.processCSA;
    % ECS
    INa_dif_ecs     = INa_dif_ecs * morph.ecsDiffSA;
    IK_dif_ecs      = IK_dif_ecs * morph.ecsDiffSA;
    ICa_dif_ecs     = ICa_dif_ecs * morph.ecsDiffSA;
    switch synCov
        case {'A','a'}
            INa_dif_ecs = 0;
            ICa_dif_ecs = 0;
            IK_dif_ecs = 0;
    end
    %% Total membrane currents
    % Presynapse
    INa_pre         = INa_l_pre + INa_nka_pre;
    IK_pre          = IK_l_pre + IK_nka_pre;
    ICa_pre         = ICa_vg + ICa_l_pre + ICa_pm_pre;
    I_pre           = INa_pre + IK_pre + ICa_pre;
    % Postsynapse
    INa_post        = INa_nmda_post + INa_nka_post + INa_l_post + INa_ampa_post;
    IK_post         = IK_nka_post + IK_l_post + IK_nmda_post;
    ICa_post        = ICa_nmda_post + ICa_l_post + ICa_pm_post;
    I_post          = INa_post + IK_post + ICa_post;
    % Astrocyte
    INa_psc_tot     = INa_ncx + INa_nka_psc + INa_l_psc + INa_pf + INaEAAT;
    IK_psc_tot      = IK_nka_psc + IK_kir + IK_l_psc + IK_pf + IKEAAT;
    ICa_psc_tot     = ICa_ncx + ICa_l_psc + ICa_pf;
    INa_psc         = INa_ncx + INa_nka_psc + INa_l_psc + INaEAAT;
    IK_psc          = IK_nka_psc + IK_kir + IK_l_psc + IKEAAT;
    ICa_psc         = ICa_ncx + ICa_l_psc;
    % ECS 
    INa_ecs         = INa_dif_ecs - INa_psc - INa_post - INa_pre;
    IK_ecs          = IK_dif_ecs - IK_psc - IK_post - IK_pre;
    ICa_ecs         = ICa_dif_ecs - ICa_psc - ICa_pre - ICa_post;
%     ICa_ecs         = ICa_dif_ecs - ICa_psc - ICa_pre;%increase NMDA 10/02/22
    %% Convert currents to fluxes
    % Presynapse
%     JNa_pre         = -INa_pre * (1 / (zNa*F*morph.volSyn));
%     JK_pre          = -IK_pre * (1 / (zK*F*morph.volSyn));
    JCa_pre         = -ICa_pre * (1 / (zCa*F*morph.volSyn));
    % Postsynapse
%     JNa_post        = -INa_post * (1 / (zNa*F*morph.volSyn));
%     JK_post         = -IK_post * (1 / (zK*F*morph.volSyn));
    JCa_post        = -ICa_post * (1 / (zCa*F*morph.volSyn));
    % Astrocyte
    JNa_psc         = -INa_psc_tot * (1 / (zNa*F*morph.pscVol));
    JK_psc          = -IK_psc_tot * (1 / (zK*F*morph.pscVol));
    JCa_psc         = -ICa_psc_tot * (1 / (zCa*F*morph.pscVol));
    % ECS 
    JNa_ecs         = -INa_ecs * (1 / (zNa*F*morph.ecsVol));
    JK_ecs          = -IK_ecs * (1 / (zK*F*morph.ecsVol));
    JCa_ecs         = -ICa_ecs * (1 / (zCa*F*morph.ecsVol));
    %% Solve for change in concentrations (using Euler's forward method)
    % Presynapse
%     Na_pre          = Na_pre + (dt * JNa_pre);
%     K_pre           = K_pre + (dt * JK_pre);
    Ca_pre          = Ca_pre + (dt * JCa_pre);
    % Postsynapse
%     Na_post         = Na_post + (dt * JNa_post);
%     K_post          = K_post + (dt * JK_post);
%     Ca_post         = Ca_post + (dt * JCa_post);%increase NMDA 10/02/22
    nmdar           = nmdar + (dt * nmda_dr); % NMDA activation
    ampar           = ampar + (dt * ampa_dr); % AMPA activation
    % Astrocyte
    Na_psc          = Na_psc + (dt * JNa_psc);
    K_psc           = K_psc + (dt * JK_psc);
    Ca_psc          = Ca_psc + (dt * JCa_psc);
    % ECS
    Na_ecs          = Na_ecs + (dt * JNa_ecs);
    K_ecs           = K_ecs + (dt * JK_ecs);
    Ca_ecs          = Ca_ecs + (dt * JCa_ecs);  
    %% Check concentration hasn't dropped too low
    if Na_pre <= 25e-9
        Na_pre = 25e-9;
    end
    if K_pre <= 25e-9
        K_pre = 25e-9;
    end
    if Ca_pre <= 25e-9
        Ca_pre = 25e-9;
    end
    if Na_post <= 25e-9
        Na_post = 25e-9;
    end
    if K_post <= 25e-9
        K_post = 25e-9;
    end
    if Ca_post <= 25e-9
        Ca_post = 25e-9;
    end
    if Na_psc <= 25e-9
        Na_psc = 25e-9;
    end
    if K_psc <= 25e-9
        K_psc = 25e-9;
    end
    if Ca_psc <= 25e-9
        Ca_psc = 25e-9;
    end
    if Na_ecs <= 25e-9
        Na_ecs = 25e-9;
    end
    if K_ecs <= 25e-9
        K_ecs = 25e-9;
    end
    if Ca_ecs <= 25e-9
        Ca_ecs = 25e-9;
    end
    %% Save results
    if saveCount == savePeriod
        % Voltage
        store.V_pre(resCount,1)         = V_pre;
        store.V_post(resCount,1)        = V_post;
        store.V_ast(resCount,1)         = V_ast;
        % Concentrations
        %   Presynapse
        store.Na_pre(resCount,1)        = Na_pre;
        store.K_pre(resCount,1)         = K_pre;
        store.Ca_pre(resCount,1)        = Ca_pre;
        %   Postsynapse
        store.Na_post(resCount,1)       = Na_post;
        store.K_post(resCount,1)        = K_post;
        store.Ca_post(resCount,1)       = Ca_post;
        %   Astrocyte
        store.Na_psc(resCount,1)        = Na_psc;
        store.K_psc(resCount,1)         = K_psc;
        store.Ca_psc(resCount,1)        = Ca_psc;
        %   ECS
        store.Na_ecs(resCount,1)        = Na_ecs;
        store.K_ecs(resCount,1)         = K_ecs;
        store.Ca_ecs(resCount,1)        = Ca_ecs;
        store.Glu_ecs(resCount,1)       = Glu_ecs;
        % Currents
        %   Presynapse
        store.INa_vg_pre(resCount,1)    = INa_vg;
        store.IK_vg_pre(resCount,1)     = IK_vg;
        store.ICa_vg_pre(resCount,1)    = ICa_vg;
        store.INa_NKA_pre(resCount,1)   = INa_nka_pre;
        store.IK_NKA_pre(resCount,1)    = IK_nka_pre;
        store.ICa_pm_pre(resCount,1)    = ICa_pm_pre;
        store.INa_l_pre(resCount,1)     = INa_l_pre;
        store.IK_l_pre(resCount,1)      = IK_l_pre;
        store.ICa_l_pre(resCount,1)     = ICa_l_pre;
        %   Postsynapse
        store.INa_nka_post(resCount,1)  = INa_nka_post;
        store.IK_nka_post(resCount,1)   = IK_nka_post;
        store.INa_nmda_post(resCount,1) = INa_nmda_post;
        store.ICa_nmda_post(resCount,1) = ICa_nmda_post;
        store.IK_nmda_post(resCount,1)  = IK_nmda_post;
        store.INa_l_post(resCount,1)    = INa_l_post;
        store.IK_l_post(resCount,1)     = IK_l_post;
        store.ICa_l_post(resCount,1)    = ICa_l_post;
        store.INa_ampa_post(resCount,1) = INa_ampa_post;
        store.ICa_pm_post(resCount,1)   = ICa_pm_post;
        %   Astrocyte
        store.INa_NCX_ast(resCount,1)   = INa_ncx;
        store.ICa_NCX_ast(resCount,1)   = ICa_ncx;
        store.INa_NKA_ast(resCount,1)   = INa_nka_psc;
        store.IK_NKA_ast(resCount,1)    = IK_nka_psc;
        store.IK_ir_ast(resCount,1)     = IK_kir;
        store.INa_l_ast(resCount,1)     = INa_l_psc;
        store.IK_l_ast(resCount,1)      = IK_l_psc;
        store.ICa_l_ast(resCount,1)     = ICa_l_psc;
        store.INa_pf_ast(resCount,1)    = INa_pf;
        store.IK_pf_ast(resCount,1)     = IK_pf;
        store.ICa_pf_ast(resCount,1)    = ICa_pf;
        store.INa_EAAT_ast(resCount,1)  = INaEAAT;
        store.IK_EAAT_ast(resCount,1)   = IKEAAT;
        %   ECS   
        store.INa_dif_ecs(resCount,1)   = INa_dif_ecs;
        store.IK_dif_ecs(resCount,1)    = IK_dif_ecs;
        store.ICa_dif_ecs(resCount,1)   = ICa_dif_ecs;
        % Other variables
        store.sim_time(resCount,1)      = i * dt;
        store.I_ext(resCount,1)         = I_stim; %store external stimulus
        store.glu_rel(resCount,1)       = GluRel;
        % Increment and reset counters
        saveCount                       = 0;
        resCount                        = resCount + 1;
    end
    store.step(i,1)     = i;% store step number
    saveCount           = saveCount + 1; % Increment save counter
    
end % End simulation

%% Plot results
figs = Plot_Model_Results_Full(store,figureFile,iStop);
%% Save workspace and figures
%save(saveFile);
savefig(figs,figureFile);
%% Display elapsed time
clc; % clear command window
fprintf(['Presimulation took ',num2str(elapsed), ' seconds']);
elapsed = toc; % get elapsed time
fprintf(['\nSimulation took ',num2str(elapsed/60), ' minutes']);

%% Helper functions
function [cell]         = Model_Morphology
%Model_Morphology Calculate the cellular morpholgy of the model.
% The volume and surface areas of each cellular compartment in the model is calculated.
% Outputs:
%   1. cell - Structure containing cellular morphology

% Process morphology
% Diameter (m)
processDiameter     = 100e-9;
% Radius (m)
processRadius       = processDiameter / 2;
% Length (m)
cell.processLength  = 2e-6;
% Cross sectional area (m^2)
cell.processCSA     = pi * processRadius^2;
cell.processVol = pi * processRadius^2 * cell.processLength * 1000;
% Perisynaptic cradle morphology
% Diameter outside (m)
diameterOutPsc      = 500e-9;
% Diameter inside (m)
diameterInPsc       = 300e-9;
% Radius outside (m)
radiusOutPsc        = diameterOutPsc / 2;
% Radius inside (m)
radiusInPsc         = diameterInPsc / 2;
% Length (m)
cell.lengthPsc           = 500e-9;
% Half cross sectional area (m^2)
csaPsc              = (pi * radiusInPsc^2) / 2;
% Half surface area (m^2)
cell.pscSA         = (2 * pi * radiusInPsc * cell.lengthPsc) / 2;
% Volume (L)
volOutPsc           = pi * radiusOutPsc^2 * cell.lengthPsc * 1000;
volInPsc            = pi * radiusInPsc^2 * cell.lengthPsc * 1000;
% Half volume of PsC
cell.pscVol        = (volOutPsc - volInPsc) / 2;
% Synapse morphology
% Diameter (m)
diameterSyn         = diameterInPsc - 30e-9;
% Radius (m)
radiusSyn           = diameterSyn / 2;
% Length (m)
lengthSyn           = cell.lengthPsc;
% Half cross sectional area (m^2)
csaSyn              = (pi * radiusSyn^2) / 2;
% Half surface area (m^2)
cell.synSA         = (2 * pi * radiusSyn * lengthSyn) / 2;
% Half volume (L)
cell.volSyn             = (pi * radiusSyn^2 * lengthSyn * 1000)/2;
% Extracellular space (ECS) morphology
% Volume (L)
cell.ecsVol        = (volInPsc - (cell.volSyn*2));
% Diffusion surface area (m^2)
cell.ecsDiffSA     = cell.lengthPsc * ((diameterInPsc - diameterSyn) + (csaPsc - csaSyn));
% Cleft
lCleft = 20e-9;
rCleft = radiusSyn;
cell.cleftVol = lCleft * rCleft^2 * pi * 1000;
cell.ecsVol = cell.ecsVol + cell.cleftVol;

% Volume ratio
cell.ecsVolRatio   = cell.ecsVol / cell.pscVol;
cell.pscVolRatio   = cell.pscVol / cell.ecsVol;
end
function [storage]      = Model_Storage(steps,save_period)
%MODEL_STORAGE Used to create storage variables for model
%   Used to create storage variables for model
%Input:
%   steps - Number of iteration steps
%   save_period - Save every x steps
%Output: 
%   storage - Structure containing storage vectors

temp = zeros(floor(steps/save_period),1); % Temporary vector
% temp1 = zeros(floor(steps),1); % Full amount of steps
% Simulation variables
storage.sim_time    = temp;
storage.I_ext       = temp;
% storage.step        = temp1;
% storage.randNo      = temp;
% storage.meanGlu     = temp;
% Voltage
storage.V_pre       = temp;
storage.V_post      = temp;
storage.V_ast       = temp;
% Concentrations
storage.Na_psc      = temp;
storage.Na_ecs      = temp;
storage.Na_pre      = temp;
storage.Na_post     = temp;
storage.K_psc       = temp;
storage.K_ecs       = temp;
storage.K_pre       = temp;
storage.K_post      = temp;
storage.Ca_pre      = temp;
storage.Ca_post     = temp;
storage.Ca_ecs      = temp;
storage.Ca_psc      = temp;
storage.Glu_ecs     = temp;
% storage.Glu_ecs     = temp1;
% Currents
% Presynapse
storage.INa_vg_pre  = temp;
storage.INa_l_pre   = temp;
storage.INa_NKA_pre = temp;
storage.IK_vg_pre   = temp;
storage.IK_l_pre    = temp;
storage.IK_NKA_pre  = temp;
storage.ICa_vg_pre  = temp;
storage.ICa_l_pre   = temp;
storage.ICa_pm_pre  = temp;
% Postsynapse
storage.INa_nmda_post = temp;
storage.INa_l_post  = temp;
storage.INa_nka_post = temp;
storage.IK_nka_post  = temp;
storage.IK_l_post  = temp;
storage.IK_nmda_post = temp;
storage.ICa_nmda_post = temp;
storage.ICa_l_post  = temp;
storage.INa_ampa_post = temp;
% storage.IK_ampa_post = temp;
storage.ICa_pm_post  = temp;
% Astrocytic
storage.INa_pf_ast  = temp;
storage.INa_l_ast   = temp;
storage.INa_NKA_ast = temp;
storage.INa_EAAT_ast= temp;
storage.INa_NCX_ast = temp;
storage.IK_ir_ast   = temp;
storage.IK_l_ast    = temp;
storage.IK_NKA_ast  = temp;
storage.IK_pf_ast   = temp;
storage.IK_EAAT_ast = temp;
storage.ICa_NCX_ast = temp;
storage.ICa_l_ast   = temp;
% storage.IPMCA_ast   = temp;
storage.ICa_pf_ast  = temp;
% ECS
storage.INa_dif_ecs = temp;
storage.IK_dif_ecs = temp;
storage.ICa_dif_ecs = temp;
% Other
storage.glu_rel     = temp;
storage.prob_rel    = temp;
end
function [dV,x,INa,IK,ICa,Glu] = Neuron_Model(V,Vrest,x,dt,Iext,I,E_Ca)
% Params
mult = 1000;
conv = 1e-2; % Conversion factor from uA/cm^2 to A/m^2
V = V * mult; % Convert voltage to mV
Vrest = Vrest * mult; % Convert voltage to mV
dt = dt * mult; % Convert dt to ms
I = I / conv;% convert I to uA/cm^2
% constants
C = 1; % uF/cm^2
E_Na = 115 + Vrest; % mV
% E_Ca = 137;
% E_Ca = E(3);
E_K = -12 + Vrest; %mV
E_Leak = 10.6 + Vrest; % mV
g_Na = 120; % mS/cm^2
g_K = 36; % mS/cm^2
g_Leak = 0.3; % mS/cm^2
g_Ca = 0.1; % maximal conductance (mS/cm^2);
% Gating variable rates
[alphaN, betaN] = n_equations(V, Vrest);
[alphaM, betaM] = m_equations(V, Vrest);
[alphaH, betaH] = h_equations(V, Vrest);
[m_inf, tau_m] = mca_equations(V);
% Conductances
gK = g_K * x(2)^4;
gNa = g_Na * x(1)^3 * x(3);
gCa = g_Ca * x(4)^3;
% Currents
INa = gNa*(V-E_Na);
IK = gK*(V-E_K);
ICa = gCa*(V-E_Ca);
ILeak = g_Leak*(V-E_Leak);
% Ionic = Iext - (INa+IK+ILeak+ICa+I);
Ionic = Iext - (INa+IK+ILeak+ICa);
% Calculate new membrane potential
dV = V + dt * Ionic * (1/C);
% Calculate new gating variables
x(1) = x(1) + (alphaM *(1-x(1)) - betaM * x(1))*dt;
x(2) = x(2) + (alphaN *(1-x(2)) - betaN * x(2))*dt;
x(3) = x(3) + (alphaH *(1-x(3)) - betaH * x(3))*dt;
x(4) = x(4) + ((m_inf - x(4))/tau_m) * dt;
% Get glutamate concentration
Glu = Trans_Rel(V);
% Conversions
dV = dV / mult; % to V
INa = INa * conv * 1e-2; % to A/m^2
IK = IK * conv * 1e-2; % to A/m^2
ICa = ICa * conv; % to A/m^2
end
function [a_m, b_m]     = m_equations(V,Vrest)
% calculate alpha m and beta m
a_m = (2.5-0.1*(V-Vrest))/(exp(2.5-0.1*(V-Vrest))-1);
b_m = 4*exp((Vrest-V)/18);
end 
function [a_n, b_n]     = n_equations(V,Vrest)
% calculate alpha n and beta n
a_n = (0.1-0.01*(V-Vrest))/(exp(1-0.1*(V-Vrest))-1);
b_n = 0.125*exp((Vrest-V)/80);
end
function [a_h, b_h]     = h_equations(V,Vrest)
% calculate alpha h and beta h 
a_h = 0.07*exp((Vrest-V)/20);
b_h = 1/(1+exp(3-0.1*(V-Vrest)));
end
function [m_inf, tau_m] = mca_equations(V)
m_inf = 1 / (1 + exp(-((V+24.758)/8.429)));
if V >= -40
    tau_m = 0.2702+1.1622*exp(-(V+22.098)^2/164.19);
else
    tau_m = 0.6923*exp((V-4.7)/1089.372);
end
end
function [Trans]        = Trans_Rel(Vpre)
%TRANS_REL Function relating presynaptic voltage with neurotransmitter
%release
% Ref: Synthesis of Models for Excitable Membranes, Synaptic Transmission..
%      by Destexhe, Mainen and Sejnowski (1994)
%Input:
%   Vpre - Presynaptic membrane potential
%Output:
%   Trans - Neurotransmitter concentration

% Params
Tmax = 1e-3; % Maximal transmitter concentration; unit: M
Vp = 2e-6; % Value when half-activation reached; unit: V
Kp = 5e-6; % Slope or steepness; unit: V
% Equation
Trans = Tmax / (1 + exp(-(Vpre-Vp)/Kp));
end
function [INa,IK]       = NKA(Kout,Nain,astrocyte)
%NKA Find current density through NKA
%Input:
%   Kout - Extracellular K concentration
%   Nain - Intracellular Na concentration
%Output:
%   INa - Na current density
%   IK - K current density

% Params
density = 3; % density multiplier
Pmax = 1.12e-6; % Maximal pump velocity; unit: mol/m^2s
Pmax = Pmax * 96485; % Pmax * F; unit: A/m^2
KdNa = 10e-3; % Half maximal velocity concentration for Na; unit: M
% KdK = 1.5e-3; % Half velocity concentration for K; unit: M
if astrocyte == 1
    KdK = 3.6e-3; % Half velocity concentration for K in PsC; unit: M
else
    KdK = 0.6e-3; % Half velocity concentration for K in Pre; unit: M
end
T1    = (Nain^1.5)/((Nain^1.5) + (KdNa^1.5));
T2    = (Kout) / (Kout + KdK);
% Find currents
if astrocyte == 1
    INKA = Pmax * T1 * T2 * density; % Higher density of NKA in PsC
else
    INKA = Pmax * T1 * T2;
end
% Stoichometry
INa   = 3 * INKA;  
IK    = -2 * INKA;
end
function [dateString]   = Get_Date_String
%Get_Date_String Function to return the current date in a formatted String
%Output:
%   dateString - Current date in String format

currentDate = datetime;%get current datetime
currentDate = datestr(currentDate);%convert to string
dateString = replace(currentDate, "-", "_");%remove hypens
dateString = replace(dateString, ":", "_");%remove colons
dateString = replace(dateString, " ", "_");%remove spaces
end
function [E]            = Nernst_Potential(out,in,z)
%Nernst_Potential Calculates the reversal potential of a single ion
%Input: 
%   out - Extracellular ionic concentration
%   in - Intracellular ionic concentration
%   z - Ionic valency
%Output:
%   E - Reversal potential

% Params
R   = 8.3145; % Ideal gas constant; unit: J/K.mol
T   = 310; % Absolute temperature; unit: K
F   = 96485.33; % Faradays constant; unit: C/mol
RT  = R*T; % Energy per mole at given temerature; unit: J/mol
zF  = z*F; % Magnitude of charge per mole;

E = (RT/zF) * log(out/in); % Calculate reversal potential
end
function [ICa]          = PMCA(Ca)
%PMCA Find current density through PMCA
%Input:
%   Ca - Intracellular Ca concentration
%Output:
%   ICa - Ca current density

% Params
Vmax = 0.2e-6; % Maximal pump velocity; unit: mol/m^2s
Vmax = Vmax * 96485; % Vmax * F; unit: A/m^2
Kd   = 0.2e-6; % Half maximal velocity concentration; unit: M
% Calculate current
ICa = Vmax * (Ca/(Kd+Ca));
end
function [INa,ICa]      = NCX(Naout,Nain,Caout,Cain,V)
%NCX Find current density through NCX
%   Ref: Wade et al. 2019 Calcium microdomain formation.
%Input:
%   Naout - Extracellular Na concentration
%   Nain - Intracellular Na concentration
%   Caout - Extracellular Ca concentration
%   Cain - Intracellular Ca concentration
%   V - Membrane potential
%Output:
%   INa - Na current density
%   ICa - Ca current density

% Params
F         = 96485; % Faraday
RT        = 310 * 8.31; % Gas constant * temperature
gamm      = 0.35; % Energy partiton
IbarINaCa = 1; % Maximal current density; unit: A/m^2  1A/m^2 = 100uA/cm^2

Nainout   = (Nain/Naout).^3;
Cainout   = Cain/Caout;
T1        = exp((gamm*F*V)/(RT));
T2        = exp(((gamm-1)*F*V)/(RT));
INCX      = IbarINaCa * ((Nainout *T1) - (Cainout * T2));
% Stoichometry
INa       = 3 * INCX;
ICa       = -INCX;
end
function [I]            = Process(out,in,z,V,V0,l)
%PROCESS Find current along the process
%Input:
%   out - Soma ionic concentration
%   in - PsC ionic concentration
%   z - Valency of ion
%   V - Astrocytic membrane potential
%   V0 - Resting potential
%   l - Length of process
%Output:
%   I - Astrocyte process current

% Params
phi = 10;
epsilon = 8.85e-12 * 0.82;
g = 0.018; % Conductance
e   = 1.602176634e-19; % Elementary charge; unit: C
kB  = 1.380649e-23; % Boltzman constant; unit: J/K
kBT = kB * 310; % unit: J
Vr = Nernst_Potential(out, in,z);
V1 = V - V0 - Vr;

T1 = g *(V1 / l);
T2 = (e * V1) / (l * pi * epsilon) ;
T3 = (phi * kBT / e) - sqrt(T2);
T4 = (-1 * e * T3) / kBT; 

T5 = exp(T4);

I = T1 * T5;
end
function [g]            = Leak_Conduct(I,out,in,V,z)
%LEAK_CONDUCT Calculates the conductance of a simple leak channel
%Input:
%   I - Current
%   out - Extracellular ionic concentration
%   in - Intracellular ionic concentration
%   V - Membrane potential
%   z - Ionic valency
%Output:
%   g - Conductance

E = Nernst_Potential(out,in,z);% Ionic reversal potential
g = -I * (1/(V-E));
end
function [I]            = Leak(g,out,in,V,z)
%LEAK Calculates simple leak current density
%Input:
%   g - Conductance
%   out - Extracellular ionic concentration
%   in - Intracellular ionic concentration
%   V - Membrane potential
%   z - Ionic valency
%Output:
%   I - Leak current

E = Nernst_Potential(out,in,z);
I = g * (V-E);
end
function [IK]           = Kir(out,in,V)
%KIR Calculates the Kir channel current density 
%Input:
%   out - Extracellular K concentration
%   in - Intracellular K concentration
%   V - Membrane potential
%Output:
%   IK - K current density

% Params
gMax = 144; % Maximal conductance; unit: pS/um2
g = gMax * sqrt(out); % Conductance
% E = -0.01; % Reversal potential; unit: V; Ref: DiFranco (2015)
% E = 0.025 * log(out/in); % Reversal potential; unit: V; Ref: Wade NCX (2019)
E = Nernst_Potential(out,in,1); % Reversal potential; unit: V; Ref: DiFranco (2015)
% Calculate current
Ikir = g * (V-E);
IK = Ikir;
end
function [I]            = ECS_Diff(out,in,z)
%ECS_Diff Calculates ionic diffusion from ECS
%Input:
%   out - Extracellular ionic concentration
%   in - Intracellular ionic concentration
%   z - Ionic valency
%Output:
%   I - Current density

% Params
g = 1; % Maximal conductance
lambda = 10; % Conductance multiplier
E = Nernst_Potential(out,in,z); % Reversal potential

I = g * lambda * E; % Current density
end
function [I,dr]         = AMPA(transmitter,r,V,gx)
%AMPA Calculates AMPA current density
%Input:
%   transmitter - Extracellular glu concentration
%   r - Receptor activation
%   V - Membrane potential
%   gx - Conductance multiplier
%Output:
%   I - Current density
%   dr - Rate of activation

reversal    = 0;        % V
alpha       = 1.1e6;   % M-1 sec-1
beta        = 190;      % sec-1
g = 0.26*gx;    % S/m2

dr=alpha * transmitter * (1-r) - (beta * r);    %Receptor Change

I = g * r * (V-reversal);
end
function [INa,ICa,IK,dr]= NMDA(transmitter,r,V,gx)
%NMDA Calculates NMDA current density
%Input:
%   transmitter - Extracellular glu concentration
%   r - Receptor activation
%   V - Membrane potential
%   gx - Conductance multiplier
%Output:
%   INa - Na current density
%   ICa - Ca current density
%   IK - K current density
%   dr - Rate of activation

E_NMDA = 0;
g = 0.18 * gx;
alpha = 7.2e4; % /M.s
beta = 6.6; % /s
dr=alpha*transmitter*(1-r)-beta*r;    %Receptor Change
MG=1;

Mg_V = 1/(1+exp(-0.062*V) * MG/3.57); % Mg block

INMDA = g * r * (V-E_NMDA) * Mg_V;
ICa = INMDA;
INa = 1*INMDA;
IK = -1*INMDA;

% INa = 0;
% ICa = 0;
end
function [INa, IK]      = EAAT2(Glu)
dens=10000*1e-12;%EAAT2 density; per m^2 (10,000 per um^2)
eff=0.5;%EAAT2 efficacy
n=1;%binding sites
F = 96485.33;% Faradays constant; unit: C/mol;
V=30*dens*eff;%max velocity; unit: mol/m^2s
km=20e-6;%concentration half v is reached (accounting for efficacy)

dv=V*(Glu^n/(km^n+Glu^n));%MM
I=dv*F;%current

% g=1;%conductance of single EAAT2
% g_max=dens*g;%maximal conductance
% I=dv*g_max*F*eff;%current

INa=-3*I;%3Na in
IK=I;%1K out
end
function [Glu,Q,R]      = Glu_Rel_Ca(Ca,Glu,Q,R,dt)
%params
k1 = 0.5;%per ms
k_1 = 1;%per ms
k2 = 1;
k3 = 10;%per ms
kT = 5;%per ms
Tmax = 40;%mM (concentration per vesicle)
n = 4;%binding sites
Ca_rest=50e-9;
Ca = (Ca-Ca_rest)*1e6;%convert Ca to uM
dt = dt*1e3;%convert dt to ms 
Glu = Glu*1e3;%convert Glu to mM

X = 1 - Q;
dQ = k1*Ca*X-k_1*Q-n*k2*Q^n; 
dR = k2*Q^n-k3*R;
dT = (Tmax*R)-(Glu*(1/kT));
Q = Q + dt * dQ;
R = R + dt * dR;
Glu = Glu + dt * dT;

Glu = Glu * 1e-3;%convert Glu to M
end

function [fig]          = Plot_Model_Results_Full(store,fileName,iStop)
%PLOT_MODEL_RESULTS Plots model results
%   Plot each result in a separate figure
%Input:
%   store - Structure containing results
%   fileName - File path to save results to
%   iStop - Time at end of plot inset
%Output:
%   fig - Structure containing result figures

import mlreportgen.report.* ;
import mlreportgen.dom.*;
t = store.sim_time; % get time
mM = 1e3;
uM = 1e6;
nM = 1e9;
fA = 1e15;
yLabel = 'I (fA)';
insLen = 100; %inset length; unit: ms

ds = (iStop*1000) -insLen; %inset start
de = ds+1+insLen; % inset end
xTick = [59.9 60];

% Add report container (required)
rpt = Report(fileName,'pdf');
% Add content to container (required)
% Types of content added here: title 
% page and table of contents reporters
titlepg = TitlePage;
titlepg.Title = 'Simulation Results';
titlepg.Author = 'Marinus Toman';
add(rpt,titlepg);
%% Plot astrocyte results
% Astrocytic currents figure
fig(1) = figure('Name','Astrocyte Currents');
fig(1).Units = 'normalized';
fig(1).OuterPosition = [0 0 1 1];

% Column 1 - Sodium
subplot(3,3,1);
plot(t,store.INa_pf_ast*fA);
ylabel(yLabel);
% xlabel('Time (s)');
title('Na^+ Thin Process');
%Inset
axes('Position',[.275 .85 .06 .06]);
plot(t(ds:de),store.INa_pf_ast(ds:de)*fA)
set(gca,'XTick',xTick);

subplot(3,3,4);
plot(t,store.INa_NKA_ast*fA);
ylabel(yLabel);
% xlabel('Time (s)');
title('Na^+ NKA');
%Inset
axes('Position',[.275 .55 .06 .06]);
plot(t(ds:de),store.INa_NKA_ast(ds:de)*fA)
set(gca,'XTick',xTick);

subplot(3,3,7);
plot(t,store.INa_NCX_ast*fA);
ylabel(yLabel);
xlabel('Time (s)');
title('Na^+ NCX');
%Inset
axes('Position',[.275 .25 .06 .06]);
plot(t(ds:de),store.INa_NCX_ast(ds:de)*fA)
set(gca,'XTick',xTick);

% Column 2 - Potassium
subplot(3,3,2);
plot(t,store.IK_pf_ast*fA);
title('K^+ Thin Process');
ylabel(yLabel);
% xlabel('Time (s)');
%Inset
axes('Position',[.56 .85 .06 .06]);
plot(t(ds:de),store.IK_pf_ast(ds:de)*fA)
set(gca,'XTick',xTick);

subplot(3,3,5);
plot(t,store.IK_NKA_ast*fA);
title('K^+ NKA');
ylabel(yLabel);
%xlabel('Time (s)');
%Inset
axes('Position',[.56 .43 .06 .06]);
plot(t(ds:de),store.IK_NKA_ast(ds:de)*fA)
set(gca,'XTick',xTick);

subplot(3,3,8);
plot(t,store.IK_ir_ast*fA);
xlabel('Time (s)');
title('K^+ Kir');
ylabel(yLabel);
%Inset
axes('Position',[.56 .13 .06 .06]);
plot(t(ds:de),store.IK_ir_ast(ds:de)*fA)
set(gca,'XTick',xTick);

% Column 3 - Post Conc
subplot(3,3,3);
plot(t,store.ICa_pf_ast*fA);
title('Ca^{2+} Thin Process');
ylabel(yLabel);
%xlabel('Time (s)');
%Inset
axes('Position',[.84 .85 .06 .06]);
plot(t(ds:de),store.ICa_pf_ast(ds:de)*fA)
set(gca,'XTick',xTick);

% subplot(3,3,6);

subplot(3,3,9);
plot(t,store.ICa_NCX_ast*fA);
xlabel('Time (s)');
title('Ca^{2+} NCX');
ylabel(yLabel);
%Inset
axes('Position',[.84 .13 .06 .06]);
plot(t(ds:de),store.ICa_NCX_ast(ds:de)*fA)
set(gca,'XTick',xTick);

% Add content to report sections (optional)
% Text and formal image added to chapter
chap = Chapter('Astrocyte Currents');
add(chap,Figure(gcf));
add(rpt,chap);

%% Plot Astrocyte and concentrations
fig(2) = figure('Name','Currents and Concentrations');
fig(2).Units = 'normalized';
fig(2).OuterPosition = [0 0 1 1];

% Column 1 - PsC
subplot(3,3,1);
plot(t,store.INa_EAAT_ast*fA);
title('Na^+ EAAT');
ylabel(yLabel);
%xlabel('Time (s)');
%Inset
axes('Position',[.275 .85 .06 .06]);
plot(t(ds:de),store.INa_EAAT_ast(ds:de)*fA)
set(gca,'XTick',xTick);

subplot(3,3,4);
plot(t,store.INa_l_ast*fA);
title('Na^+ Leak PsC');
ylabel(yLabel);
%Inset
axes('Position',[.275 .55 .06 .06]);
plot(t(ds:de),store.INa_l_ast(ds:de)*fA)
set(gca,'XTick',xTick);

subplot(3,3,7);
plot(t,store.ICa_l_ast*fA);
title('Ca Leak PsC');
xlabel('Time (s)');
ylabel(yLabel);
%Inset
axes('Position',[.275 .13 .06 .06]);
plot(t(ds:de),store.ICa_l_ast(ds:de)*fA)
set(gca,'XTick',xTick);

% Col 2
subplot(3,3,2);
plot(t,store.IK_EAAT_ast*fA);
title('K^+ EAAT');
ylabel(yLabel);
%xlabel('Time (s)');
%Inset
axes('Position',[.56 .85 .06 .06]);
plot(t(ds:de),store.IK_EAAT_ast(ds:de)*fA)
set(gca,'XTick',xTick);

subplot(3,3,5);
plot(t,store.IK_l_ast*fA);
title('K Leak PsC');
ylabel(yLabel);
%Inset
axes('Position',[.56 .43 .06 .06]);
plot(t(ds:de),store.IK_l_ast(ds:de)*fA)
set(gca,'XTick',xTick);

subplot(3,3,8);
plot(t,store.Glu_ecs*mM);
xlabel('Time (s)');
ylabel('[Glu] (mM)');
title('[Glu]_{ECS}');
%Inset
axes('Position',[.56 .25 .06 .06]);
plot(t(ds:de),store.Glu_ecs(ds:de)*mM)
set(gca,'XTick',xTick);

% Column 3 - Pre 
subplot(3,3,3);
plot(t,store.Na_post*mM);
title('[Na^{+}]_{Post}');
ylabel('[Na^{+}] (mM)');
% xlabel('Time (s)');
%Inset
axes('Position',[.84 .85 .06 .06]);
plot(t(ds:de),store.Na_post(ds:de)*mM)
set(gca,'XTick',xTick);

subplot(3,3,6);
plot(t,store.K_post*mM);
title('[K^{+}]_{Post}');
ylabel('[K^{+}] (mM)');
% xlabel('Time (s)');
%Inset
axes('Position',[.84 .55 .06 .06]);
plot(t(ds:de),store.K_post(ds:de)*mM)
set(gca,'XTick',xTick);

subplot(3,3,9);
plot(t,store.Ca_post*uM);
xlabel('Time (s)');
ylabel('[Ca^{2+}] (\muM)');
title('[Ca^{2+}]_{Post}');
%Inset
axes('Position',[.84 .25 .06 .06]);
plot(t(ds:de),store.Ca_post(ds:de)*uM)
set(gca,'XTick',xTick);

% Add content to report sections (optional)
% Text and formal image added to chapter
chap = Chapter('Astro Currents & Concentrations');
add(chap,Figure(gcf));
add(rpt,chap);

%% Plot concentrations
fig(3) = figure('Name','Concentrations');
fig(3).Units = 'normalized';
fig(3).OuterPosition = [0 0 1 1];

% Column 1 - PsC
subplot(3,3,1);
plot(t,store.Na_psc*mM);
title('[Na^+]_{PsC}');
ylabel('[Na^+] (mM)');
% xlabel('Time (s)');
%Inset
axes('Position',[.275 .85 .06 .06]);
plot(t(ds:de),store.Na_psc(ds:de)*mM)
set(gca,'XTick',xTick);

subplot(3,3,4);
plot(t,store.K_psc*mM);
ylabel('[K^+] (mM)');
% xlabel('Time (s)');
title('[K^+]_{PsC}');
%Inset
axes('Position',[.275 .55 .06 .06]);
plot(t(ds:de),store.K_psc(ds:de)*mM)
set(gca,'XTick',xTick);

subplot(3,3,7);
plot(t,store.Ca_psc*nM);
xlabel('Time (s)');
ylabel('[Ca^{2+}] (nM)');
title('[Ca^{2+}]_{PsC}');
%Inset
axes('Position',[.275 .25 .06 .06]);
plot(t(ds:de),store.Ca_psc(ds:de)*nM)
set(gca,'XTick',xTick);

% Column 2 - ECS
subplot(3,3,2);
plot(t,store.Na_ecs*mM);
title('[Na^+]_{ECS}');
% xlabel('Time (s)');
ylabel('[Na^+] (mM)');
%Inset
axes('Position',[.56 .73 .06 .06]);
plot(t(ds:de),store.Na_ecs(ds:de)*mM)
set(gca,'XTick',xTick);

subplot(3,3,5);
plot(t,store.K_ecs*mM);
title('[K^+]_{ECS}');
% xlabel('Time (s)');
ylabel('[K^+] (mM)');
%Inset
axes('Position',[.56 .55 .06 .06]);
plot(t(ds:de),store.K_ecs(ds:de)*mM)
set(gca,'XTick',xTick);

subplot(3,3,8);
plot(t,store.Ca_ecs*mM);
xlabel('Time (s)');
ylabel('[Ca^{2+}] (mM)');
title('[Ca^{2+}]_{ECS}');
%Inset
axes('Position',[.56 .13 .06 .06]);
plot(t(ds:de),store.Ca_ecs(ds:de)*mM)
set(gca,'XTick',xTick);

% Column 3 - Pre 
subplot(3,3,3);
plot(t,store.Na_pre*mM);
title('[Na^{+}]_{Pre}');
ylabel('[Na^{+}] (mM)');
% xlabel('Time (s)');
%Inset
axes('Position',[.84 .85 .06 .06]);
plot(t(ds:de),store.Na_pre(ds:de)*mM)
set(gca,'XTick',xTick);

subplot(3,3,6);
plot(t,store.K_pre*mM);
title('[K^{+}]_{Pre}');
ylabel('[K^{+}] (mM)');
% xlabel('Time (s)');
%Inset
axes('Position',[.84 .55 .06 .06]);
plot(t(ds:de),store.K_pre(ds:de)*mM)
set(gca,'XTick',xTick);

subplot(3,3,9);
plot(t,store.Ca_pre*uM);
xlabel('Time (s)');
ylabel('[Ca^{2+}] (\muM)');
title('[Ca^{2+}]_{Pre}');
%Inset
axes('Position',[.84 .25 .06 .06]);
plot(t(ds:de),store.Ca_pre(ds:de)*uM)
set(gca,'XTick',xTick);

% Add content to report sections (optional)
% Text and formal image added to chapter
chap = Chapter('Concentrations');
add(chap,Figure(gcf));
add(rpt,chap);
%% Plot Pre currents 
fig(4) = figure('Name','Presynaptic Currents');
fig(4).Units = 'normalized';
fig(4).OuterPosition = [0 0 1 1];

% Column 1 - Na
subplot(3,3,1);
plot(t,store.INa_vg_pre*fA);
title('Na^+ VG');
ylabel(yLabel);
% xlabel('Time (s)');
%Inset
axes('Position',[.275 .85 .06 .06]);
plot(t(ds:de),store.INa_vg_pre(ds:de)*fA)
set(gca,'XTick',xTick);

subplot(3,3,4);
plot(t,store.INa_NKA_pre*fA);
ylabel(yLabel);
% xlabel('Time (s)');
title('Na^+ NKA');
%Inset
axes('Position',[.275 .55 .06 .06]);
plot(t(ds:de),store.INa_NKA_pre(ds:de)*fA)
set(gca,'XTick',xTick);

subplot(3,3,7);
plot(t,store.INa_l_pre*fA);
xlabel('Time (s)');
ylabel(yLabel);
title('Na^{+} Leak');
%Inset
axes('Position',[.275 .25 .06 .06]);
plot(t(ds:de),store.INa_l_pre(ds:de)*fA)
set(gca,'XTick',xTick);

% Column 2 - K
subplot(3,3,2);
plot(t,store.IK_vg_pre*fA);
title('K^+ VG');
% xlabel('Time (s)');
ylabel(yLabel);
%Inset
axes('Position',[.56 .85 .06 .06]);
plot(t(ds:de),store.IK_vg_pre(ds:de)*fA)
set(gca,'XTick',xTick);

subplot(3,3,5);
plot(t,store.IK_NKA_pre*fA);
title('K^+ NKA');
% xlabel('Time (s)');
ylabel(yLabel);
%Inset
axes('Position',[.56 .43 .06 .06]);
plot(t(ds:de),store.IK_NKA_pre(ds:de)*fA)
set(gca,'XTick',xTick);

subplot(3,3,8);
plot(t,store.IK_l_pre*fA);
xlabel('Time (s)');
ylabel(yLabel);
title('K^{+} Leak');
%Inset
axes('Position',[.56 .25 .06 .06]);
plot(t(ds:de),store.IK_l_pre(ds:de)*fA)
set(gca,'XTick',xTick);

% Column 3 - Ca 
subplot(3,3,3);
plot(t,store.ICa_vg_pre*fA);
title('Ca^{2+} VG');
ylabel(yLabel);
% xlabel('Time (s)');
%Inset
axes('Position',[.84 .85 .06 .06]);
plot(t(ds:de),store.ICa_vg_pre(ds:de)*fA)
set(gca,'XTick',xTick);

subplot(3,3,6);
plot(t,store.ICa_pm_pre*fA);
title('Ca^{2+} PMCA');
ylabel(yLabel);
% xlabel('Time (s)');
%Inset
axes('Position',[.84 .55 .06 .06]);
plot(t(ds:de),store.ICa_pm_pre(ds:de)*fA)
set(gca,'XTick',xTick);

subplot(3,3,9);
plot(t,store.ICa_l_pre*fA);
xlabel('Time (s)');
ylabel(yLabel);
title('Ca^{2+} Leak');
%Inset
axes('Position',[.84 .25 .06 .06]);
plot(t(ds:de),store.ICa_l_pre(ds:de)*fA)
set(gca,'XTick',xTick);

% Add content to report sections (optional)
% Text and formal image added to chapter
chap = Chapter('Presynaptic Currents');
add(chap,Figure(gcf));
add(rpt,chap);

%% Plot Post currents 
fig(5) = figure('Name','Postsynaptic Currents');
fig(5).Units = 'normalized';
fig(5).OuterPosition = [0 0 1 1];

% Column 1 - Na
subplot(3,3,1);
plot(t,store.INa_ampa_post*fA);
title('Na^+ AMPA');
ylabel(yLabel);
% xlabel('Time (s)');
%Inset
axes('Position',[.275 .85 .06 .06]);
plot(t(ds:de),store.INa_ampa_post(ds:de)*fA)
set(gca,'XTick',xTick);

subplot(3,3,4);
plot(t,store.INa_nmda_post*fA);
ylabel(yLabel);
% xlabel('Time (s)');
title('Na^+ NMDA');
%Inset
axes('Position',[.275 .55 .06 .06]);
plot(t(ds:de),store.INa_nmda_post(ds:de)*fA)
set(gca,'XTick',xTick);

subplot(3,3,7);
plot(t,store.INa_nka_post*fA);
xlabel('Time (s)');
ylabel(yLabel);
title('Na^{+} NKA');
%Inset
axes('Position',[.275 .25 .06 .06]);
plot(t(ds:de),store.INa_nka_post(ds:de)*fA)
set(gca,'XTick',xTick);

% Column 2 - K
% subplot(3,3,2);

subplot(3,3,5);
plot(t,store.IK_nmda_post*fA);
title('K^+ NMDA');
% xlabel('Time (s)');
ylabel(yLabel);
%Inset
axes('Position',[.56 .55 .06 .06]);
plot(t(ds:de),store.IK_nmda_post(ds:de)*fA)
set(gca,'XTick',xTick);

subplot(3,3,8);
plot(t,store.IK_nka_post*fA);
xlabel('Time (s)');
ylabel(yLabel);
title('K^{+} NKA');
%Inset
axes('Position',[.56 .13 .06 .06]);
plot(t(ds:de),store.IK_nka_post(ds:de)*fA)
set(gca,'XTick',xTick);

% Column 3 - Ca 
subplot(3,3,3);
plot(t,store.ICa_l_post*fA);
title('Ca^{2+} Leak');
ylabel(yLabel);
% xlabel('Time (s)');
%Inset
axes('Position',[.84 .85 .06 .06]);
plot(t(ds:de),store.ICa_l_post(ds:de)*fA)
set(gca,'XTick',xTick);

subplot(3,3,6);
plot(t,store.ICa_nmda_post*fA);
title('Ca^{2+} NMDA');
ylabel(yLabel);
% xlabel('Time (s)');
%Inset
axes('Position',[.84 .55 .06 .06]);
plot(t(ds:de),store.ICa_nmda_post(ds:de)*fA)
set(gca,'XTick',xTick);

subplot(3,3,9);
plot(t,store.ICa_pm_post*fA);
xlabel('Time (s)');
ylabel(yLabel);
title('Ca^{2+} PMCA');
%Inset
axes('Position',[.84 .25 .06 .06]);
plot(t(ds:de),store.ICa_pm_post(ds:de)*fA)
set(gca,'XTick',xTick);

% Add content to report sections (optional)
% Text and formal image added to chapter
chap = Chapter('Postsynaptic Currents');
add(chap,Figure(gcf));
add(rpt,chap);

%% Plot ECS currents 
fig(6) = figure('Name','ECS');
fig(6).Units = 'normalized';
fig(6).OuterPosition = [0 0 1 1];

% Column 1 - ECS Diffusion
subplot(3,3,1);
plot(t,store.INa_dif_ecs*fA);
title('Na^+ ECS Dif');
ylabel(yLabel);
% xlabel('Time (s)');
%Inset
axes('Position',[.275 .73 .06 .06]);
plot(t(ds:de),store.INa_dif_ecs(ds:de)*fA)
set(gca,'XTick',xTick);

subplot(3,3,4);
plot(t,store.IK_dif_ecs*fA);
ylabel(yLabel);
% xlabel('Time (s)');
title('K^+ ECS Dif');
%Inset
axes('Position',[.275 .55 .06 .06]);
plot(t(ds:de),store.IK_dif_ecs(ds:de)*fA)
set(gca,'XTick',xTick);

subplot(3,3,7);
plot(t,store.ICa_dif_ecs*fA);
xlabel('Time (s)');
ylabel(yLabel);
title('Ca^{2+} ECS Dif');
%Inset
axes('Position',[.275 .13 .06 .06]);
plot(t(ds:de),store.ICa_dif_ecs(ds:de)*fA)
set(gca,'XTick',xTick);

% Column 2 - Post leak
subplot(3,3,2);
plot(t,store.INa_l_post*fA);
title('Na^+ Leak Post');
% xlabel('Time (s)');
ylabel(yLabel);
%Inset
axes('Position',[.56 .85 .06 .06]);
plot(t(ds:de),store.INa_l_post(ds:de)*fA)
set(gca,'XTick',xTick);

subplot(3,3,5);
plot(t,store.IK_l_post*fA);
title('K^+ Leak Post');
% xlabel('Time (s)');
ylabel(yLabel);
%Inset
axes('Position',[.56 .43 .06 .06]);
plot(t(ds:de),store.IK_l_post(ds:de)*fA)
set(gca,'XTick',xTick);

subplot(3,3,8);
plot(t,store.I_ext);
xlabel('Time (s)');
ylabel('I (nA)');
title('External Stimulus');
%Inset
axes('Position',[.56 .25 .06 .06]);
plot(t(ds:de),store.I_ext(ds:de))
set(gca,'XTick',xTick);

% Column 3 - Voltages 
subplot(3,3,3);
plot(t,store.V_pre*mM);
title("Presynaptic Membrane Potential");
% xlabel('Time (s)');
ylabel('V_{Pre} (mV)');
%Inset
axes('Position',[.84 .85 .06 .06]);
plot(t(ds:de),store.V_pre(ds:de)*mM)
set(gca,'XTick',xTick);

subplot(3,3,6);
plot(t,store.V_post*mM);
title("Postsynaptic Membrane Potential");
% xlabel('Time (s)');
ylabel('V_{Post} (mV)');
%Inset
axes('Position',[.84 .55 .06 .06]);
plot(t(ds:de),store.V_post(ds:de)*mM)
set(gca,'XTick',xTick);

subplot(3,3,9);
plot(t,store.V_ast*mM);
title("Astrocytic Potential");
xlabel('Time (s)');
ylabel('V_{Ast} (mV)');
%Inset
axes('Position',[.84 .25 .06 .06]);
plot(t(ds:de),store.V_ast(ds:de)*mM)
set(gca,'XTick',xTick);

% Add content to report sections (optional)
% Text and formal image added to chapter
chap = Chapter('ECS Currents');
add(chap,Figure(gcf));
add(rpt,chap);

% Close the report (required)
close(rpt);
end