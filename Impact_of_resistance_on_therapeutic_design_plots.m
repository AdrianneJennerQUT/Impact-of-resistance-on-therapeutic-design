%% colours:
healthy = "#0071BD"; healthyRGB = [0, 113, 189];
cancer = "#F3AE1A"; cancerRGB = [243, 174, 26];
mutant1 = "#D95318"; mutant1RGB = [217, 83, 24];
mutant2 = "#78AB2F"; mutant2RGB = [120, 171, 47];

treatedcancer = '#1ac8f3'; treatedcancerRGB = [26, 200, 243];
treatedmutant1 = '#b918d9'; treatedmutant1RGB = [185, 24, 217];

%% panel 3
%% 3AB
N = 1000;          % total number of cells at each timestep
M = 10000;          % total number of timesteps
pc = 0.1;         % starting proportion for cancerous cells
Nc = floor(pc*N); % starting number of cancerous cells
Nh = N - Nc;      % starting number of healthy cells
Nc_vec = zeros(M+1,1); Nc_vec(1) = Nc;    % vector containing the number of cells of type 1 over time

fh = 1;         % fitness value for healthy cells
fc = fh;         % fitness value for cancerous cells

for n = 1:M
    Nc_n = Nc; Nh_n = Nh; % store current N1 and N2 for probability calculations
    
    % select a cell type to proliferate based on fitness
    r = rand; % generate random number to be used for probabilities
    if r < Nc_n*fc/(Nc_n*fc+Nh_n*fh) % if random value is less than the probability of cancer proliferating
        Nc = Nc + 1; % increase the number of cells of cancer cells
    else
        Nh = Nh + 1; % increase the number of cells of healthy cells
    end
    
    % select a cell type to die
    r = rand; % generate a new random number for new probabilities
    if r < Nc_n/N % if random value is less than the probability of cancer dying
        Nc = Nc - 1; % decrease the number of cancer cells
    else
        Nh = Nh - 1; % decrease the number of healthy cells
    end
    
    Nc_vec(n+1) = Nc; % input new value to vector
end

figure
hold on
plot(0:M,1-Nc_vec/N,"Color",healthy)
plot(0:M,Nc_vec/N,"Color",cancer)
xlabel('Time (cell divisions)')
ylabel('Proportion of cells')
legend('Healthy cells','Cancerous cells')
xlim([0,M])
ylim([0,1])
set(gca,'FontSize',18)


%% 3C
f1_vec = 0:0.1:10; % fitness values for cell type 1

N = 100;          % total number of cells at each timestep
M = 100;          % total number of timesteps
p1 = 0.5;         % starting proportion for cells of type 1
f2 = 5;           % fitness of type 2 cells
n_total = 10000;   % number of simulations

sims_mat = zeros(n_total, length(f1_vec)); % matrix to store simulation information

progress = 0;
for f1_n = 1:length(f1_vec)
    for n = 1:n_total % run n_total simulations for each f1 value
        N1 = floor(p1*N); % starting number of type 1 cells
        N2 = N - N1;      % starting number of type 2 cells

        for i = 1:M
            N1_i = N1; N2_i = N2; % store current N1 and N2 for probability calculations

            % select a cell type to proliferate based on fitness
            r = rand; % generate random number to be used for probabilities
            if r < N1_i*f1_vec(f1_n)/(N1_i*f1_vec(f1_n)+N2_i*f2) % if random value is less than the probability of type 1 proliferating
                N1 = N1 + 1; % increase the number of cells of type 1
            else
                N2 = N2 + 1; % increase the number of cells of type 2
            end

            % select a cell type to die
            r = rand; % generate a new random number for new probabilities
            if r < N1_i/N % if random value is less than the probability of type 1 dying
                N1 = N1 - 1; % decrease the number of cells of type 1
            else
                N2 = N2 - 1; % decrease the number of cells of type 2
            end
        end
        progress = progress + 1;
        fprintf('%.3f%% complete\n', progress/length(f1_vec)/n_total*100)

        sims_mat(n, f1_n) = N1/N; % update matrix with simulation value for the current f1
    end
end

options.handle = gcf;
options.error = 'std';
options.alpha = 0.5;
options.line_width = 2;

options.color_line = cancerRGB./255;
options.color_area = cancerRGB./255;
plot_areaerrorbar(sims_mat,options)

xlabel('Cancer fitness')
ylabel('Final cancer proportion')
set(gca,'XTick', [20 40 60 80 100])
set(gca,'XTickLabel', [2 4 6 8 10])
xlim([0,length(f1_vec)])
ylim([0,1])
set(gca,'FontSize',18)


%% 3D
N = 100;          % total number of cells at each timestep
M_vec = 65:500;  % total number of timesteps
fc_vec = [0:0.1:20,21:500];   % fitness value for cancer cells, proportional to proliferation rate
fh = 1;               % fitness value for healthy cells
n_total = 1000;   % number of simulations
tumor_thres = 0.5; % threshold for percentage of cancer cells to be considered tumorous

fc_tumorvec = zeros(size(M_vec)); % vector containing the fitnesses that allow for tumor growth for each M

progress = 0;
for M_n = 1:length(M_vec)
    for fc_n = 1:length(fc_vec)
        avg = 0;
        for n = 1:n_total
            Nc = 1;           % number of cancer cells starts at 1
            Nh = N - Nc;      % starting number of healthy cells

            for i = 1:M_vec(M_n)
                Nc_i = Nc; Nh_i = Nh; % store current Nc and Nh for probability calculations

                % select a cell type to proliferate based on fitness
                r = rand; % generate random number to be used for probabilities
                if r < Nc_i*fc_vec(fc_n)/(Nc_i*fc_vec(fc_n)+Nh_i*fh) % if random value is less than the probability of type 1 proliferating
                    Nc = Nc + 1; % increase the number of cells of type 1
                else
                    Nh = Nh + 1; % increase the number of cells of type 2
                end

                % select a cell type to die
                r = rand; % generate a new random number for new probabilities
                if r < Nc_i/N % if random value is less than the probability of type 1 dying
                    Nc = Nc - 1; % decrease the number of cells of type 1
                else
                    Nh = Nh - 1; % decrease the number of cells of type 2
                end
            end
            avg = avg + Nc/n_total; % update average value
        end

        if avg/N > tumor_thres
            fc_tumorvec(M_n) = fc_vec(fc_n); % save the fitness value that allowed the cell population to reached the threshold
            break
        end

        progress = progress + 1;
        fprintf('%.2f%% complete.\n', progress/length(fc_vec)/length(M_vec)*100)
    end
end

fc_tumorvec(fc_tumorvec==0) = NaN; % remove fitnesses that didnt allow population to reach the threshold

figure
plot(M_vec, fc_tumorvec)
%title('Fitness value that allows the growth of a tumor as the total number of cell divisions increase')
xlabel('Total time (cell divisions)')
ylabel('Min. fitness for growth')
xlim([0 M_vec(end)])
set(gca,'FontSize',18)
%ylim([0 10])

axes('position',[.55 .25 .25 .25],'FontSize',15)
box on % put box around new pair of axes
M_range = 150:500;
plot(M_range,fc_tumorvec(M_range-M_vec(1))) % plot on new axes
axis tight



%% panel 4
%% 4A
N = 1000;
M = 20000;
fi = 1.1; % fit cells
fj = 1; % less fit cells
fk = 3; % very fit mutated cells
mk = 1/100; % percentage chance for k mutation given type i reproduces
Ni_0 = 0.1*N;
Nk_0 = 0;
n_total = 1000;

Ni_mat = zeros(n_total, M+1); % matrix to store simulation information
Nk_mat = Ni_mat;
tunnel_percent = 0; % initially 0% tunnel

for i = 1:n_total
    Ni = Ni_0;
    Nj = N - Ni_0 - Nk_0;
    Nk = Nk_0;

    Ni_mat(i,1) = Ni/N;
    Nk_mat(i,1) = Nk/N;
    
    for n = 1:M
        % store current values for probability calcs
        Ni_n = Ni; Nj_n = Nj; Nk_n = Nk;

        r = rand;
        if r < Ni_n*fi/(Ni_n*fi+Nj_n*fj+Nk_n*fk)
            if rand < mk
                Nk = Nk + 1; % mutate if random value is less than mk
            else
                Ni = Ni + 1;
            end
        elseif r >= Ni_n*fi/(Ni_n*fi+Nj_n*fj+Nk_n*fk) && r < (Ni_n*fi+Nj_n*fj)/(Ni_n*fi+Nj_n*fj+Nk_n*fk)
            Nj = Nj + 1;
        else
            Nk = Nk + 1;
        end

        r = rand;
        if r < Ni_n/N
            Ni = Ni - 1;
        elseif r >= Ni_n/N && r < (Ni_n+Nj_n)/N
            Nj = Nj - 1;
        else
            Nk = Nk - 1;
        end

        Ni_mat(i,n+1) = Ni/N; Nk_mat(i,n+1) = Nk/N;
    end
end

options.handle = gcf;
options.error = 'std';
options.alpha = 0.5;
options.line_width = 2;

options.color_line = cancerRGB./255;
options.color_area = cancerRGB./255;
plot_areaerrorbar(Ni_mat,options)

hold on
options.color_line = mutant1RGB./255;
options.color_area = mutant1RGB./255;
plot_areaerrorbar(Nk_mat,options)

hold on
options.color_line = healthyRGB./255;
options.color_area = healthyRGB./255;
plot_areaerrorbar(1-Ni_mat-Nk_mat,options)

xlabel('Time (cell divisions)')
ylabel('Proportion of cells')
xlim([0,M])
ylim([0,1])
set(gca,'FontSize',18)

%% 4B
N = 1000;
M = 20000;
fj = 1;
mk = 1/1000;
Ni_0 = 0.1*N;
Nk_0 = 0;
n_total = 100;

fi_vec = 0:0.01:5;
fk_vec = 0:0.01:5;

fix_i_mat = zeros(length(fi_vec),length(fk_vec)); % matrix for % of the time that type i fixates for each fi,fk
fix_k_mat = fix_i_mat;

progress = 0;

for fi_n = 1:length(fi_vec)
    for fk_n = 1:length(fk_vec)
        fix_i = 0; % number of times type i fixates
        fix_k = 0;
        fi = fi_vec(fi_n); fk = fk_vec(fk_n);

        for i = 1:n_total
            Ni = Ni_0;
            Nj = N - Ni_0 - Nk_0;
            Nk = Nk_0;
            
            for n = 1:M
                % store current values for probability calcs
                Ni_n = Ni; Nj_n = Nj; Nk_n = Nk;
        
                r = rand;
                if r < Ni_n*fi/(Ni_n*fi+Nj_n*fj+Nk_n*fk)
                    if rand < mk
                        Nk = Nk + 1; % mutate if random value is less than mk
                    else
                        Ni = Ni + 1;
                    end
                elseif r >= Ni_n*fi/(Ni_n*fi+Nj_n*fj+Nk_n*fk) && r < (Ni_n*fi+Nj_n*fj)/(Ni_n*fi+Nj_n*fj+Nk_n*fk)
                    Nj = Nj + 1;
                else
                    Nk = Nk + 1;
                end
        
                r = rand;
                if r < Ni_n/N
                    Ni = Ni - 1;
                elseif r >= Ni_n/N && r < (Ni_n+Nj_n)/N
                    Nj = Nj - 1;
                else
                    Nk = Nk - 1;
                end

                if Ni == N
                    fix_i = fix_i + 1; 
                    break
                end
                if Nk == N
                    fix_k = fix_k + 1; 
                    break
                end
            end
            progress = progress + 1;
            fprintf('%.3f%% complete\n',progress/n_total/length(fi_vec)/length(fk_vec)*100)
        end
        fix_i_mat(fi_n,fk_n) = fix_i/n_total;
        fix_k_mat(fi_n,fk_n) = fix_k/n_total;
    end
end

% save
save('fix_i_mat.mat','fix_i_mat'); save('fix_k_mat.mat','fix_k_mat'); 

%% load
fi_vec = 0:0.01:5;
fk_vec = 0:0.01:5;
fix_i_mat_file = matfile('fix_i_mat.mat'); fix_i_mat = fix_i_mat_file.fix_i_mat;
fix_k_mat_file = matfile('fix_k_mat.mat'); fix_k_mat = fix_k_mat_file.fix_k_mat;

% plot
RGBpic = zeros(length(fi_vec),length(fk_vec),3);
RGBpic(:,:,1) = (fix_i_mat == fix_k_mat)*healthyRGB(1) + (fix_i_mat > fix_k_mat)*cancerRGB(1) + (fix_i_mat < fix_k_mat)*mutant1RGB(1);
RGBpic(:,:,2) = (fix_i_mat == fix_k_mat)*healthyRGB(2) + (fix_i_mat > fix_k_mat)*cancerRGB(2) + (fix_i_mat < fix_k_mat)*mutant1RGB(2);
RGBpic(:,:,3) = (fix_i_mat == fix_k_mat)*healthyRGB(3) + (fix_i_mat > fix_k_mat)*cancerRGB(3) + (fix_i_mat < fix_k_mat)*mutant1RGB(3);


figure
imshow(flipud(RGBpic)./255)
axis on
xlabel('Mutant cancer fitness')
ylabel('Original cancer fitness')
xticks([1 100 200 300 400 500]);
yticks([1 100 200 300 400 500]);
set(gca,'XTick', [1 100 200 300 400 500])
set(gca,'XTickLabel', [0 1 2 3 4 5])
set(gca,'YTick', [1 100 200 300 400])
set(gca,'YTickLabel', [5 4 3 2 1])
set(gca,'FontSize',18)



%% 4DE
% parameters
N = 1000;
M = 20000;
fi = 1.1; % fit cells
fj = 1; % less fit cells
fk = 1.5; % very fit mutated cells
mk = 1/100; % percentage chance for k mutation given type i reproduces
Ni_0 = 0.1*N;
Nk_0 = 0;
eff_i = 0.8;
eff_k = 0.4;
t_start = 5000;
t_length = 3000; % 0 for 4D, 5000 for 4E

Ni = Ni_0;
Nj = N - Ni_0 - Nk_0;
Nk = Nk_0;

Ni_vec = zeros(1,M+1);
Nk_vec = zeros(1,M+1);
% simulation
for n = 1:M
    % store current values for probability calcs
    Ni_n = Ni; Nj_n = Nj; Nk_n = Nk;

    if n >= t_start && n < t_start+t_length
        efi = fi*(1-eff_i);
        efk = fk*(1-eff_k);
    else
        efi = fi;
        efk = fk;
    end
    
    r = rand;
    if r < Ni_n*efi/(Ni_n*efi+Nj_n*fj+Nk_n*efk)
        if rand < mk
            Nk = Nk + 1; % mutate if random value is less than mk
        else
            Ni = Ni + 1;
        end
    elseif r >= Ni_n*efi/(Ni_n*efi+Nj_n*fj+Nk_n*efk) && r < (Ni_n*efi+Nj_n*fj)/(Ni_n*efi+Nj_n*fj+Nk_n*efk)
        Nj = Nj + 1;
    else
        Nk = Nk + 1;
    end
    
    r = rand;
    if r < Ni_n/N
        Ni = Ni - 1;
    elseif r >= Ni_n/N && r < (Ni_n+Nj_n)/N
        Nj = Nj - 1;
    else
        Nk = Nk - 1;
    end
    
    Ni_vec(n+1) = Ni; Nk_vec(n+1) = Nk;
end
Ni_vec = Ni_vec/N; Nk_vec = Nk_vec/N;

figure
a = area(0:M, [Nk_vec', Ni_vec', (1 - Ni_vec - Nk_vec)']);
a(1).FaceColor = mutant1;
a(2).FaceColor = cancer;
a(3).FaceColor = healthy;
if t_length ~= 0
    hold on
    xline(t_start+1,'k--')
    xline(t_start+t_length+1,'k--')
end
xlabel('Time (cell divisions)')
ylabel('Proportion of population')
%legend('Very fit mutated cell', 'Fit cells', 'Unfit cells')
ylim([0 1])
set(gca,'FontSize',18)

%% 4F (violin plot from final results of DE)
% parameters
N = 1000;
M = 20000;
fi = 1.1; % fit cells
fj = 1; % less fit cells
fk = 1.5; % very fit mutated cells
mk = 1/100; % percentage chance for k mutation given type i reproduces
Ni_0 = 0.1*N;
Nk_0 = 0;
eff_i = 0.8;
eff_k = 0.4;
t_start = 5000;
t_length = 3000;

n_total = 1000;

Ni_end = zeros(1,n_total);
Nk_end = Ni_end;
Ni_t_end = Ni_end;
Nk_t_end = Ni_end;

for i = 1:n_total
    Ni = Ni_0;
    Nj = N - Ni_0 - Nk_0;
    Nk = Nk_0;
    Ni_t = Ni;
    Nj_t = Nj;
    Nk_t = Nk;
    % simulation
    for n = 1:M
        % store current values for probability calcs
        Ni_n = Ni; Nj_n = Nj; Nk_n = Nk;
        Ni_t_n = Ni_t; Nj_t_n = Nj_t; Nk_t_n = Nk_t;
    
        if n >= t_start && n < t_start+t_length
            efi = fi*(1-eff_i);
            efk = fk*(1-eff_k);
        else
            efi = fi;
            efk = fk;
        end
        
        % no treatment sim
        r = rand;
        if r < Ni_n*fi/(Ni_n*fi+Nj_n*fj+Nk_n*fk)
            if rand < mk
                Nk = Nk + 1;
            else
                Ni = Ni + 1;
            end
        elseif r >= Ni_n*fi/(Ni_n*fi+Nj_n*fj+Nk_n*fk) && r < (Ni_n*fi+Nj_n*fj)/(Ni_n*fi+Nj_n*fj+Nk_n*fk)
            Nj = Nj + 1;
        else
            Nk = Nk + 1;
        end
        r = rand;
        if r < Ni_n/N
            Ni = Ni - 1;
        elseif r >= Ni_n/N && r < (Ni_n+Nj_n)/N
            Nj = Nj - 1;
        else
            Nk = Nk - 1;
        end

        % treatment sim
        r = rand;
        if r < Ni_t_n*efi/(Ni_t_n*efi+Nj_t_n*fj+Nk_t_n*efk)
            if rand < mk
                Nk_t = Nk_t + 1; % mutate if random value is less than mk
            else
                Ni_t = Ni_t + 1;
            end
        elseif r >= Ni_t_n*efi/(Ni_t_n*efi+Nj_t_n*fj+Nk_t_n*efk) && r < (Ni_t_n*efi+Nj_t_n*fj)/(Ni_t_n*efi+Nj_t_n*fj+Nk_t_n*efk)
            Nj_t = Nj_t + 1;
        else
            Nk_t = Nk_t + 1;
        end
        r = rand;
        if r < Ni_t_n/N
            Ni_t = Ni_t - 1;
        elseif r >= Ni_t_n/N && r < (Ni_t_n+Nj_t_n)/N
            Nj_t = Nj_t - 1;
        else
            Nk_t = Nk_t - 1;
        end
    end
    Ni_end(i) = Ni/N; Nk_end(i) = Nk/N;
    Ni_t_end(i) = Ni_t/N; Nk_t_end(i) = Nk_t/N;
end

figure
vs = violinplot([Ni_end',Nk_end',Ni_t_end',Nk_t_end'],{'Untreated cancer','Untreated mutant','Treated cancer','Treated mutant'});
xlabel('Cancer type')
ylabel('Proportion of population')
set(gca,'FontSize',18)




%% panel 5
%% 5B-F, S3A-E,S4AB
% A: t_start = 5000, B: t_start = 2000;, C: t_start2 =
% t_start+t_length+5000, D: st_weakness = 5, E: drug_conc = @(n,ts) exp(log(0.5)*time_per_n/half_life*(n-ts)) (half_life = 5000)
% 3-cell sim
N = 1000;
M = 40000;
fi = 1.1; % fit cells
fj = 1; % less fit cells
fk = 1.5; % very fit mutated cells
mk = 1/100; % percentage chance for k mutation given type i reproduces
Ni_0 = 0.1*N;
Nk_0 = 0;

t_start = 5000; % starting time for treatment
t_length = 5000; % number of cell divisions to apply treatment over
t_start2 = 0;%t_start + t_length + 5000;
eff_i = 0.8; % how effective treatment is on cell type i (0 to 1)
eff_k = 0.4; % how effective treatment is on cell type k (0 to 1)
eff_i2 = 0.4;
eff_k2 = 0.8;
st_weakness = 1; % how much weaker sustained treatment dosage is to normal

time_per_n = 1; % t/n, rate of time per cell division n
half_life = 5000; % t, same time units as in time_per_n
drug_conc = @(n,ts) 1;%exp(log(0.5)*time_per_n/half_life*(n-ts)); % drug concentration as a ratio of its initial amount, with starting time ts
t_off = 0; % drug concentration that treatment is considered 'turned off'
plot_conc = true; % true if plotting drug concentration

n_total = 1000;

Ni_end = zeros(1,n_total);
Nk_end = Ni_end;
Ni_t_end = Ni_end;
Nk_t_end = Ni_end;

Ni_mat = zeros(n_total, M+1); % matrix to store simulation information
Nk_mat = Ni_mat;

progress = 0;
for i = 1:n_total
    Ni = Ni_0; Ni_t = Ni;
    Nj = N - Ni_0 - Nk_0; Nj_t = Nj;
    Nk = Nk_0; Nk_t = Nk;
    
    for n = 1:M
        if n >= t_start && n < t_start+t_length*st_weakness && drug_conc(n,t_start) == 1
            efi = fi*(1-eff_i/st_weakness);
            efk = fk*(1-eff_k/st_weakness);
        elseif t_start2 > t_start && n >= t_start2 && n < t_start2+t_length*st_weakness && drug_conc(n,t_start) == 1
            efi = fi*(1-eff_i2/st_weakness);
            efk = fk*(1-eff_k2/st_weakness);
        elseif drug_conc(n,t_start) < 1 && drug_conc(n,t_start) > t_off
            efi = fi*(1-eff_i*drug_conc(n,t_start));
            efk = fk*(1-eff_k*drug_conc(n,t_start));
        else
            efi = fi;
            efk = fk;
        end
    
        % store current values for probability calcs
        Ni_n = Ni; Nj_n = Nj; Nk_n = Nk;
        Ni_t_n = Ni_t; Nj_t_n = Nj_t; Nk_t_n = Nk_t;
    
        % without treatment
        r = rand;
        if r < Ni_n*fi/(Ni_n*fi+Nj_n*fj+Nk_n*fk)
            if rand < mk
                Nk = Nk + 1; % mutate if random value is less than mk
            else
                Ni = Ni + 1;
            end
        elseif r >= Ni_n*fi/(Ni_n*fi+Nj_n*fj+Nk_n*fk) && r < (Ni_n*fi+Nj_n*fj)/(Ni_n*fi+Nj_n*fj+Nk_n*fk)
            Nj = Nj + 1;
        else
            Nk = Nk + 1;
        end
    
        r = rand;
        if r < Ni_n/N
            Ni = Ni - 1;
        elseif r >= Ni_n/N && r < (Ni_n+Nj_n)/N
            Nj = Nj - 1;
        else
            Nk = Nk - 1;
        end

        % with treatment
        r = rand;
        if r < Ni_t_n*efi/(Ni_t_n*efi+Nj_t_n*fj+Nk_t_n*efk)
            if rand < mk
                Nk_t = Nk_t + 1; % mutate if random value is less than mk
            else
                Ni_t = Ni_t + 1;
            end
        elseif r >= Ni_t_n*efi/(Ni_t_n*efi+Nj_t_n*fj+Nk_t_n*efk) && r < (Ni_t_n*efi+Nj_t_n*fj)/(Ni_t_n*efi+Nj_t_n*fj+Nk_t_n*efk)
            Nj_t = Nj_t + 1;
        else
            Nk_t = Nk_t + 1;
        end
    
        r = rand;
        if r < Ni_t_n/N
            Ni_t = Ni_t - 1;
        elseif r >= Ni_t_n/N && r < (Ni_t_n+Nj_t_n)/N
            Nj_t = Nj_t - 1;
        else
            Nk_t = Nk_t - 1;
        end
        
        Ni_mat(i,n+1) = Ni_t/N; Nk_mat(i,n+1) = Nk_t/N;
    end

    Ni_end(i) = Ni/N; Nk_end(i) = Nk/N;
    Ni_t_end(i) = Ni_t/N; Nk_t_end(i) = Nk_t/N;

    progress = progress + 1;
    fprintf('%.2f%% complete\n', progress/n_total*100)
end


% plot

figure
options.handle = gcf;
options.error = 'std';
options.alpha = 0.5;
options.line_width = 2;

options.color_line = cancerRGB./255;
options.color_area = cancerRGB./255;
plot_areaerrorbar(Ni_mat,options)

hold on
options.color_line = mutant1RGB./255;
options.color_area = mutant1RGB./255;
plot_areaerrorbar(Nk_mat,options)

hold on
options.color_line = healthyRGB./255;
options.color_area = healthyRGB./255;
plot_areaerrorbar(1-Ni_mat-Nk_mat,options)

hold on
xline(t_start+1,'k--')
if drug_conc(t_start+1000,t_start) == 1
    xline(t_start+t_length*st_weakness+1,'k--')
end
if t_start2 > t_start
    xline(t_start2+1,'k--')
    xline(t_start2+t_length*st_weakness+1,'k--')
end

% plotting drug concentration
if plot_conc && drug_conc(t_start+1000,t_start) ~= 1
    plot([0,t_start,t_start:0.1:min(log(t_off)*half_life/time_per_n/log(0.5),M),M],[1,1,drug_conc(t_start:0.1:min(log(t_off)*half_life/time_per_n/log(0.5),M),t_start),0],'b--')
end

xlabel('Time (cell divisions)')
ylabel('Proportion of cells')
xlim([0,M])
ylim([0,1])
set(gca,'FontSize',18)



figure
plot(Nk_mat',"Color",mutant1)
hold on
plot(Ni_mat',"Color",cancer)
xline(t_start+1,'k--')
if drug_conc(t_start+1000,t_start) == 1
    xline(t_start+t_length*st_weakness+1,'k--')
end
if t_start2 > t_start
    xline(t_start2+1,'k--')
    xline(t_start2+t_length*st_weakness+1,'k--')
end
if plot_conc && drug_conc(t_start+1000,t_start) ~= 1
    plot([0,t_start,t_start:0.1:min(log(t_off)*half_life/time_per_n/log(0.5),M),M],[1,1,drug_conc(t_start:0.1:min(log(t_off)*half_life/time_per_n/log(0.5),M),t_start),0],'b--')
end
xlabel('Time (cell divisions)')
ylabel('Proportion of cells')
xlim([0,M])
ylim([0,1])
set(gca,'FontSize',18)



figure
vs = violinplot([Ni_end',Nk_end',Ni_t_end',Nk_t_end'],{'Untreated cancer','Untreated mutant','Treated cancer','Treated mutant'});
xlabel('Cancer type')
ylabel('Proportion of population')
set(gca,'FontSize',18)


%% 5GH, S2BE
percent_ifile = matfile('percent_i_test2.mat'); percent_i = percent_ifile.percent_i;
percent_kfile = matfile('percent_k_test2.mat'); percent_k = percent_kfile.percent_k;

mk_vec = 0:2/10000:1/100; % range of mutation rates
eff_k_vec = 0:0.01:1; % range of efficacy of sustained treatment on mutated cancer

figure
surf(eff_k_vec, mk_vec, percent_i)
colormap turbo
shading interp
% colorbar
clim([0 1])
%title(['Chance for type i cancer to reach ', num2str(T*100), '% of N'])
xlabel('Efficacy against mutant cells')
ylabel('Mutation rate')
zlabel('Probability')
set(gca,'FontSize',18)

figure
surf(eff_k_vec, mk_vec, percent_k)
colormap turbo
shading interp
% colorbar
clim([0 1])
%title(['Chance for type k cancer to reach ', num2str(T*100), '% of N'])
xlabel('Efficacy against mutant cells')
ylabel('Mutation rate')
zlabel('Probability')
set(gca,'FontSize',18)

percent_ifile = matfile('percent_i_normaltest2.mat'); percent_i = percent_ifile.percent_i;
percent_kfile = matfile('percent_k_normaltest2.mat'); percent_k = percent_kfile.percent_k;

figure
surf(eff_k_vec, mk_vec, percent_i)
colormap turbo%summer
shading interp
%colorbar
clim([0 1])
%title(['Chance for type i cancer to reach ', num2str(T*100), '% of N'])
xlabel('Efficacy against mutant cells')
ylabel('Mutation rate')
zlabel('Probability')
set(gca,'FontSize',18)

figure
surf(eff_k_vec, mk_vec, percent_k)
colormap turbo%summer
shading interp
colorbar
clim([0 1])
%title(['Chance for type k cancer to reach ', num2str(T*100), '% of N'])
xlabel('Efficacy against mutant cells')
ylabel('Mutation rate')
zlabel('Probability')
set(gca,'FontSize',18)





%% panel 6
%% 6A
data = csvread('paper_data.csv'); % raw data that the prisoners dilemma paper stole that i am now stealing
paper_fit = csvread('paper_data_fit.csv'); % fit made by prisoners dilemma paper

x = data(:,1) - min(data(:,1)); % start x values at 0
y = data(:,2);

N = 1000; % given in prisoner dilemma paper
p1 = 1/N; % seems like plot starts with small number of cells

method = 'exact';
error_best = inf;
for Mperx = 8500:9500
    for fifj = 1:0.01:1.2
        Ni_n = Ni_nFunc(N,floor(x*Mperx),floor(p1*N),fifj,method); % calculate expected Ni using current fi/fj
        error = sqrt(mean((y - Ni_n).^2)); % RMSE between function and data
        if error < error_best % if current error is least so far, save current information
            Mperx_best = Mperx;
            fifj_best = fifj;
            error_best = error;
        end
    end
end

% plot

%results:
% Mperx_best = 8998;
% fifj_best = 1.1;

figure
scatter(x*Mperx_best,y, 10, 'filled')
hold on
plot(0:max(x*Mperx_best),Ni_nFunc(N,0:max(x*Mperx_best),floor(p1*N),fifj_best,method),'k','LineWidth',2)
xlabel('Time (cell divisions)')
ylabel('Proportion of population')
legend('Raw data', 'Model fit')
set(gca,'FontSize',18)


%% S5A, SI plot confirming RMSE
x = data(:,1) - min(data(:,1)); % start x values at 0
y = data(:,2);

N = 1000; % given in prisoner dilemma paper
p1 = 1/N; % seems like plot starts with small number of cells

Mperx_best = 8998; % best M per x unit, as found last week
fifj_vec = 0.8:0.001:2;
error_vec = zeros(size(fifj_vec));

method = 'exact';
error_best = inf;

for fifj_n = 1:length(fifj_vec)
    Ni_n = Ni_nFunc(N,floor(x*Mperx_best),floor(p1*N),fifj_vec(fifj_n),method); % calculate expected Ni using current fi/fj
    error_vec(fifj_n) = sqrt(mean((y - Ni_n).^2)); % RMSE between function and data
end

figure
plot(fifj_vec, error_vec)
xlabel('Fitness ratio f_c/f_h')
ylabel('RMSE')
set(gca,'FontSize',18)




%% 6B
fi_best = 1.035; % from fit in S5D, and assuming no mutant from result in 6A

% fitting treatment data
eff_i_vec = 0:0.001:1; % unknown eff_i_vec, range of values to select from
x_treat_n = floor(x_treat*24/time_per_n); % x data in cell divisions'

% for waning treatment:
time_per_n = 8/N; % 8 hours per cell for ehrlich cells to divide ('growth kinetics of the G2-phase of ehrlich ascites tumour cells, ...')
half_life = 5.8*24; % mean half-life of 5.8 days (https://www.drugs.com/medical-answers/long-herceptin-stay-body-3541312/)
drug_conc = @(n,ts) exp(log(0.5)*time_per_n/half_life*(n-ts)); % drug concentration as a ratio of its initial amount, with starting time ts
t_off = 0; % drug concentration that treatment is considered 'turned off'


fi = fi_best;
progress = 0;
error_best = inf;
Ni_avg = 0;
for eff_i = eff_i_vec
    Ni_vec = zeros(1,M+1); Ni_vec(1) = Ni_0;
    for n = 1:M
        if n < 7*24/time_per_n
            efi = fi*(1-eff_i*drug_conc(n,0));
        elseif n < 14*24/time_per_n
            efi = fi*(1-eff_i*drug_conc(n,7*24/time_per_n));
        elseif n < 21*24/time_per_n
            efi = fi*(1-eff_i*drug_conc(n,14*24/time_per_n));
        elseif drug_conc(n,21*24/time_per_n) > t_off
            efi = fi*(1-eff_i*drug_conc(n,21*24/time_per_n));
        else
            efi = fi;
        end   
        fifj = efi/fj;

        Ni_vec(n+1) = Ni_vec(n)*(1 - 1/N + fifj/(Ni_vec(n)*(fifj-1)+N));
    end
    %progress = progress + 1;
    %fprintf('%.3f%% complete\n', progress/length(eff_i_vec)*100)

    error = sqrt(mean((y_treat' - Ni_vec(x_treat_n+1)).^2)); % RMSE between function and data
    if error < error_best % if current error is least so far, save current information
        Ni_best_t = Ni_vec;
        eff_i_best = eff_i;
        error_best = error;
    end
end


% get eff_i_best = 0.055, or 5.5% effective (error_best = 46.9616)


% run n_total sims with this model fit
fi = fi_best; % as found before when fitting data
eff_i = eff_i_best;

n_total = 100;
Ni_mat = zeros(n_total,M+1); Ni_mat(:,1) = Ni_0;
Ni_t_mat = Ni_mat;

progress = 0;
for i = 1:n_total
    Ni = Ni_0; Nj = N - Ni;
    Ni_t = Ni_0; Nj_t = N - Ni_t;
    for n = 1:M
        Ni_n = Ni; Nj_n = Nj;
        Ni_t_n = Ni_t; Nj_t_n = Nj_t;
        if n < 7*24/time_per_n
            efi = fi*(1-eff_i*drug_conc(n,0));
        elseif n < 14*24/time_per_n
            efi = fi*(1-eff_i*drug_conc(n,7*24/time_per_n));
        elseif n < 21*24/time_per_n
            efi = fi*(1-eff_i*drug_conc(n,14*24/time_per_n));
        elseif drug_conc(n,21*24/time_per_n) > t_off
            efi = fi*(1-eff_i*drug_conc(n,21*24/time_per_n));
        else
            efi = fi;
        end
        
        % vehicle
        if rand < Ni_n*fi/(Ni_n*fi+Nj_n*fj)
            Ni = Ni + 1;
        else
            Nj = Nj + 1;
        end
        if rand < Ni_n/N
            Ni = Ni - 1;
        else
            Nj = Nj - 1;
        end
        Ni_mat(i,n+1) = Ni;

        % treatment
        if rand < Ni_t_n*efi/(Ni_t_n*efi+Nj_t_n*fj)
            Ni_t = Ni_t + 1;
        else
            Nj_t = Nj_t + 1;
        end
        if rand < Ni_t_n/N
            Ni_t = Ni_t - 1;
        else
            Nj_t = Nj_t - 1;
        end
        Ni_t_mat(i,n+1) = Ni_t;
    end
    %progress = progress + 1;
    %fprintf('%.2f%% complete\n',progress/n_total*100)
end


figure

errorbar(x_growth_actual*24/time_per_n,y_growth_actual,growth_error,'r.','LineWidth',1.5,'MarkerSize',12)
hold on
errorbar(x_treat_actual*24/time_per_n,y_treat_actual,treatment_error,'b.','LineWidth',1.5,'MarkerSize',12)
hold on

options.handle = gcf;
options.error = 'std';
options.alpha = 0.5;
options.line_width = 2;

options.color_line = cancerRGB./255;
options.color_area = cancerRGB./255;
plot_areaerrorbar(Ni_mat,options)

hold on
options.color_line = treatedcancerRGB./255;
options.color_area = treatedcancerRGB./255;
plot_areaerrorbar(Ni_t_mat,options)

hold on
plot((0:M),Ni_best,'--','Color',cancer)
plot((0:M),Ni_best_t,'--','Color',treatedcancer)

xline(7/time_per_n*24+1,'k--'); xline(14/time_per_n*24+1,'k--'); xline(21/time_per_n*24+1,'k--');

xlabel('Time (days)')
ylabel('Cancer volume (mm^3)')
%legend('Vehicle data','Treatment data','','Vehicle fit','','Treatment fit','Expected')
xticks([0 25 50 75 100 125]*24/time_per_n);
set(gca,'XTick', [0 25 50 75 100 125]*24/time_per_n)
set(gca,'XTickLabel', [0 25 50 75 100 125])
set(gca,'FontSize',18)
box on

%% S5BC, confirming treatment model fits
fi_vec = 0.8:0.0001:1.5;
error_v = zeros(size(fi_vec)); % errors for vehicle data
eff_i_vec = 0:0.0001:0.2;
error_t = zeros(size(eff_i_vec)); % treatment

progress = 0;
for fi_n = 1:length(fi_vec)
    fifj = fi_vec(fi_n)/fj;
    Ni_vec = zeros(1,M+1); Ni_vec(1) = Ni_0;
    for n = 1:M
        Ni_vec(n+1) = Ni_vec(n)*(1 - 1/N + fifj/(Ni_vec(n)*(fifj-1)+N));
    end
    progress = progress + 1;
    fprintf('%.3f%% complete\n', progress/length(fi_vec)*100)

    error_v(fi_n) = sqrt(mean((y_growth' - Ni_vec(x_growth_n+1)).^2)); % RMSE between function and data
end

progress = 0;
for eff_i_n = 1:length(eff_i_vec)
    eff_i = eff_i_vec(eff_i_n);
    Ni_vec = zeros(1,M+1); Ni_vec(1) = Ni_0;
    for n = 1:M
        if n < 7*24/time_per_n
            efi = fi*(1-eff_i*drug_conc(n,0));
        elseif n < 14*24/time_per_n
            efi = fi*(1-eff_i*drug_conc(n,7*24/time_per_n));
        elseif n < 21*24/time_per_n
            efi = fi*(1-eff_i*drug_conc(n,14*24/time_per_n));
        elseif drug_conc(n,21*24/time_per_n) > t_off
            efi = fi*(1-eff_i*drug_conc(n,21*24/time_per_n));
        else
            efi = fi;
        end
        fifj = efi/fj;

        Ni_vec(n+1) = Ni_vec(n)*(1 - 1/N + fifj/(Ni_vec(n)*(fifj-1)+N));
    end
    progress = progress + 1;
    fprintf('%.3f%% complete\n', progress/length(eff_i_vec)*100)

    error_t(eff_i_n) = sqrt(mean((y_treat' - Ni_vec(x_treat_n+1)).^2)); % RMSE between function and data
end


figure
plot(fi_vec, error_v)
xlabel('Fitness ratio f_c/f_h')
ylabel('RMSE')
set(gca,'FontSize',18)
xlim([0.8 1.5])

figure
plot(eff_i_vec, error_t)
xlabel('Drug efficacy eff_c')
ylabel('RMSE')
set(gca,'FontSize',18)



%% 6C
fi = fi_best;
eff_i = eff_i_best;

N = 5000; % assumed number of cells
fj = 1; % healthy cells
fifj = fi/fj; % ratio of cancerous to healthy cells
Ni_0 = y_growth(1); % starting cells

time_per_n = 8/N; % 8 hours per cell for ehrlich cells to divide ('growth kinetics of the G2-phase of ehrlich ascites tumour cells, ...')
t_start = 20*24/time_per_n; % arbitrary time to start treatment
t_length = 50*24/time_per_n; % arbitrary treatment length
st_weakness = 3; % weakness of sustained treatment compared to normal
st_length = st_weakness*t_length; % time taken for sustained treatment (low dose, applied enough times to give the same dosage as normal treatment)

M = floor(3*t_start+st_length); % arbitrary cell divisions for simulation length

time_per_n = 8/N; % 8 hours per cell for ehrlich cells to divide ('growth kinetics of the G2-phase of ehrlich ascites tumour cells, ...')
half_life = 5.8*24; % mean half-life of 5.8 days (https://www.drugs.com/medical-answers/long-herceptin-stay-body-3541312/)
drug_conc = @(n,ts) exp(log(0.5)*time_per_n/half_life*(n-ts)); % drug concentration as a ratio of its initial amount, with starting time ts


Ni_vec = zeros(1, M+1); Ni_vec(1) = Ni_0; % simulation without treatment
Ni_t_vec = Ni_vec; % simulation with normal treatment
Ni_st_vec = Ni_vec; % simulation with sustained treatment 

Ni = Ni_0; Nj = N-Ni_0;
Ni_t = Ni; Nj_t = Nj;
Ni_st = Ni; Nj_st = Nj;

progress = 0;
for n = 1:M
    Ni_n = Ni; Nj_n = Nj; % store current N1 and N2 for probability calculations
    
    % select a cell type to proliferate based on fitness
    r = rand; % generate random number to be used for probabilities
    if r < Ni_n*fifj/(Ni_n*fifj+Nj_n) % if random value is less than the probability of type i proliferating
        Ni = Ni + 1; % increase the number of cells of type i
    else
        Nj = Nj + 1; % increase the number of cells of type j
    end
    
    % select a cell type to die
    r = rand; % generate a new random number for new probabilities
    if r < Ni_n/N % if random value is less than the probability of type i dying
        Ni = Ni - 1; % decrease the number of cells of type i
    else
        Nj = Nj - 1; % decrease the number of cells of type j
    end
    
    Ni_vec(n+1) = Ni; % input new value to vector

    if n < t_start
        Ni_t = Ni; Nj_t = Nj; % keep same data for treatment sims if treatment hasn't started
        Ni_st = Ni; Nj_st = Nj;
        Ni_t_vec(n+1) = Ni;
        Ni_st_vec(n+1) = Ni;
    elseif n >= t_start % if treatment has started, set effective fitness
        if n < t_start+t_length
            efi_t = fi*(1-eff_i);
            efi_st = fi*(1-eff_i/st_weakness);
        elseif n < t_start+st_length % normal treatment starts waning
            efi_t = fi*(1-eff_i*drug_conc(n,t_start+t_length)); 
            efi_st = fi*(1-eff_i/st_weakness);
        else % if treatment has ended, treatment wanes
            efi_t = fi*(1-eff_i*drug_conc(n,t_start+t_length));
            efi_st = fi*(1-eff_i/st_weakness*drug_conc(n,t_start+st_length));
        end
        efifj_t = efi_t/fj;
        efifj_st = efi_st/fj;

        Ni_t_n = Ni_t; Nj_t_n = Nj_t;
        Ni_st_n = Ni_st; Nj_st_n = Nj_st;

        if rand < Ni_t_n*efifj_t/(Ni_t_n*efifj_t+Nj_t_n) % if random value is less than the probability of type i proliferating
            Ni_t = Ni_t + 1; % increase the number of cells of type i
        else
            Nj_t = Nj_t + 1; % increase the number of cells of type j
        end

        if rand < Ni_t_n/N % if random value is less than the probability of type i dying
            Ni_t = Ni_t - 1; % decrease the number of cells of type i
        else
            Nj_t = Nj_t - 1; % decrease the number of cells of type j
        end

        Ni_t_vec(n+1) = Ni_t;

        if rand < Ni_st_n*efifj_st/(Ni_st_n*efifj_st+Nj_st_n) % if random value is less than the probability of type i proliferating
            Ni_st = Ni_st + 1; % increase the number of cells of type i
        else
            Nj_st = Nj_st + 1; % increase the number of cells of type j
        end

        if rand < Ni_st_n/N % if random value is less than the probability of type i dying
            Ni_st = Ni_st - 1; % decrease the number of cells of type i
        else
            Nj_st = Nj_st - 1; % decrease the number of cells of type j
        end

        Ni_st_vec(n+1) = Ni_st;
    end

    %progress = progress + 1;
    %fprintf('%.3f%% complete\n', progress/M*100)
end


% recursive
Ni_vec2 = zeros(1, M+1); Ni_vec2(1) = Ni_0; % simulation without treatment
Ni_t_vec2 = Ni_vec2; % simulation with normal treatment
Ni_st_vec2 = Ni_vec2; % simulation with sustained treatment 
for n = 1:M
    Ni_vec2(n+1) = Ni_vec2(n)*(1 - 1/N + fifj/(Ni_vec2(n)*(fifj-1)+N));
    if n < t_start
        Ni_t_vec2(n+1) = Ni_vec2(n+1);
        Ni_st_vec2(n+1) = Ni_vec2(n+1);
    elseif n >= t_start % if treatment has started, set effective fitness
        if n < t_start+t_length
            efi_t = fi*(1-eff_i);
            efi_st = fi*(1-eff_i/st_weakness);
        elseif n < t_start+st_length % normal treatment starts waning
            efi_t = fi*(1-eff_i*drug_conc(n,t_start+t_length)); 
            efi_st = fi*(1-eff_i/st_weakness);
        else % if treatment has ended, treatment wanes
            efi_t = fi*(1-eff_i*drug_conc(n,t_start+t_length));
            efi_st = fi*(1-eff_i/st_weakness*drug_conc(n,t_start+st_length));
        end
        efifj_t = efi_t/fj;
        efifj_st = efi_st/fj;

        Ni_t_vec2(n+1) = Ni_t_vec2(n)*(1 - 1/N + efifj_t/(Ni_t_vec2(n)*(efifj_t-1)+N));
        Ni_st_vec2(n+1) = Ni_st_vec2(n)*(1 - 1/N + efifj_st/(Ni_st_vec2(n)*(efifj_st-1)+N));
    end
end


% plot
figure
plot((0:M)*time_per_n/24,Ni_vec/N,'Color',cancer)
hold on
plot((0:M)*time_per_n/24,Ni_t_vec/N,'Color',treatedcancer)
plot((0:M)*time_per_n/24,Ni_st_vec/N,'Color',treatedmutant1)
plot((0:M)*time_per_n/24,Ni_vec2/N,'--','Color',cancer)
plot((0:M)*time_per_n/24,Ni_t_vec2/N,'--','Color',treatedcancer)
plot((0:M)*time_per_n/24,Ni_st_vec2/N,'--','Color',treatedmutant1)
xline((t_start+1)*time_per_n/24,'k--')
xline((t_start+t_length+1)*time_per_n/24,'k--')
xline((t_start+st_length+1)*time_per_n/24,'k--')
xlabel('Time (days)')
ylabel('Proportion of population')
legend('No treatment','Normal treatment','Sustained treatment', 'Expected')%,'Expected without treatment','Expected normal treatment','Expected sustained treatment')
set(gca,'FontSize',18)



%% other SI plots
%% S1ABD (3AB with mean / standard deviation)
N = 1000;          % total number of cells at each timestep
M = 100000;          % total number of timesteps
pc = 0.1;         % starting proportion for cancer cells
fh = 1;           % fitness of healthy cells
fc = fh;          % fitness of cancer cells
n_total = 1000;   % number of simulations

sims_mat = zeros(n_total, M+1); % matrix to store simulation information

progress = 0;
for i = 1:n_total % run n_total simulations
    Nc = floor(pc*N); % starting number of cancer cells
    Nh = N - Nc;      % starting number of healthy cells
    sims_mat(i,1) = Nc/N; % initialise the current simulation

    for n = 1:M
        Nc_n = Nc; Nh_n = Nh; % store current Nc and Nh for probability calculations

        % select a cell type to proliferate based on fitness
        r = rand; % generate random number to be used for probabilities
        if r < Nc_n*fc/(Nc_n*fc+Nh_n*fh) % if random value is less than the probability of cancer proliferating
            Nc = Nc + 1; % increase the number of cancer cells
        else
            Nh = Nh + 1; % increase the number of healthy cells
        end

        % select a cell type to die
        r = rand; % generate a new random number for new probabilities
        if r < Nc_n/N % if random value is less than the probability of cancer dying
            Nc = Nc - 1; % decrease the number of cancer cells
        else
            Nh = Nh - 1; % decrease the number of healthy cells
        end
        sims_mat(i, n+1) = Nc/N; % update matrix with simulation value for the current fc
    end
    progress = progress + 1;
    fprintf('%.3f%% complete\n', progress/n_total*100)
end

figure
options.handle = gcf;
options.error = 'std';
options.alpha = 0.5;
options.line_width = 2;

options.color_line = cancerRGB./255;
options.color_area = cancerRGB./255;
plot_areaerrorbar(sims_mat,options)

hold on
options.color_line = healthyRGB./255;
options.color_area = healthyRGB./255;
plot_areaerrorbar(1-sims_mat)

xlabel('Time (cell divisions)')
ylabel('Proportion of cells')
xlim([0,M])
ylim([0,1])
set(gca,'FontSize',18)


%% S1E
f1_vec = 0:0.1:10; % fitness values for cell type 1
avg_vec = zeros(size(f1_vec)); % average values of N1 for each f1

N = 100;          % total number of cells at each timestep
M_vec = [10,50,100,150,200,300,500,1000];     % total number of timesteps
p1 = 0.5;         % starting proportion for cells of type 1
f2 = 5;           % fitness of type 2 cells
n_total = 1000;   % number of simulations

figure
hold on

for M_n = 1:length(M_vec)
    for f1_n = 1:length(f1_vec)
        avg = 0; % initialise average result
        for n = 1:n_total % take an average over n_total simulations
            N1 = floor(p1*N); % starting number of type 1 cells
            N2 = N - N1;      % starting number of type 2 cells

            for i = 1:M_vec(M_n)
                N1_i = N1; N2_i = N2; % store current N1 and N2 for probability calculations

                % select a cell type to proliferate based on fitness
                r = rand; % generate random number to be used for probabilities
                if r < N1_i*f1_vec(f1_n)/(N1_i*f1_vec(f1_n)+N2_i*f2) % if random value is less than the probability of type 1 proliferating
                    N1 = N1 + 1; % increase the number of cells of type 1
                else
                    N2 = N2 + 1; % increase the number of cells of type 2
                end

                % select a cell type to die
                r = rand; % generate a new random number for new probabilities
                if r < N1_i/N % if random value is less than the probability of type 1 dying
                    N1 = N1 - 1; % decrease the number of cells of type 1
                else
                    N2 = N2 - 1; % decrease the number of cells of type 2
                end
            end

            avg = avg + N1/n_total; % update average value
        end

        avg_vec(f1_n) = avg; % add average value of N1 for the current f1
    end
    
    plot(f1_vec, avg_vec/N)
end

hold off
legend('M = ' + string(M_vec))
%title(['Average percentage of cells of type 1 for different fitness values, where f_{2}=',num2str(f2)])
xlabel('Cancer fitness')
ylabel('Final cancer proportion')
ylim([0, 1])
set(gca,'FontSize',18)


%% S1E, plot of number of sims vs average Ni_M
N = 1000;
M = 10000;
fi = 1.1;
fj = 1;
Ni_0 = 0.1*N;

n_totvec = 1:2000;
Ni_M = zeros(1,length(n_totvec));

progress = 0;
for n_total = n_totvec
    for i = 1:n_total
        Ni = Ni_0; Nj = N - Ni;
        for n = 1:M
            Ni_n = Ni;
            if rand < Ni_n*fi/(Ni_n*fi+Nj*fj)
                Ni = Ni + 1;
            else
                Nj = Nj + 1;
            end
            if rand < Ni_n/N
                Ni = Ni - 1;
            else
                Nj = Nj - 1;
            end
        end
        Ni_M(n_total) = Ni_M(n_total) + Ni/N/n_total;

        progress = progress + 1;
        fprintf('%.2f%%\n',progress/sum(n_totvec)*100)
    end
end

figure
plot(n_totvec,Ni_M)
xlabel('Total simulations')
ylabel('Average final cancer proportion')
set(gca,'FontSize',18)




%% S2A
% 3-cell sim
N = 1000;
M = 20000;
fi = 1.1; % fit cells
fj = 1; % less fit cells
fk = 1.5; % very fit mutated cells
mk = 1/100; % percentage chance for k mutation given type i reproduces
Ni_0 = 0.1*N;
Nk_0 = 0;

t_start = 5000; % cell division to start treatment
t_length = 3000; % number of cell divisions to apply treatment over
eff_i = 1; % how effective treatment is on cell type i (0 to 1)
eff_k = 0.5; % how effective treatment is on cell type k (0 to 1)


% simualtion
Ni = Ni_0;
Nj = N - Ni_0 - Nk_0;
Nk = Nk_0;

Ni_vec = zeros(1,M+1); Ni_vec(1) = Ni;
Nk_vec = zeros(1,M+1); Nk_vec(1) = Nk;


for n = 1:M
    % store current values for probability calcs
    Ni_n = Ni; Nj_n = Nj; Nk_n = Nk;

    r = rand;
    if r < Ni_n*fi/(Ni_n*fi+Nj_n*fj+Nk_n*fk)
        if rand < mk
            Nk = Nk + 1; % mutate if random value is less than mk
        else
            Ni = Ni + 1;
        end
    elseif r >= Ni_n*fi/(Ni_n*fi+Nj_n*fj+Nk_n*fk) && r < (Ni_n*fi+Nj_n*fj)/(Ni_n*fi+Nj_n*fj+Nk_n*fk)
        Nj = Nj + 1;
    else
        Nk = Nk + 1;
    end

    r = rand;
    if r < Ni_n/N
        Ni = Ni - 1;
    elseif r >= Ni_n/N && r < (Ni_n+Nj_n)/N
        Nj = Nj - 1;
    else
        Nk = Nk - 1;
    end

    Ni_vec(n+1) = Ni; Nk_vec(n+1) = Nk;

    if n == t_start-1 % if right before the starting time for treatment
        Ni_t = Ni; Nj_t = Nj; Nk_t = Nk;
        Ni_t_vec = Ni_vec; Nk_t_vec = Nk_vec;
    elseif n >= t_start % run a new simulation where treatment was applied
        if n < t_start + t_length
            efi = fi*(1-eff_i); % effective fi, lower if treatment applied
            efk = fk*(1-eff_k);
        else
            efi = fi;
            efk = fk;
        end
        
        Ni_t_n = Ni_t; Nj_t_n = Nj_t; Nk_t_n = Nk_t;

        r = rand;
        if r < Ni_t_n*efi/(Ni_t_n*efi+Nj_t_n*fj+Nk_t_n*efk)
            if rand < mk
                Nk_t = Nk_t + 1; % mutate if random value is less than mk
            else
                Ni_t = Ni_t + 1;
            end
        elseif r >= Ni_t_n*efi/(Ni_t_n*efi+Nj_t_n*fj+Nk_t_n*efk) && r < (Ni_t_n*efi+Nj_t_n*fj)/(Ni_t_n*efi+Nj_t_n*fj+Nk_t_n*efk)
            Nj_t = Nj_t + 1;
        else
            Nk_t = Nk_t + 1;
        end
    
        r = rand;
        if r < Ni_t_n/N
            Ni_t = Ni_t - 1;
        elseif r >= Ni_t_n/N && r < (Ni_t_n+Nj_t_n)/N
            Nj_t = Nj_t - 1;
        else
            Nk_t = Nk_t - 1;
        end
    
        Ni_t_vec(n+1) = Ni_t; Nk_t_vec(n+1) = Nk_t;
    end
end


figure
plot(0:M,Ni_t_vec/N,"Color",treatedcancer)
hold on
plot(0:M,Nk_t_vec/N,"Color",treatedmutant1)
plot(0:M,Ni_vec/N,"Color",cancer)
plot(0:M,Nk_vec/N,"Color",mutant1)

hold on
xline(t_start+1, 'k--')
xline(t_start+t_length+1, 'k--')
xlabel('Time (cell divisions)')
ylabel('Proportion of population')
legend('Untreated original cancer','Untreated mutated cancer','Treated original cancer','Treated mutated cancer')
ylim([0 1])
set(gca,'FontSize',18)


%% S2C
N = 1000;
M = 30000;
fi = 1.1; % fit cells
fj = 1; % less fit cells
fk = 3; % very fit mutated cells
mk = 1/1000; % percentage chance for k mutation given type i reproduces
Ni_0 = 0.1*N;
Nk_0 = 0;

t_start = 5000; % cell division to start treatment
t_length = 5000; % number of cell divisions to apply treatment over
eff_i = 0.8; % how effective treatment is on cell type i (0 to 1)
eff_k = 0.5; % how effective treatment is on cell type k (0 to 1)


% example of 3-cell simualtion for control tumour growth
Ni = Ni_0;
Nj = N - Ni_0 - Nk_0;
Nk = Nk_0;
Ni_vec = zeros(1,M+1); Ni_vec(1) = Ni;
Nk_vec = zeros(1,M+1); Nk_vec(1) = Nk;

for n = 1:M
    % store current values for probability calcs
    Ni_n = Ni; Nj_n = Nj; Nk_n = Nk;

    r = rand;
    if r < Ni_n*fi/(Ni_n*fi+Nj_n*fj+Nk_n*fk)
        if rand < mk
            Nk = Nk + 1; % mutate if random value is less than mk
        else
            Ni = Ni + 1;
        end
    elseif r >= Ni_n*fi/(Ni_n*fi+Nj_n*fj+Nk_n*fk) && r < (Ni_n*fi+Nj_n*fj)/(Ni_n*fi+Nj_n*fj+Nk_n*fk)
        Nj = Nj + 1;
    else
        Nk = Nk + 1;
    end

    r = rand;
    if r < Ni_n/N
        Ni = Ni - 1;
    elseif r >= Ni_n/N && r < (Ni_n+Nj_n)/N
        Nj = Nj - 1;
    else
        Nk = Nk - 1;
    end

    Ni_vec(n+1) = Ni; Nk_vec(n+1) = Nk;
end



% now for the treatment
[Ni_t_vec3, Nk_t_vec3] = treat(N,fi,fj,Ni_vec,t_start,t_length,eff_i,Nk_vec,fk,mk,eff_k);
Ni_t_vec2 = treat(N,fi,fj,Ni_vec,t_start,t_length,eff_i,Nk_vec);

figure
plot(0:M,Ni_t_vec2/N,'k',0:M,(Ni_t_vec3+Nk_t_vec3)/N,'r')
hold on
xline(t_start+1, 'k--')
xline(t_start+t_length+1, 'k--')
xlabel('Time (cell divisions)')
ylabel('Proportion of population')
legend('Intra-clonally homogenous','Intra-clonally heterogeneous')
set(gca,'FontSize',18)



%% S2D (same data as 4B)
fi_vec = 0:0.01:5;
fk_vec = 0:0.01:5;
fix_i_mat_file = matfile('fix_i_mat.mat'); fix_i_mat = fix_i_mat_file.fix_i_mat;
fix_k_mat_file = matfile('fix_k_mat.mat'); fix_k_mat = fix_k_mat_file.fix_k_mat;

figure
ax1 = axes;
contour(fix_i_mat)
ax2 = axes;
contour(fix_k_mat)
linkaxes([ax1,ax2])
ax2.Visible = 'off';
% ax2.XTick = [];
% ax2.YTick = [];

colormap(ax1,'copper')
colormap(ax2,flipud(autumn))
set([ax1,ax2],'Position',[0.26 0.2 0.55 0.72])
colorbar(ax1,'Position',[0.08 0.2 0.0675 0.72])
colorbar(ax2,'Position',[0.85 0.2 0.0675 0.72])

% shading interp
clim(ax1,[0 1])
clim(ax2,[0 1])
xlabel(ax1, 'Mutant cancer fitness')
ylabel(ax1, 'Original cancer fitness')
xticks(ax1, [1 100 200 300 400 500]);
yticks(ax1, [100 200 300 400 500]);
set(ax1,'XTick', [1 100 200 300 400 500])
set(ax1,'XTickLabel', [0 1 2 3 4 5])
set(ax1,'YTick', [100 200 300 400 500])
set(ax1,'YTickLabel', [1 2 3 4 5])
set(ax1,'FontSize',18)
set(ax2,'FontSize',18)




%% S5D
% from 'Targeting HER2-Positive Breast Cancer with Trastuzumab-DM1, an Antibody–Cytotoxic Drug Conjugate'
growth_data = csvread('growth_data.csv'); % raw data for normal cancer growth
treatment_data = csvread('treatment_data.csv'); % raw data for treatment of herceptin
growth_errorvalues = csvread('growth_error.csv');
treatment_errorvalues = csvread('treatment_error.csv');
growth_actual = csvread('growth_data_actual.csv');
treatment_actual = csvread('treatment_data_actual.csv');
% known/assumed parameters
x_growth = growth_data(1:end-1,1); y_growth = growth_data(1:end-1,2);
x_treat = treatment_data(:,1); y_treat = treatment_data(:,2);
x_growth_actual = growth_actual(:,1); y_growth_actual = growth_actual(:,2);
x_treat_actual = treatment_actual(:,1); y_treat_actual = treatment_actual(:,2);
growth_error = abs(y_growth_actual-growth_errorvalues(:,2));
treatment_error = abs(y_treat_actual-treatment_errorvalues(:,2));

% consider single cancer clone
N = 5000; % 1mm^3 of cancer volume related to one cell, guessed based on max of data
fi_vec = 1:0.001:2; % unknown fi, range of values to select from
fj = 1; % healthy cells
Ni_0 = y_growth(1); % starting cells in data

time_per_n = 8/N; % 8 hours per cell for ehrlich cells to divide ('growth kinetics of the G2-phase of ehrlich ascites tumour cells, ...')
dosage = 15; % mg, dosage given in paper
effect_time = 9.35 + (66-9.35)*dosage/500; % time that dosage is effective for normal treatment, in days (https://www.drugs.com/medical-answers/long-herceptin-stay-body-3541312/)
t_start = 0; % start treatment immediately
%t_length = effect_time*24/time_per_n; % based off the effective time and time_per_n for ehrlich cells
t_length = 28*24/time_per_n; % ~28 days total for treatment

M = 130*24/time_per_n; % simulation length as long as data observation
x_growth_n = floor(x_growth*24/time_per_n); % x data in cell divisions

% fit fitness first using growth data
progress = 0;
error_best = inf;
Ni_avg = 0;
for fi = fi_vec
    fifj = fi/fj;
    Ni_vec = zeros(1,M+1); Ni_vec(1) = Ni_0;
    for n = 1:M
        Ni_vec(n+1) = Ni_vec(n)*(1 - 1/N + fifj/(Ni_vec(n)*(fifj-1)+N));
    end
    progress = progress + 1;
    fprintf('%.3f%% complete\n', progress/length(fi_vec)*100)

    error = sqrt(mean((y_growth' - Ni_vec(x_growth_n+1)).^2)); % RMSE between function and data
    if error < error_best % if current error is least so far, save current information
        Ni_best = Ni_vec;
        fi_best = fi;
        error_best = error;
    end
end


% fitting treatment data
eff_i_vec = 0:0.001:1; % unknown eff_i_vec, range of values to select from
x_treat_n = floor(x_treat*24/time_per_n); % x data in cell divisions

fi = fi_best;
progress = 0;
error_best = inf;
Ni_avg = 0;
for eff_i = eff_i_vec
    Ni_vec = zeros(1,M+1); Ni_vec(1) = Ni_0;
    for n = 1:M
        if n < t_length
            efi = fi*(1-eff_i);
        else
            efi = fi;
        end   
        fifj = efi/fj;

        Ni_vec(n+1) = Ni_vec(n)*(1 - 1/N + fifj/(Ni_vec(n)*(fifj-1)+N));
    end
    %progress = progress + 1;
    %fprintf('%.3f%% complete\n', progress/length(eff_i_vec)*100)

    error = sqrt(mean((y_treat' - Ni_vec(x_treat_n+1)).^2)); % RMSE between function and data
    if error < error_best % if current error is least so far, save current information
        Ni_best_t = Ni_vec;
        eff_i_best = eff_i;
        error_best = error;
    end
end

%
% results from fit  
% fi_best = 1.035;
% eff_i_best = 0.043;
% error_best = 39.9997

% run n_total sims with this model fit
fi = fi_best;
eff_i = eff_i_best;

n_total = 100;
Ni_mat = zeros(n_total,M+1); Ni_mat(:,1) = Ni_0;
Ni_t_mat = Ni_mat;

progress = 0;
for i = 1:n_total
    Ni = Ni_0; Nj = N - Ni;
    Ni_t = Ni_0; Nj_t = N - Ni_t;
    for n = 1:M
        Ni_n = Ni; Nj_n = Nj;
        Ni_t_n = Ni_t; Nj_t_n = Nj_t;
        if n < t_length
            efi = fi*(1-eff_i);
        else
            efi = fi;
        end   
        
        % vehicle
        if rand < Ni_n*fi/(Ni_n*fi+Nj_n*fj)
            Ni = Ni + 1;
        else
            Nj = Nj + 1;
        end
        if rand < Ni_n/N
            Ni = Ni - 1;
        else
            Nj = Nj - 1;
        end
        Ni_mat(i,n+1) = Ni;

        % treatment
        if rand < Ni_t_n*efi/(Ni_t_n*efi+Nj_t_n*fj)
            Ni_t = Ni_t + 1;
        else
            Nj_t = Nj_t + 1;
        end
        if rand < Ni_t_n/N
            Ni_t = Ni_t - 1;
        else
            Nj_t = Nj_t - 1;
        end
        Ni_t_mat(i,n+1) = Ni_t;
    end
    %progress = progress + 1;
    %fprintf('%.2f%% complete\n',progress/n_total*100)
end


figure

errorbar(x_growth_actual*24/time_per_n,y_growth_actual,growth_error,'r.','LineWidth',1.5,'MarkerSize',12)
hold on
errorbar(x_treat_actual*24/time_per_n,y_treat_actual,treatment_error,'b.','LineWidth',1.5,'MarkerSize',12)
hold on

options.handle = gcf;
options.error = 'std';
options.alpha = 0.5;
options.line_width = 2;

options.color_line = cancerRGB./255;
options.color_area = cancerRGB./255;
plot_areaerrorbar(Ni_mat,options)

hold on
options.color_line = treatedcancerRGB./255;
options.color_area = treatedcancerRGB./255;
plot_areaerrorbar(Ni_t_mat,options)

hold on
plot((0:M),Ni_best,'--','Color',cancer)
plot((0:M),Ni_best_t,'--','Color',treatedcancer)

xline((t_start+1),'k--')
xline((t_start+t_length+1),'k--')

xlabel('Time (days)')
ylabel('Cancer volume (mm^3)')
%legend('Vehicle data','Treatment data','','Vehicle fit','','Treatment fit','Expected')
xticks([0 25 50 75 100 125]*24/time_per_n);
set(gca,'XTick', [0 25 50 75 100 125]*24/time_per_n)
set(gca,'XTickLabel', [0 25 50 75 100 125])
set(gca,'FontSize',18)


%% S5E
% increased efficacy: analagous to higher dosage
eff_i = eff_i_best*2;
Ni_t_test1 = zeros(1,M+1); Ni_t_test1(1) = Ni_0;
for n = 1:M
    if n < t_length
        efi = fi*(1-eff_i);
    else
        efi = fi;
    end   
    fifj = efi/fj;

    Ni_t_test1(n+1) = Ni_t_test1(n)*(1 - 1/N + fifj/(Ni_t_test1(n)*(fifj-1)+N));
end

% sustained treatment: treatment length and efficacy proportional to st_weakness
eff_i = eff_i_best;
st_weakness = 2;
st_length = 28*st_weakness*24/time_per_n; % 28*st_weakness days
Ni_t_test2 = zeros(1,M+1); Ni_t_test2(1) = Ni_0;
for n = 1:M
    if n < st_length
        efi = fi*(1-eff_i/st_weakness);
    else
        efi = fi;
    end   
    fifj = efi/fj;

    Ni_t_test2(n+1) = Ni_t_test2(n)*(1 - 1/N + fifj/(Ni_t_test2(n)*(fifj-1)+N));
end


figure
errorbar(x_treat_actual,y_treat_actual,treatment_error,'b.','LineWidth',1.5,'MarkerSize',12)
hold on
plot((0:M)*time_per_n/24,Ni_best_t,'Color',treatedcancer)
plot((0:M)*time_per_n/24,Ni_t_test1,'Color',treatedmutant1)
plot((0:M)*time_per_n/24,Ni_t_test2,'Color',mutant2)
hold on
xline(t_length*time_per_n/24,'b--')
xline(st_length*time_per_n/24,'r--')
xlabel('Time (days)')
ylabel('Cancer volume (mm^3)')
legend('Treatment data','Treatment fit','Higher dosage','Sustained treatment')
set(gca,'FontSize',18)




%% FUNCTIONS
% treatment function for 2- or 3-cell models
function [Ni_t_vec, Nk_t_vec] = treat(N,fi,fj,Ni_vec,t_start,t_length,eff_i,Nk_vec,fk,mk,eff_k)
    % outputs vectors of treated cancers Ni_t_vec (normal cancer) and
    %  Nk_t_vec (mutated cancer) given the normal parameters, and the
    %  control cancer growths without treatment Ni_vec and Nk_vec, and the
    %  starting time for treatment t_start, the length of treatment
    %  t_length, the efficacy of treatment of type i and k cancers, eff_i
    %  and eff_k

    if nargin == 7
        % set default values for 3rd type when using a 2-cell model
        Nk_vec = zeros(1,length(Ni_vec)+1);
        fk = 1;
        mk = 0;
        eff_k = 0;
    elseif nargin == 8
        Ni_vec = Ni_vec + Nk_vec; % combine them if treating as a 2-cell model
        Nk_vec = zeros(1,length(Ni_vec)+1);
        fk = 1;
        mk = 0;
        eff_k = 0;
    elseif nargin ~= 11
        error('Please enter the correct arguments.')
    end
    
    M = length(Ni_vec) - 1;

    Ni = Ni_vec(t_start);
    Nj = N - Ni_vec(t_start) - Nk_vec(t_start);
    Nk = Nk_vec(t_start);
    
    Ni_t_vec = zeros(1,M+1); Ni_t_vec(1:t_start) = Ni_vec(1:t_start);
    Nk_t_vec = zeros(1,M+1); Nk_t_vec(1:t_start) = Nk_vec(1:t_start);
    
    
    for n = t_start:M
        % change effective fitnesses based on treatment
        if n < t_start + t_length
            efi = fi*(1-eff_i); % effective fi, lower if treatment applied
            efk = fk*(1-eff_k);
        else
            efi = fi;
            efk = fk;
        end

        % store current values for probability calcs
        Ni_n = Ni; Nj_n = Nj; Nk_n = Nk;
    
        r = rand;
        if r < Ni_n*efi/(Ni_n*efi+Nj_n*fj+Nk_n*efk)
            if rand < mk
                Nk = Nk + 1; % mutate if random value is less than mk
            else
                Ni = Ni + 1;
            end
        elseif r >= Ni_n*efi/(Ni_n*efi+Nj_n*fj+Nk_n*efk) && r < (Ni_n*efi+Nj_n*fj)/(Ni_n*efi+Nj_n*fj+Nk_n*efk)
            Nj = Nj + 1;
        else
            Nk = Nk + 1;
        end
    
        r = rand;
        if r < Ni_n/N
            Ni = Ni - 1;
        elseif r >= Ni_n/N && r < (Ni_n+Nj_n)/N
            Nj = Nj - 1;
        else
            Nk = Nk - 1;
        end
    
        Ni_t_vec(n+1) = Ni; Nk_t_vec(n+1) = Nk;
    end
end

function Ni_n = Ni_nFunc(N,M_vec,Ni_0,fifj,method)
    % Calculates the expected percentage of number of type i cells
    %  (Ni_next) after simulated cell divisions, 
    % Given the total number of cells (N), the number of cell divisions 
    %  (M), the initial number of type i cells (Ni_0), and the ratio 
    %  between type i fitness and type j fitness, fi/fj (fifj), and the
    %  method used to calculate/estimate Ni_n (method) out of 'exact' for
    %  the exact recursive solution, 'TS' for first order taylor series
    %  approximation, 'replace' for second replacement approximation,
    %  'infinite' for the infinite fi model (independent of fifj), or 'sim'
    %  for moran simulation
    
    Ni_M = zeros(1,max(M_vec)+1) + Ni_0; % use this as a template for iterations going up to max cell divisions
    Ni_n = zeros(size(M_vec));
    
    if strcmp(method, 'exact')
        for i = 2:max(M_vec)+1
            Ni_M(i) = Ni_M(i-1)*(1 - 1/N + fifj/(Ni_M(i-1)*(fifj-1)+N));
        end
        for M_n = 1:length(M_vec)
            Ni_n(M_n) = Ni_M(M_vec(M_n)+1);
        end
        
    elseif strcmp(method, 'TS')
        if fifj ~= 1
            A = 1 - 1/N + N*fifj/(N/2*(fifj-1)+N)^2;
            alpha = Ni_0/(1-A)*(1/N - fifj/(Ni_0*(fifj-1)+N));
            beta = Ni_0 - alpha;
            Ni_n = alpha*A.^M_vec + beta;
        end
        
    elseif strcmp(method, 'replace')
        for M_n = 1:length(M_vec)
            %Ni_n(M_n) = Ni_0*prod(1-1/N+fifj./((Ni_0+(fifj-1)/(2*fifj+2)*(0:M_vec(M_n)-1))*(fifj-1)+N));
            Ni_n(M_n) = N - (N-Ni_0)*prod(1-1/N+1/fifj./((N-Ni_0+(1/fifj-1)/(2/fifj+2)*(0:M_vec(M_n)-1))*(1/fifj-1)+N));
        end
    elseif strcmp(method, 'infinite')
        Ni_n = (Ni_0-N)*((N-1)/N).^M_vec + N;

    elseif strcmp(method, 'sim')
        Ni_M = zeros(1, max(M_vec)+1);
        ntotal = 1000;
        for n = 1:ntotal
            N1_vec = zeros(1, max(M_vec)+1); N1_vec(1) = Ni_0;
            N1 = Ni_0; N2 = N-Ni_0;
            for i = 1:max(M_vec)
                N1_i = N1; N2_i = N2; % store current N1 and N2 for probability calculations
                
                % select a cell type to proliferate based on fitness
                r = rand; % generate random number to be used for probabilities
                if r < N1_i*fifj/(N1_i*fifj+N2_i) % if random value is less than the probability of type 1 proliferating
                    N1 = N1 + 1; % increase the number of cells of type 1
                else
                    N2 = N2 + 1; % increase the number of cells of type 2
                end
                
                % select a cell type to die
                r = rand; % generate a new random number for new probabilities
                if r < N1_i/N % if random value is less than the probability of type 1 dying
                    N1 = N1 - 1; % decrease the number of cells of type 1
                else
                    N2 = N2 - 1; % decrease the number of cells of type 2
                end
                
                N1_vec(i+1) = N1; % input new value to vector
            end
            Ni_M = Ni_M + N1_vec/ntotal;
        end

        for M_n = 1:length(M_vec)
            Ni_n(M_n) = Ni_M(M_vec(M_n)+1);
        end
        
    else % if no proper method inputted, throw error
        error('Please select ''exact'', ''TS'', ''replace'', or ''infinite'' for the method.')
    end
    
    Ni_n = Ni_n/N; % complete as a percentage of population
end

% FUNCTION ^^
