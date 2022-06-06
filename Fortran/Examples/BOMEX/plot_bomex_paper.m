clear;
close all;
clc;

fileID = fopen('output', 'w');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Target solution wit nz = 601
% CFL= 0.1 with physics
%
DW600  = load('bomex_data/bomexweno600_1.dat');

lw = 8;
ms = 10;
fs = 60;
fs2 = 20
figure;clf
start_idx=240;
end_idx = 240;

for i=start_idx:end_idx 
    fprintf('i = %d \n', i);
    for k=1:5
      if(k<5)
        subplot(2,2, k)
      else
        figure
      end
      %axis square
      %grid on 
      nz=600;
      s=(i-1)*nz+1;e=s+nz-1;
      if(k==5 )
        plot(DW600(s:e,1+k)*1e+3, DW600(s:e, 1)*1e-3, '-k', 'LineWidth', lw, 'MarkerSize', ms)
        ylim([0.5 2])
      else
        plot(DW600(s:e,1+k), DW600(s:e, 1)*1e-3, '-k', 'LineWidth', lw, 'MarkerSize', ms)
      end
     
      ylabel('z (km)')
      if(k==1)
        xlabel('K')
        title('Potential Temperature')
      elseif(k==2)
        title('Water Vapor Mixing Ratio')
        xlabel('kg/kg')
      elseif(k==3)
        title('U Wind')
        xlabel('m/s')
      elseif(k==4)
        title('V Wind')
        xlabel('m/s')
      elseif(k==5)
        %title('Cloud Mixing Ratio')
        xlabel('g/kg')
      end
      set(gca, 'FontSize', fs)

      hold off;

    end
    pause
end

%
DSTD  = load('bomex_data/bomexweno600_1Standard.dat');

lw = 8;
ms = 10;
fs = 60;
figure;clf

for i=start_idx:end_idx 
    fprintf('i = %d \n', i)
    for k=1:5
      if(k<5)
        subplot(2,2, k)
      else
        figure
      end
      %axis square
      %grid on 
      nz=600;
      s=(i-1)*nz+1;e=s+nz-1;
      if(k==5 )
        plot(DW600(s:e,1+k)*1e+3, DW600(s:e, 1)*1e-3, '-k', 'LineWidth', lw, 'MarkerSize', ms)
        hold on
        plot(DSTD(s:e,1+k)*1e+3,  DSTD(s:e, 1)*1e-3,  '-b',  'LineWidth', lw, 'MarkerSize', ms)
        legend('Target', 'Standard', 'Interpreter', 'latex')
        ylim([0.7 2])
        hold off
        %title('Cloud Mixing Ratio')
        xlabel('g/kg')
        set(gca, 'FontSize', fs)
        axes('Position',[.7 .7 .2 .2])
        box on
        plot(DW600(s:e,1+k)*1e+3, DW600(s:e, 1)*1e-3, '-k', 'LineWidth', lw, 'MarkerSize', ms)
        hold on
        plot(DSTD(s:e,1+k)*1e+3,  DSTD(s:e, 1)*1e-3,  '-b',  'LineWidth', lw, 'MarkerSize', ms)
        ylim([1.5 1.55])
        xlim([-0.005, 0.005])
        hold off
      else
        plot(DW600(s:e,1+k), DW600(s:e, 1)*1e-3, '-k', 'LineWidth', lw, 'MarkerSize', ms)
        hold on
        plot(DSTD(s:e,1+k),  DSTD(s:e, 1)*1e-3,  '-b',  'LineWidth', lw, 'MarkerSize', ms)
        legend('Target', 'Standard', 'Interpreter', 'latex')
      end
     
      if(k==1)
        xlabel('K')
        title('Potential Temperature')
      elseif(k==2)
        title('Water Vapor Mixing Ratio')
        xlabel('kg/kg')
      elseif(k==3)
        title('U Wind')
        xlabel('m/s')
      elseif(k==4)
        title('V Wind')
        xlabel('m/s')
      elseif(k==5)
        %title('Cloud Mixing Ratio')
        %xlabel('g/kg')
      end
      if(k~=5)
      set(gca, 'FontSize', fs)
      ylabel('z (km)')
      else
      set(gca, 'FontSize', 20)
      end

      hold off;

    end
    %pause
end

DCLIP  = load('bomex_data/bomexweno600_1Clipping.dat');
figure;clf
for i=start_idx:end_idx 
    fprintf('i = %d \n', i)
    for k=1:5
      if(k<5)
        subplot(2,2, k)
      else
        figure
      end
      %axis square
      %grid on 
      nz=600;
      s=(i-1)*nz+1;e=s+nz-1;
      if(k==5 )
        plot(DW600(s:e,1+k)*1e+3, DW600(s:e, 1)*1e-3, '-k', 'LineWidth', lw, 'MarkerSize', ms)
        hold on
        plot(DCLIP(s:e,1+k)*1e+3,  DCLIP(s:e, 1)*1e-3,  '-b',  'LineWidth', lw, 'MarkerSize', ms)
        legend('Target', 'Clipping', 'Interpreter', 'latex')
        ylim([0.7 2])
        hold off
        %title('Cloud Mixing Ratio')
        xlabel('g/kg')
        set(gca, 'FontSize', fs)
        axes('Position',[.7 .7 .2 .2])
        box on
        plot(DW600(s:e,1+k)*1e+3, DW600(s:e, 1)*1e-3, '-k', 'LineWidth', lw, 'MarkerSize', ms)
        hold on
        plot(DCLIP(s:e,1+k)*1e+3,  DCLIP(s:e, 1)*1e-3,  '-b',  'LineWidth', lw, 'MarkerSize', ms)
        ylim([1.57 1.6])
        xlim([-0.005, 0.005])
        hold off
      else
        plot(DW600(s:e,1+k), DW600(s:e, 1)*1e-3, '-k', 'LineWidth', lw, 'MarkerSize', ms)
        hold on
        plot(DCLIP(s:e,1+k),  DCLIP(s:e, 1)*1e-3,  '-b',  'LineWidth', lw, 'MarkerSize', ms)
        legend('Target', 'Clipping', 'Interpreter', 'latex')
      end
     
      if(k==1)
        xlabel('K')
        title('Potential Temperature')
      elseif(k==2)
        title('Water Vapor Mixing Ratio')
        xlabel('kg/kg')
      elseif(k==3)
        title('U Wind')
        xlabel('m/s')
      elseif(k==4)
        title('V Wind')
        xlabel('m/s')
      elseif(k==5)
        %title('Cloud Mixing Ratio')
        %xlabel('g/kg')
      end
      if(k~=5)
      set(gca, 'FontSize', fs)
      ylabel('z (km)')
      else
      set(gca, 'FontSize', 20)
      end

      hold off;

    end
    %pause
end

DPCHIP  = load('bomex_data/bomexweno600_1PCHIP.dat');
figure;clf
for i=start_idx:end_idx 
    fprintf('i = %d \n', i)
    for k=1:5
      if(k<5)
        subplot(2,2, k)
      else
        figure
      end
      %axis square
      %grid on 
      nz=600;
      s=(i-1)*nz+1;e=s+nz-1;
      if(k==5 )
        plot(DW600(s:e,1+k)*1e+3, DW600(s:e, 1)*1e-3, '-k', 'LineWidth', lw, 'MarkerSize', ms)
        hold on
        plot(DPCHIP(s:e,1+k)*1e+3,  DPCHIP(s:e, 1)*1e-3,  '-b',  'LineWidth', lw, 'MarkerSize', ms)
        legend('Target', 'PCHIP', 'Interpreter', 'latex')
        ylim([0.7 2])
        hold off
        %title('Cloud Mixing Ratio')
        xlabel('g/kg')
        set(gca, 'FontSize', fs)
        axes('Position',[.7 .7 .2 .2])
        box on
        plot(DW600(s:e,1+k)*1e+3, DW600(s:e, 1)*1e-3, '-k', 'LineWidth', lw, 'MarkerSize', ms)
        hold on
        plot(DPCHIP(s:e,1+k)*1e+3,  DPCHIP(s:e, 1)*1e-3,  '-b',  'LineWidth', lw, 'MarkerSize', ms)
        ylim([1.35 1.45])
        xlim([-0.005, 0.005])
        hold off
      else
        plot(DW600(s:e,1+k), DW600(s:e, 1)*1e-3, '-k', 'LineWidth', lw, 'MarkerSize', ms)
        hold on
        plot(DPCHIP(s:e,1+k),  DPCHIP(s:e, 1)*1e-3,  '-b',  'LineWidth', lw, 'MarkerSize', ms)
        legend('Target', 'PCHIP', 'Interpreter', 'latex')
      end
     
      if(k==1)
        xlabel('K')
        title('Potential Temperature')
      elseif(k==2)
        title('Water Vapor Mixing Ratio')
        xlabel('kg/kg')
      elseif(k==3)
        title('U Wind')
        xlabel('m/s')
      elseif(k==4)
        title('V Wind')
        xlabel('m/s')
      elseif(k==5)
        %title('Cloud Mixing Ratio')
        %xlabel('g/kg')
      end
      if(k~=5)
      set(gca, 'FontSize', fs)
      ylabel('z (km)')
      else
      set(gca, 'FontSize', 20)
      end

      hold off;

    end
    %pause
end


DDBI5S1  = load('bomex_data/bomexweno600_1DBI5st=1eps0=1.0e-5eps1=1.0e-5.dat');
DDBI7S1  = load('bomex_data/bomexweno600_1DBI7st=1eps0=1.0e-5eps1=1.0e-5.dat');
DDBI5S2  = load('bomex_data/bomexweno600_1DBI5st=2eps0=1.0e-5eps1=1.0e-5.dat');
DDBI7S2  = load('bomex_data/bomexweno600_1DBI7st=2eps0=1.0e-5eps1=1.0e-5.dat');
DDBI5S3  = load('bomex_data/bomexweno600_1DBI5st=3eps0=1.0e-5eps1=1.0e-5.dat');
DDBI7S3  = load('bomex_data/bomexweno600_1DBI7st=3eps0=1.0e-5eps1=1.0e-5.dat');
for i=start_idx:end_idx 
      fprintf('i = %d \n', i)
      nz=600;
      s=(i-1)*nz+1;e=s+nz-1;

      figure;clf
      plot(DW600(s:e,6)*1e+3, DW600(s:e, 1)*1e-3, '-k', 'LineWidth', lw, 'MarkerSize', ms)
      hold on
      plot(DDBI5S1(s:e,6)*1e+3,  DDBI5S1(s:e, 1)*1e-3,  '-b',  'LineWidth', lw, 'MarkerSize', ms)
      plot(DDBI7S1(s:e,6)*1e+3,  DDBI7S1(s:e, 1)*1e-3,  '-r',  'LineWidth', lw, 'MarkerSize', ms)
      legend('Target', 'DBI $$\mathcal{P}_{5}$$', 'DBI $$\mathcal{P}_{7}$$', 'Interpreter', 'latex')
      ylim([0.7 2])
      hold off
      %title('Cloud Mixing Ratio')
      xlabel('g/kg')
      ylabel('z (km)')
      set(gca, 'FontSize', fs)
      axes('Position',[.7 .7 .2 .2])
      box on
      plot(DW600(s:e,6)*1e+3, DW600(s:e, 1)*1e-3, '-k', 'LineWidth', lw, 'MarkerSize', ms)
      hold on
      plot(DDBI5S1(s:e,6)*1e+3,  DDBI5S1(s:e, 1)*1e-3,  '-b',  'LineWidth', lw, 'MarkerSize', ms)
      plot(DDBI7S1(s:e,6)*1e+3,  DDBI7S1(s:e, 1)*1e-3,  '-r',  'LineWidth', lw, 'MarkerSize', ms)
      ylim([1.35 1.45])
      xlim([-0.005, 0.005])
      hold off
      set(gca, 'FontSize', fs2)

      figure;clf
      plot(DW600(s:e,6)*1e+3, DW600(s:e, 1)*1e-3, '-k', 'LineWidth', lw, 'MarkerSize', ms)
      hold on
      plot(DDBI5S2(s:e,6)*1e+3,  DDBI5S2(s:e, 1)*1e-3,  '-b',  'LineWidth', lw, 'MarkerSize', ms)
      plot(DDBI7S2(s:e,6)*1e+3,  DDBI7S2(s:e, 1)*1e-3,  '-r',  'LineWidth', lw, 'MarkerSize', ms)
      legend('Target', 'DBI $$\mathcal{P}_{5}$$', 'DBI $$\mathcal{P}_{7}$$', 'Interpreter', 'latex')
      ylim([0.7 2])
      hold off
      %title('Cloud Mixing Ratio')
      xlabel('g/kg')
      ylabel('z (km)')
      set(gca, 'FontSize', fs)
      axes('Position',[.7 .7 .2 .2])
      box on
      plot(DW600(s:e,6)*1e+3, DW600(s:e, 1)*1e-3, '-k', 'LineWidth', lw, 'MarkerSize', ms)
      hold on
      plot(DDBI5S2(s:e,6)*1e+3,  DDBI5S2(s:e, 1)*1e-3,  '-b',  'LineWidth', lw, 'MarkerSize', ms)
      plot(DDBI7S2(s:e,6)*1e+3,  DDBI7S2(s:e, 1)*1e-3,  '-r',  'LineWidth', lw, 'MarkerSize', ms)
      ylim([1.4 1.50])
      xlim([-0.005, 0.005])
      hold off
      set(gca, 'FontSize', fs2)

      figure;clf
      plot(DW600(s:e,6)*1e+3, DW600(s:e, 1)*1e-3, '-k', 'LineWidth', lw, 'MarkerSize', ms)
      hold on
      plot(DDBI5S3(s:e,6)*1e+3,  DDBI5S3(s:e, 1)*1e-3,  '-b',  'LineWidth', lw, 'MarkerSize', ms)
      plot(DDBI7S3(s:e,6)*1e+3,  DDBI7S3(s:e, 1)*1e-3,  '-r',  'LineWidth', lw, 'MarkerSize', ms)
      legend('Target', 'DBI $$\mathcal{P}_{5}$$', 'DBI $$\mathcal{P}_{7}$$', 'Interpreter', 'latex')
      ylim([0.7 2])
      hold off
      %title('Cloud Mixing Ratio')
      xlabel('g/kg')
      ylabel('z (km)')
      set(gca, 'FontSize', fs)
      axes('Position',[.7 .7 .2 .2])
      box on
      plot(DW600(s:e,6)*1e+3, DW600(s:e, 1)*1e-3, '-k', 'LineWidth', lw, 'MarkerSize', ms)
      hold on
      plot(DDBI5S3(s:e,6)*1e+3,  DDBI5S3(s:e, 1)*1e-3,  '-b',  'LineWidth', lw, 'MarkerSize', ms)
      plot(DDBI7S3(s:e,6)*1e+3,  DDBI7S3(s:e, 1)*1e-3,  '-r',  'LineWidth', lw, 'MarkerSize', ms)
      ylim([1.4 1.45])
      xlim([-0.005, 0.005])
      hold off
      set(gca, 'FontSize', fs2)



end

DPPI5S1  = load('bomex_data/bomexweno600_1PPI5st=1eps0=1.0e-5eps1=1.0e-5.dat');
DPPI7S1  = load('bomex_data/bomexweno600_1PPI7st=1eps0=1.0e-5eps1=1.0e-5.dat');
DPPI5S2  = load('bomex_data/bomexweno600_1PPI5st=2eps0=1.0e-5eps1=1.0e-5.dat');
DPPI7S2  = load('bomex_data/bomexweno600_1PPI7st=2eps0=1.0e-5eps1=1.0e-5.dat');
DPPI5S3  = load('bomex_data/bomexweno600_1PPI5st=3eps0=1.0e-5eps1=1.0e-5.dat');
DPPI7S3  = load('bomex_data/bomexweno600_1PPI7st=3eps0=1.0e-5eps1=1.0e-5.dat');
figure;clf
for i=start_idx:end_idx 
      fprintf('i = %d \n', i)
      nz=600;
      s=(i-1)*nz+1;e=s+nz-1;

      figure;clf
      plot(DW600(s:e,6)*1e+3, DW600(s:e, 1)*1e-3, '-k', 'LineWidth', lw, 'MarkerSize', ms)
      hold on
      plot(DPPI5S1(s:e,6)*1e+3,  DPPI5S1(s:e, 1)*1e-3,  '-b',  'LineWidth', lw, 'MarkerSize', ms)
      plot(DPPI7S1(s:e,6)*1e+3,  DPPI7S1(s:e, 1)*1e-3,  '-r',  'LineWidth', lw, 'MarkerSize', ms)
      legend('Target', 'PPI $$\mathcal{P}_{5}$$', 'PPI $$\mathcal{P}_{7}$$', 'Interpreter', 'latex')
      ylim([0.7 2])
      hold off
      %title('Cloud Mixing Ratio')
      xlabel('g/kg')
      ylabel('z (km)')
      set(gca, 'FontSize', fs)
      axes('Position',[.7 .7 .2 .2])
      box on
      plot(DW600(s:e,6)*1e+3, DW600(s:e, 1)*1e-3, '-k', 'LineWidth', lw, 'MarkerSize', ms)
      hold on
      plot(DPPI5S1(s:e,6)*1e+3,  DPPI5S1(s:e, 1)*1e-3,  '-b',  'LineWidth', lw, 'MarkerSize', ms)
      plot(DPPI7S1(s:e,6)*1e+3,  DPPI7S1(s:e, 1)*1e-3,  '-r',  'LineWidth', lw, 'MarkerSize', ms)
      ylim([1.4 1.45])
      xlim([-0.005, 0.005])
      hold off
      set(gca, 'FontSize', fs2)

      figure;clf
      plot(DW600(s:e,6)*1e+3, DW600(s:e, 1)*1e-3, '-k', 'LineWidth', lw, 'MarkerSize', ms)
      hold on
      plot(DPPI5S2(s:e,6)*1e+3,  DPPI5S2(s:e, 1)*1e-3,  '-b',  'LineWidth', lw, 'MarkerSize', ms)
      plot(DPPI7S2(s:e,6)*1e+3,  DPPI7S2(s:e, 1)*1e-3,  '-r',  'LineWidth', lw, 'MarkerSize', ms)
      legend('Target', 'PPI $$\mathcal{P}_{5}$$', 'PPI $$\mathcal{P}_{7}$$', 'Interpreter', 'latex')
      ylim([0.7 2])
      hold off
      %title('Cloud Mixing Ratio')
      xlabel('g/kg')
      ylabel('z (km)')
      set(gca, 'FontSize', fs)
      axes('Position',[.7 .7 .2 .2])
      box on
      plot(DW600(s:e,6)*1e+3, DW600(s:e, 1)*1e-3, '-k', 'LineWidth', lw, 'MarkerSize', ms)
      hold on
      plot(DPPI5S2(s:e,6)*1e+3,  DPPI5S2(s:e, 1)*1e-3,  '-b',  'LineWidth', lw, 'MarkerSize', ms)
      plot(DPPI7S2(s:e,6)*1e+3,  DPPI7S2(s:e, 1)*1e-3,  '-r',  'LineWidth', lw, 'MarkerSize', ms)
      ylim([1.5 1.9])
      xlim([-0.005, 0.005])
      hold off
      set(gca, 'FontSize', fs2)

      figure;clf
      plot(DW600(s:e,6)*1e+3, DW600(s:e, 1)*1e-3, '-k', 'LineWidth', lw, 'MarkerSize', ms)
      hold on
      plot(DPPI5S3(s:e,6)*1e+3,  DPPI5S3(s:e, 1)*1e-3,  '-b',  'LineWidth', lw, 'MarkerSize', ms)
      plot(DPPI7S3(s:e,6)*1e+3,  DPPI7S3(s:e, 1)*1e-3,  '-r',  'LineWidth', lw, 'MarkerSize', ms)
      legend('Target', 'PPI $$\mathcal{P}_{5}$$', 'PPI $$\mathcal{P}_{7}$$', 'Interpreter', 'latex')
      ylim([0.7 3])
      hold off
      %title('Cloud Mixing Ratio')
      xlabel('g/kg')
      ylabel('z (km)')
      set(gca, 'FontSize', fs)
      axes('Position',[.7 .7 .2 .2])
      box on
      plot(DW600(s:e,6)*1e+3, DW600(s:e, 1)*1e-3, '-k', 'LineWidth', lw, 'MarkerSize', ms)
      hold on
      plot(DPPI5S3(s:e,6)*1e+3,  DPPI5S3(s:e, 1)*1e-3,  '-b',  'LineWidth', lw, 'MarkerSize', ms)
      plot(DPPI7S3(s:e,6)*1e+3,  DPPI7S3(s:e, 1)*1e-3,  '-r',  'LineWidth', lw, 'MarkerSize', ms)
      ylim([1.8 2.5])
      xlim([-0.005, 0.005])
      hold off
      set(gca, 'FontSize', fs2)

      max_qc_pchip = max(DPCHIP(s:e, 6)*1e+3);
      max_qc_clip  = max(DCLIP(s:e, 6)*1e+3);
      max_qc_std   = max(DSTD(s:e, 6)*1e+3);
      max_qc       = max(DW600(s:e, 6)*1e+3);

      max_qc_ppi5_s1 = max(DPPI5S1(s:e, 6)*1e+3);
      max_qc_ppi7_s1 = max(DPPI7S1(s:e, 6)*1e+3);
      max_qc_ppi5_s2 = max(DPPI5S2(s:e, 6)*1e+3);
      max_qc_ppi7_s2 = max(DPPI7S2(s:e, 6)*1e+3);
      max_qc_ppi5_s3 = max(DPPI5S3(s:e, 6)*1e+3);
      max_qc_ppi7_s3 = max(DPPI7S3(s:e, 6)*1e+3);

      max_qc_dbi5_s1 = max(DDBI5S1(s:e, 6)*1e+3);
      max_qc_dbi7_s1 = max(DDBI7S1(s:e, 6)*1e+3);
      max_qc_dbi5_s2 = max(DDBI5S2(s:e, 6)*1e+3);
      max_qc_dbi7_s2 = max(DDBI7S2(s:e, 6)*1e+3);
      max_qc_dbi5_s3 = max(DDBI5S3(s:e, 6)*1e+3);
      max_qc_dbi7_s3 = max(DDBI7S3(s:e, 6)*1e+3);


      qc_pchip = trapz(DPCHIP(s:e, 1), DPCHIP(s:e, 6)*1e+3);
      qc_clip  = trapz(DCLIP(s:e, 1),  DCLIP(s:e, 6)*1e+3);
      qc_std   = trapz(DSTD(s:e, 1),   DSTD(s:e, 6)*1e+3);
      qc       = trapz(DW600(s:e, 1),  DW600(s:e, 6)*1e+3);

      qc_ppi5_s1 = trapz(DPPI5S1(s:e, 1), DPPI5S1(s:e, 6)*1e+3);
      qc_ppi7_s1 = trapz(DPPI7S1(s:e, 1), DPPI7S1(s:e, 6)*1e+3);
      qc_ppi5_s2 = trapz(DPPI5S2(s:e, 1), DPPI5S2(s:e, 6)*1e+3);
      qc_ppi7_s2 = trapz(DPPI7S2(s:e, 1), DPPI7S2(s:e, 6)*1e+3);
      qc_ppi5_s3 = trapz(DPPI5S1(s:e, 1), DPPI5S3(s:e, 6)*1e+3);
      qc_ppi7_s3 = trapz(DPPI7S1(s:e, 1), DPPI7S3(s:e, 6)*1e+3);

      qc_dbi5_s1 = trapz(DDBI5S1(s:e, 1), DDBI5S1(s:e, 6)*1e+3);
      qc_dbi7_s1 = trapz(DDBI7S1(s:e, 1), DDBI7S1(s:e, 6)*1e+3);
      qc_dbi5_s2 = trapz(DDBI5S2(s:e, 1), DDBI5S2(s:e, 6)*1e+3);
      qc_dbi7_s2 = trapz(DDBI7S2(s:e, 1), DDBI7S2(s:e, 6)*1e+3);
      qc_dbi5_s3 = trapz(DDBI5S3(s:e, 1), DDBI5S3(s:e, 6)*1e+3);
      qc_dbi7_s3 = trapz(DDBI7S3(s:e, 1), DDBI7S3(s:e, 6)*1e+3);

     
     fprintf(fileID,  '%.2f   & %.2f   & %.2f   & %.2f  \n', max_qc, max_qc_std, max_qc_clip, max_qc_pchip)
     fprintf(fileID,  '%.2f   & %.2f   & %.2f   & %.2f  \n', qc,     qc_std,     qc_clip,     qc_pchip)
     fprintf(fileID, '\n\n') 
    
     fprintf(fileID,  '%.2f   & %.2f   & %.2f   & %.2f  & %.2f   & %.2f   & %.2f  & %.2f   & %.2f \n', ...
              max_qc_dbi5_s1, max_qc_dbi7_s1, max_qc_dbi5_s2, max_qc_dbi7_s2, max_qc_dbi5_s3, max_qc_dbi7_s3)
     fprintf(fileID,  '%.2f   & %.2f   & %.2f   & %.2f  & %.2f   & %.2f   & %.2f  & %.2f   & %.2f \n', ...
              qc_dbi5_s1, qc_dbi7_s1, qc_dbi5_s2, qc_dbi7_s2, qc_dbi5_s3, qc_dbi7_s3)
     fprintf(fileID, '\n') 
     fprintf(fileID,  '%.2f   & %.2f   & %.2f   & %.2f  & %.2f   & %.2f   & %.2f  & %.2f   & %.2f \n', ...
              max_qc_ppi5_s1, max_qc_ppi7_s1, max_qc_ppi5_s2, max_qc_ppi7_s2, max_qc_ppi5_s3, max_qc_ppi7_s3)
     fprintf(fileID,  '%.2f   & %.2f   & %.2f   & %.2f  & %.2f   & %.2f   & %.2f  & %.2f   & %.2f \n', ...
              qc_ppi5_s1, qc_ppi7_s1, qc_ppi5_s2, qc_ppi7_s2, qc_ppi5_s3, qc_ppi7_s3)

     if(qc_pchip > 2*qc)
       fprintf(fileID, 'qc_ pchip is  %.2f  qc \n',  qc_pchicp/qc)
     else
       fprintf(fileID, 'qc_ pchip is  %.2f % more than  qc \n',  qc_pchip/qc*100- 100)
     end 
     %
     if(qc_std > 2*qc)
       fprintf(fileID, 'qc_std is  %.2f  qc \n',  qc_std/qc)
     else
       fprintf(fileID, 'qc_std is  %.2f % more than  qc \n',  qc_std/qc*100- 100)
     end 
     %
     if(qc_clip > 2*qc)
       fprintf(fileID, 'qc_clip is  %.2f  qc \n',  qc_clip/qc)
     else
       fprintf(fileID, 'qc_clip is  %.2f % more than  qc \n',  qc_clip/qc*100- 100)
     end 
     %
     if(qc_dbi5_s1 > 2*qc)
       fprintf(fileID, 'qc_dbi5_s1 is  %.2f  qc \n',  qc_dbi5_s1/qc)
     else
       fprintf(fileID, 'qc_dbi5_s1 is  %.2f % more than  qc \n',  qc_dbi5_s1/qc*100- 100)
     end 
     %
     if(qc_dbi7_s1 > 2*qc)
       fprintf(fileID, 'qc_dbi7_s1 is  %.2f  qc \n',  qc_dbi7_s1/qc)
     else
       fprintf(fileID, 'qc_dbi7_s1 is  %.2f % more than  qc \n',  qc_dbi7_s1/qc*100- 100)
     end
     %
     if(qc_dbi5_s2 > 2*qc)
       fprintf(fileID, 'qc_dbi5_s2 is  %.2f  qc \n',  qc_dbi5_s2/qc)
     else
       fprintf(fileID, 'qc_dbi5_s2 is  %.2f % more than  qc \n',  qc_dbi5_s2/qc*100- 100)
     end 
     %
     if(qc_dbi7_s2 > 2*qc)
       fprintf(fileID, 'qc_dbi7_s2 is  %.2f  qc \n',  qc_dbi7_s2/qc)
     else
       fprintf(fileID, 'qc_dbi7_s2 is  %.2f % more than  qc \n',  qc_dbi7_s2/qc*100- 100)
     end
     %
     if(qc_dbi5_s3 > 2*qc)
       fprintf(fileID, 'qc_dbi5_s3 is  %.2f  qc \n',  qc_dbi5_s3/qc)
     else
       fprintf(fileID, 'qc_dbi5_s3 is  %.2f % more than  qc \n',  qc_dbi5_s3/qc*100- 100)
     end 
     %
     if(qc_dbi7_s3 > 2*qc)
       fprintf(fileID, 'qc_dbi7_s3 is  %.2f  qc \n',  qc_dbi7_s3/qc)
     else
       fprintf(fileID, 'qc_dbi7_s3 is  %.2f % more than  qc \n',  qc_dbi7_s3/qc*100- 100)
     end
     %
     if(qc_ppi5_s1 > 2*qc)
       fprintf(fileID, 'qc_ppi5_s1 is  %.2f  qc \n',  qc_ppi5_s1/qc)
     else
       fprintf(fileID, 'qc_ppi5_s1 is  %.2f % more than  qc \n',  qc_ppi5_s1/qc*100- 100)
     end 
     %
     if(qc_ppi7_s1 > 2*qc)
       fprintf(fileID, 'qc_ppi7_s1 is  %.2f  qc \n',  qc_ppi7_s1/qc)
     else
       fprintf(fileID, 'qc_ppi7_s1 is  %.2f % more than  qc \n',  qc_ppi7_s1/qc*100- 100)
     end
     %
     if(qc_dbi5_s2 > 2*qc)
       fprintf(fileID, 'qc_ppi5_s2 is  %.2f  qc \n',  qc_ppi5_s2/qc)
     else
       fprintf(fileID, 'qc_ppi5_s2 is  %.2f % more than  qc \n',  qc_ppi5_s2/qc*100- 100)
     end 
     %
     if(qc_ppi7_s2 > 2*qc)
       fprintf(fileID, 'qc_ppi7_s2 is  %.2f  qc \n',  qc_ppi7_s2/qc)
     else
       fprintf(fileID, 'qc_ppi7_s2 is  %.2f % more than  qc \n',  qc_ppi7_s2/qc*100- 100)
     end
     %
     if(qc_ppi5_s3 > 2*qc)
       fprintf(fileID, 'qc_ppi5_s3 is  %.2f  qc \n',  qc_ppi5_s3/qc)
     else
       fprintf(fileID, 'qc_ppi5_s3 is  %.2f % more than  qc \n',  qc_ppi5_s3/qc*100- 100)
     end 
     %
     if(qc_dbi7_s3 > 2*qc)
       fprintf(fileID, 'qc_ppi7_s3 is  %.2f  qc \n',  qc_ppi7_s3/qc)
     else
       fprintf(fileID, 'qc_ppi7_s3 is  %.2f % more than  qc \n',  qc_ppi7_s3/qc*100- 100)
     end






end


fclose(fileID);

%
%
%%
