% visual inspection

clear

load nom.dat

s=size(nom);
no_files = s(1);


for i = 1:no_files
    
    
    yyyy = nom(i,1);
    mm = nom(i,2);
    dd = nom(i,3);
    hh = nom(i,4);
    mi  = nom(i,5);
    
    % transfert number into string for filename
    
    YYYY = num2str(yyyy);
    
    if mm < 10
        MM = ['0' num2str(mm)];
    else
        MM = num2str(mm);
    end
    
    if dd < 10
        DD = ['0' num2str(dd)];
    else
        DD = num2str(dd);
    end
        
    if hh < 10
        HH = ['0' num2str(hh)];
    else
        HH = num2str(hh);
    end

    if mi < 10
        MI = ['0' num2str(mi)];
    else
        MI = num2str(mi);
    end    
    
    fname = [YYYY '-' MM '-' DD '_' HH 'h' MI '.bin'];
    
    
    if (yyyy==1993 & mm > 5)
    
    file = load(fname);
    
    II = find(T<1);
    A = find(b==i); %indices de a correspondants au profil i
    length(A)
    DENS = sw_dens(S(a(A),i), T(a(A),i), pres(a(A))); %density for bin into CIL core
    H(i) = cp*sum(DENS.*T(a(A), i))*dz; %Hest content of CIL for this profile
    
    
    % visual test for validation!
    subplot(1,2,1)
    I = find(file(:,2)~=-99);
    plot(file(I,2), file(I,1))
    set(gca, 'ydir', 'reverse')
    xlabel('temperature')
    ylabel('depth(m)')
    axis([-2 10 0 280])
    title(fname)
    
    subplot(1,2,2)
    I = find(file(:,3)~=-99);
    plot(file(I,3), file(I,1))
    set(gca, 'ydir', 'reverse')
    xlabel('temperature')
    axis([22 34 0 280])
    pause
    end
    
end
