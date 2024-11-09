function [f SP]=data_spectrum(data,fc,tag1,tag2)

% data:(channel*sample) si possono dare più segnali contemporaneamente in input
%per fare ciò si deve organizzare la matrice dei dati avendo i sample nelle
%colonne e i canali (cioè i segnali) nelle righe
% fc sample frequency
% tag1: ['e'] energy density spectrum, ['p'] power density spectrum
% tag2: [0] no-mean correction, [1] mean correction

SP=[];

if tag2==1
    %data=detrend(data','constant');
    data=data-mean(data,2).*ones(size(data));
    %data=(data')*(10^-3);
end

for i=1:size(data,1)
    
    l=size(data,2);
   
    switch tag1
        
        case 'e'
            y=fft(data(i,:)); %basic step 1
            %ritorna una serie di numeri complessi (a+jb) di cui dobbiamo
            %fare il modulo quadro (|a^2+b^2|)
         
            y(1)=[]; 
            %removing the first sample of the fft
            %is useful to remove the continuos, so the mean
            n=floor(length(y)/2); %basic step 2
            %dimezziamo il numero di campioni per shannon, nelle freq sono
            %importanti solo la prima metà dei campioni, floor dovuto a
            %quando abbiamo un numero dispari di campioni
            
            % spectrum amplitude
            %m = abs(y(1:n)).^2;
            m = (abs(y(1:n))).^2/n; %basic step 3
            %fa modulo solo dei campioni da 1 fino a n cioè la prima metà
            %modulo quadro, ecco perché dopo aver fatto abs fa .^2, cioè eleva al quadrato il vettore.                                                            
            % Invece il /n serve a normalizzare l'ampiezza rispetto alla lunghezza del segnale 
            % (perché se ho più segnali (ma con la stessa freq di campionamento, quindi segnale più lungo) l'ampiezza viene maggiore)
            SP=[SP, m'];
            
        case 'p'
            %by computing the fourier transform of the (auto)correlation function
            %we obtain the power spectrum
            cc=xcorr(data(i,:),data(i,:));
            y=fft(cc);
            
            y(1)=[];
            n=floor(length(y)/2);
            
            % spectrum amplitude
            %m = abs(y(1:n)).^2;
            %m=abs(y(n/2:n));
            m = abs(y(1:n))/n;
            SP=[SP, m']; %SP è una matrice ricorsiva che serve
            %a fare lo spettro di più segnali
            %quando si passa un segnale solo, contiene i valori dell'asse y
            %cioè i moduli, le ampiezze
            
    end
    
end

% frequency window
f =((1:n)*fc/2)/n; %basic step 4
%la freq va da 1 fino a n (metà dei campioni) 
%moltiplichiamo per la sample freq e dividiamo per due, per il teo di shannon e poi
%dividiamo per n, come riportato nelle slide.

figure
plot(f,SP)
ylabel('Spectrum');
xlabel('Hertz(Hz)');
xlim([0 max(f)]);
%xlim([0 3]);
ylim([min(min(SP)) max(max(SP))]);
 %ylim([0 100]);

end

%per utilizzare dataspectrum con più segnali in colonna
%DATA=[ys;yp];    2x500 (2 segnali da 500 campioni per esempio)
%[f SP]=data_spectrum(DATA,Fs,'e',1)