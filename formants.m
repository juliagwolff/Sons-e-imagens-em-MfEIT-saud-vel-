clear all 
close all

%fucntion [F T W] = formants(x,fs,frames,overlap)

[x,fs]=audioread('7saudavel_900k_0.wav');

frames = 55; %tempo janelamento em ms
overlap = 35; %tempo sobreposição em ms


x=x(:,1); % if the sound is stereo then reduce for mono
N = length(x);

nframe = round(frames  * fs / 1000); %convert time to samples
noverlap = round(overlap * fs / 1000);

 
w = hamming(nframe); %hamming window applied in the frames;
L = 2^nextpow2(nframe); %size of FFT

%% energy for threshold
framese = 20; %tempo janelamento em ms

nframee = round(framese  * fs / 1000); % convert ms to samples

we = hamming(nframee); %hamming window for the frames of energy;

jan=we.^2;
ener=x.^2;
En=conv(ener(:,1),jan);

winlen=nframee;

temp = (1:length(En))/fs;
threshold=max(En)/(db2mag(40)); %value of the threshold


%% Formants


nlpc = (fs/1000)+2;

pos = 1; i = 1;
while (pos+nframe < N)
     frame = x(pos:pos+nframe-1);
     cn1= frame(:) .* (w(:));
     cn2= fft(cn1, L);
     C(:,i) = cn2;
     [A, err] = lpc(cn1,nlpc);
        
     rts = roots(A);

     rts = rts(imag(rts)>=0);
     angz = atan2(imag(rts),real(rts));

     [frqs,indices] = sort(angz.*(fs/(2*pi)));
     bw = -1/2*(fs/(2*pi))*log(abs(rts(indices)));

     nn = 1;
     for kk = 1:length(frqs)
        if (frqs(kk) > 90 && bw(kk) <400)
            formant(nn) = frqs(kk);
                if formant(nn)>5000
                    formant(nn)=NaN;
                end
                if En(pos)<threshold
                formant(nn)=NaN;
                end
                
                F(nn,i)=formant(nn);
                W(nn,i)=bw(kk);
                nn = nn+1;
                    
        end
     end
    
     if(i==10)
        [pxx, freqp] = periodogram(cn1,[],4096,fs);
        cn1f = abs(fft(cn1,L));
        freqf=(0:(L/2)-1)*(fs/L);
        
        %figure
        est = filter([0 -A(2:end)],1,cn1);
        
        

        

     end
     
     
     
     
     pos = pos + (nframe - noverlap);
     i = i + 1;
end

T = (round(nframe/2):(nframe-noverlap):N-1-round(nframe/2))/fs;







%D =(abs(C((temp1:temp2),(1:end))));
G = (abs(C));

G1 = log(G(1:end/2,:)) - 1000;
G2 = G1(1,:);

espectrogram = surf(T,freqf,G1,'EdgeColor','interp');
axis xy; 
axis tight;
colormap('gray')
view(0,90);
title('Espectrograma')
%set(gca,'Ydir','reverse')
xlabel('Tempo(s)')
ylabel('Frequ�ncia (s)')
zlabel('Amplitude')
ylim([0 fs/2])
uistack(espectrogram,'bottom')

