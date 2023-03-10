    	    
% Generisanje asinhrnonog signala od vise korisnika (MS/UE/Terminala), 
% sa vise Data kanala po korisniku i prijem na strani BS
%periodogram

clear all
close all
clc

%rand(10013,1);

Nuser = 8;          % Broj korisnika
ChByUser = 2;       % Broj Data/Signaling kanala po korisniku
ChByI = 1;          % Broj Data/Signaling kanala u I-grani
ChByQ = 1;          % Broj Data/Signaling kanala u Q-grani
MAIPowerVector(1:Nuser,1) = [0.25 0.125 0.0625 0.25 0.0625 0.0375 0.125 0.875];  % Raspodela snage po korisnicima na ulazu BS - Default vrednost je za slucaj idealne kontrole snage
                                                
NoDemodUser = 3;    % Redni broj korisnika ciji se signal demodulise                                            
                                

ChipByFrame = 16000;            % Broj chip-ova po frejmu 
FramePeriod = 1;             % Trajanje frejma - 1s
FrameRate = 1/FramePeriod;      % Protok frejmova - 100 frames/s
Vchip = FrameRate*ChipByFrame;  % Protok chip-ova 3.84Mchip/s
Tchip = 1/Vchip;                % Trajanje chip-a
TimeSlotsByFrame = 5;          % Broj vremenskih slotova u frejmu

DataBitsByTimeSlotVector(1:3,1) = [50 100 200];                % Broj bita u kanalu po vremenskom slotu za svih 7 brzina signaliziranja
DataBitsByFrameVector(1:3,1) = TimeSlotsByFrame*DataBitsByTimeSlotVector;   % Broj bita u kanalu po frejmu za svih 7 brzina signaliziranja
DataRateVector(1:3,1) = DataBitsByFrameVector*FrameRate;    % Bitski protok po kanalu za svih 7 brzina signaliziranja
SpreadingFactorVector(1:3,1) = [64 32 16];      % Spreading Factor (SF) za korisnika - za svih 7 brzina signaliziranja

BitRatesByUserFlag(1:Nuser,1) = [3 2 1 3 1 2 1 2];    % Zadavanje Protoka/SF po korisniku - redni broj u vektorima iznad

BitsByFrame(1:Nuser,1) = DataBitsByFrameVector(BitRatesByUserFlag(1:Nuser,1));  % Broj bita u svakom kanalu po frejmu za svakog korisnika
ChipByBit(1:Nuser,1) = ChipByFrame./BitsByFrame(1:Nuser,1);                     % Broj chip-ova po bitu za svakog korisnika

Nframe = 100;     % Ukupan broj frejmova koji se prenosi

% Generisanje Channelizing (OVSF) sekvenci za pojedine kanale za korisnike 

ChannelizingCodesMatrix(1:Nuser,1:ChipByFrame*Nframe,1:ChByUser) = zeros(Nuser,ChipByFrame*Nframe,ChByUser);
index = 0;
for user = 1:Nuser
    % Generisanje ortogonalnog skupa OVSF sekvenci - zavise od broja kanala za korisnika
    for chan = 1:ChByUser
        SF = SpreadingFactorVector(BitRatesByUserFlag(user,1),1);
        hOVSF = comm.OVSFCode('SamplesPerFrame',ChipByFrame*Nframe,'SpreadingFactor',SF,'Index',index);    
        seq = 0;
        index = index + 1;
        seq(1:ChipByFrame*Nframe,1) = step(hOVSF);
        ChannelizingCodesMatrix(user,1:ChipByFrame*Nframe,chan) = seq;
    end;    
end;
clear hgld user scs scsmatrix chan hOVSF seq

% Generisanje Data/Signaling bita - ne posmatramo protokole pa signalzacione bite generisemo slucajno

% Matrica bita po korisnicima i kanalima na nivou chip-a
DataBitChipMatrix(1:Nuser,1:ChipByFrame*Nframe,1:ChByUser) = zeros(Nuser,ChipByFrame*Nframe,ChByUser);    

for user = 1:Nuser
    data = 0;
    data(1:BitsByFrame(user,1)*Nframe,1:ChByUser) = (rand(BitsByFrame(user,1)*Nframe,ChByUser)>0.5);
    pattern(1:ChipByBit(user,1),1) = ones(ChipByBit(user,1),1);
    for chan = 1:ChByUser
        chipovi = 0; biti = 0;         
        biti(1:ChipByBit(user,1):ChipByFrame*Nframe,1) = data(1:BitsByFrame(user)*Nframe,chan);
        chipovi = conv(biti,pattern);
        DataBitChipMatrix(user,1:ChipByFrame*Nframe,chan) = chipovi(1:ChipByFrame*Nframe,1);        
    end;    
end;

tosafrejm = (0:Tchip*ChipByBit(NoDemodUser,1):(ChipByFrame-1)*Tchip)/Tchip/ChipByBit(NoDemodUser,1);
figure(11)
hold on
stem(tosafrejm,DataBitChipMatrix(NoDemodUser,1:ChipByBit(NoDemodUser,1):ChipByFrame,1),'bo')
stem(tosafrejm,DataBitChipMatrix(NoDemodUser,1:ChipByBit(NoDemodUser,1):ChipByFrame,2),'rx')
axis([0 BitsByFrame(NoDemodUser,1) -0.5 1.2]);
xlabel('Vreme [t/Tbita]')
ylabel('Vremenski oblik signala [V]')
legend('Informacioni biti - Kanal 1', 'Informacioni biti - Kanal 2')
grid on
hold off

clear user chan data pattern chipovi biti

% Modulacija i sirenje spektra za svakog korisnika

Symbols(1:Nuser,1:ChipByFrame*Nframe) = zeros(Nuser,ChipByFrame*Nframe);

for user = 1:Nuser
    Ichips(1:ChipByFrame*Nframe) = zeros(ChipByFrame*Nframe,1);
    Qchips(1:ChipByFrame*Nframe) = zeros(ChipByFrame*Nframe,1);
    
    for ichan = 1:ChByI
        Ichips(1:ChipByFrame*Nframe) = Ichips(1:ChipByFrame*Nframe)+ (2*DataBitChipMatrix(user,1:ChipByFrame*Nframe,ichan)-1).*ChannelizingCodesMatrix(user,1:ChipByFrame*Nframe,ichan);
    end;
    for qchan = 1:ChByQ
        Qchips(1:ChipByFrame*Nframe) = Qchips(1:ChipByFrame*Nframe)+ (2*DataBitChipMatrix(user,1:ChipByFrame*Nframe,qchan+ChByI)-1).*ChannelizingCodesMatrix(user,1:ChipByFrame*Nframe,qchan+ChByI);
    end;
    Symbols(user,1:ChipByFrame*Nframe) = Ichips(1:ChipByFrame*Nframe)+1i*Qchips(1:ChipByFrame*Nframe);
end;

tosachip = (0:Tchip:8*ChipByFrame/BitsByFrame(NoDemodUser,1)*Tchip-Tchip)/Tchip;

figure(13)
plot(tosachip,real(Symbols(NoDemodUser,1:8*ChipByFrame/BitsByFrame(NoDemodUser,1))),'b');
axis([0 8*ChipByFrame/BitsByFrame(NoDemodUser,1)-1 -2.5 2.5]);
xlabel('Vreme [t/Tchipa]')
ylabel('Vrednosti nakon mnozenja sa OVSF kodovima [V]')
legend('Kanal 1')
grid on
hold off

figure(14)
plot(tosachip,imag(Symbols(NoDemodUser,1:8*ChipByFrame/BitsByFrame(NoDemodUser,1))),'b');
axis([0 8*ChipByFrame/BitsByFrame(NoDemodUser,1)-1 -2.5 2.5]);
xlabel('Vreme [t/Tchipa]')
ylabel('Vrednosti nakon mnozenja sa OVSF kodovima [V]')
legend('Kanal 2')
grid on
hold off

clear ichan qchan

Prx = mean(abs(sqrt(MAIPowerVector(NoDemodUser,1))));
PrxdB = 10*log10(Prx);
%-------------------------------------------------------------------------------------------------------------


% Signal na ulazu u prijemnik BS
Signal(1,1:ChipByFrame*Nframe) = zeros(1,ChipByFrame*Nframe);

tau_vector(1:Nuser,1) = zeros(Nuser,1);      % Vektor kasnjenja signala svih korisnika u chip-ovima
teta_vector(1:Nuser,1) = pi/4*zeros(Nuser,1); % Fazni offset signala korisnika (za zeljenog je greska sinhrnizacije faze)

PhaseOffset(1:Nuser,1) = exp(-1i*teta_vector(1:Nuser,1)); % Complex fazni offset signala korisnika (za zeljenog je greska sinhrnizacije faze)
Nuser = 8;

for user = 1:Nuser
    Signal(1,(tau_vector(user,1)+1):ChipByFrame*Nframe) = Signal(1,(tau_vector(user,1)+1):ChipByFrame*Nframe) + sqrt(2*MAIPowerVector(user,1))*Symbols(user,(tau_vector(user,1)+1):ChipByFrame*Nframe)*PhaseOffset(user,1);
end;

%Signal(1,ChipByFrame*Nframe/4:ChipByFrame*Nframe-ChipByFrame*Nframe/16) = Signal(1,ChipByFrame*Nframe/4:ChipByFrame*Nframe-ChipByFrame*Nframe/4) + 1;
%U liniji iznad se dodaje spoljna interferencija u određenom segmentu.


% Sum na ulazu u prijemnik izrazen kao odnos S/N za signal zeljenog korisnika (ne za celokupni signal)
PeVector(1:4,1) = [1/2 0.3 0.001 0];
PeFlag = 3;  %izbor verovatnoce greske
SNRdB = 10*log10(erfcinv(2*PeVector(PeFlag)));
NoiseLevel_dBw = PrxdB - SNRdB + 10*log10(ChipByBit(NoDemodUser,1)/2);       % Nivo kompleksnog AWGN na ulazu u prijemnik
NoiseLevel_dBw = 10;
SNRdB = PrxdB - NoiseLevel_dBw + 10*log10(ChipByBit(NoDemodUser,1)/2);
ComplexNoise(1,1:ChipByFrame*Nframe) = wgn(1,ChipByFrame*Nframe, NoiseLevel_dBw, 1, 'dBw', 'complex');  % Odbirci kompleksnog AWGN na prijemu
%ComplexNoise(1,1:ChipByFrame*Nframe) = zeros(1,ChipByFrame*Nframe);
% Ukupan signal + AWGN
SignalAWGN(1,1:ChipByFrame*Nframe)= Signal(1,1:ChipByFrame*Nframe) + ComplexNoise(1,1:ChipByFrame*Nframe);

% Demodulacija signala #1 korisnika za sve kanale

% Kompresija sa Screambling sekvencom - CDMA razdvajanje signala vise korisnika
SigDespreaded_I(1,1:ChipByFrame*Nframe) = real(SignalAWGN(1,1:ChipByFrame*Nframe));
SigDespreaded_Q(1,1:ChipByFrame*Nframe) = imag(SignalAWGN(1,1:ChipByFrame*Nframe));

% Izdvajanje kanala koriscenjem Channalizing sekvenci - CDM demultipleksiranje kanala istog korisnika
Channels(1:ChipByFrame*Nframe,1:ChByUser) = zeros(ChipByFrame*Nframe,ChByUser);
for ichan = 1:ChByI
    Channels(1:ChipByFrame*Nframe,ichan) = SigDespreaded_I(1,1:ChipByFrame*Nframe).*ChannelizingCodesMatrix(NoDemodUser,1:ChipByFrame*Nframe,ichan);
end;
for qchan = 1:ChByQ
    Channels(1:ChipByFrame*Nframe,qchan+ChByI) = SigDespreaded_Q(1,1:ChipByFrame*Nframe).*ChannelizingCodesMatrix(NoDemodUser,1:ChipByFrame*Nframe,qchan+ChByI);
end;

% Integrator sa rasterecenjem po periodu signaliziranja
ErrorsByChan(1:ChByUser,1) = zeros(ChByUser,1);
DemodBits(1:BitsByFrame(NoDemodUser)*Nframe,1:ChByUser) = zeros(BitsByFrame(NoDemodUser)*Nframe,ChByUser);

for chan = 1:ChByUser
    matrica = 0;
    matrica(1:ChipByBit(NoDemodUser,1),1:BitsByFrame(NoDemodUser)*Nframe) = reshape(Channels(1:ChipByFrame*Nframe,chan),ChipByBit(NoDemodUser,1),BitsByFrame(NoDemodUser)*Nframe);
    DemodBits(1:BitsByFrame(NoDemodUser)*Nframe,chan) = (sum(matrica)>0);
    ErrorsByChan(chan,1) = sum(abs((DemodBits(1:BitsByFrame(NoDemodUser)*Nframe,chan))'-DataBitChipMatrix(NoDemodUser,1:ChipByBit(NoDemodUser,1):ChipByFrame*Nframe,chan)));
end;
clear matrica

Nuser
ChByUser
MAIPowerVector'                              
ChipByFrame
FramePeriod
FrameRate
Vchip
Tchip
TimeSlotsByFrame

DataBitsByTimeSlot = DataBitsByTimeSlotVector(BitRatesByUserFlag)'
DataBitsByFrame = DataBitsByFrameVector(BitRatesByUserFlag)'
DataRate = DataRateVector(BitRatesByUserFlag)'
SpreadingFactor = SpreadingFactorVector(BitRatesByUserFlag)'

'Verovatnoca greske po kanalima 1 - 2'
BER = (ErrorsByChan/BitsByFrame(NoDemodUser,1)/Nframe)'
% Prikaz spektara i vremenskih oblika pojedinih signala pri prijemu signala NoDemodUser
SignalOne(1,1:ChipByFrame*Nframe) = zeros(1,ChipByFrame*Nframe);
close all
for user = 1:Nuser
 
    
    'Korisnik #: '
    user
    
    SignalOne(1,(tau_vector(user,1)+1):ChipByFrame*Nframe)= sqrt(MAIPowerVector(user,1))*Symbols(user,(tau_vector(user,1)+1):ChipByFrame*Nframe)*PhaseOffset(user,1);
    tosachip = (0:Tchip:8*ChipByFrame/BitsByFrame(NoDemodUser,1)*Tchip-Tchip)/Tchip;
    figure(1)
    plot(tosachip,real(SignalOne(1:8*ChipByFrame/BitsByFrame(NoDemodUser))),'b')
    axis([0 8*ChipByFrame/BitsByFrame(NoDemodUser,1)-1 -inf inf])
    xlabel('Vreme [t/Tchipa]')
    ylabel('Vremenski oblik signala [V]')
    grid on
    hold off
    
    figure(2)
    plot(tosachip,imag(SignalOne(1:8*ChipByFrame/BitsByFrame(NoDemodUser))),'b')
    axis([0 8*ChipByFrame/BitsByFrame(NoDemodUser,1)-1 -inf inf])
    xlabel('Vreme [t/Tchipa]')
    ylabel('Vremenski oblik signala [V]')
    grid on
    hold off
   
    % Kompresija sa Screambling sekvencom - CDMA razdvajanje signala vise korisnika
    SigOneDespreaded_I(1,1:ChipByFrame*Nframe) = real(SignalOne(1,1:ChipByFrame*Nframe));
    SigOneDespreaded_Q(1,1:ChipByFrame*Nframe) = imag(SignalOne(1,1:ChipByFrame*Nframe));
    
    % Izdvajanje kanala koriscenjem Channalizing sekvenci - CDM demultipleksiranje kanala istog korisnika
    ChannelsOne(1:ChipByFrame*Nframe,1:ChByUser) = zeros(ChipByFrame*Nframe,ChByUser);
    for ichan = 1:ChByI
        ChannelsOne(1:ChipByFrame*Nframe,ichan) = SigOneDespreaded_I(1,1:ChipByFrame*Nframe).*ChannelizingCodesMatrix(NoDemodUser,1:ChipByFrame*Nframe,ichan);
    end;
    for qchan = 1:ChByQ
        ChannelsOne(1:ChipByFrame*Nframe,qchan+ChByI) = SigOneDespreaded_Q(1,1:ChipByFrame*Nframe).*ChannelizingCodesMatrix(NoDemodUser,1:ChipByFrame*Nframe,qchan+ChByI);
    end;
    
    figure(5)
    plot(tosachip,ChannelsOne(1:8*ChipByFrame/BitsByFrame(NoDemodUser),1),'b')
    axis([0 8*ChipByFrame/BitsByFrame(NoDemodUser,1)-1 -inf inf])
    xlabel('Vreme [t/Tchipa]')
    ylabel('Vremenski oblici nakon mnozenja sa Channalizing sekvencama - I grana [V]')
    grid on
    hold off
    
    figure(6)
    plot(tosachip,ChannelsOne(1:8*ChipByFrame/BitsByFrame(NoDemodUser),2),'b')
    axis([0 8*ChipByFrame/BitsByFrame(NoDemodUser,1)-1 -inf inf])
    xlabel('Vreme [t/Tchipa]')
    ylabel('Vremenski oblici nakon mnozenja sa Channalizing sekvencama - Q grana [V]')
    grid on
    hold off
        
    DemodOne(1:BitsByFrame(NoDemodUser)*Nframe,1:ChByUser) = zeros(BitsByFrame(NoDemodUser)*Nframe,ChByUser);

    for chan = 1:ChByUser
        matrica = 0;
        matrica(1:ChipByBit(NoDemodUser,1),1:BitsByFrame(NoDemodUser)*Nframe) = reshape(ChannelsOne(1:ChipByFrame*Nframe,chan),ChipByBit(NoDemodUser,1),BitsByFrame(NoDemodUser)*Nframe);
        DemodOne(1:BitsByFrame(NoDemodUser)*Nframe,chan) = (sum(matrica) > 0);    
    end;
    
    tosafrejm = (0:Tchip*ChipByBit(NoDemodUser,1):(ChipByFrame-1)*Tchip)/Tchip/ChipByBit(NoDemodUser,1);
    figure(7)
    stem(tosafrejm, DemodOne(1:BitsByFrame,1),'bo')
    axis([0 BitsByFrame(NoDemodUser,1) -inf inf]);
    xlabel('Vreme [t/Tbita]')
    ylabel('Odbirci simbola nakon integracije - I grana [V]')
    grid on
    hold off
        
    figure(8)
    stem(tosafrejm, DemodOne(1:BitsByFrame,2),'bo')
    axis([0 BitsByFrame(NoDemodUser,1) -inf inf]);
    xlabel('Vreme [t/Tbita]')
    ylabel('Odbirci simbola nakon integracije - Q grana [V]')
    grid on
    hold off
    

    % Proracun spektra 

    Nfft = 4096;
    
    FT_ulaziI_periodogram(1:Nfft,1) = zeros(Nfft,1);
    FT_ulaziQ_periodogram(1:Nfft,1) = zeros(Nfft,1);
    
    n = max(size(SignalOne))/Nfft;
    for i = 1:n
        FT_ulazI(1:Nfft,1) = fft(real(SignalOne((i-1)*Nfft+1:i*Nfft)),Nfft);
        FT_ulazQ(1:Nfft,1) = fft(imag(SignalOne((i-1)*Nfft+1:i*Nfft)),Nfft);
        FT_ulaziI_periodogram = FT_ulaziI_periodogram + abs(FT_ulazI).^2/(Nfft);
        FT_ulaziQ_periodogram = FT_ulaziQ_periodogram + abs(FT_ulazQ).^2/(Nfft);
    end;
    
    FT_ulaziI_periodogram = FT_ulaziI_periodogram/n;
    FT_ulaziQ_periodogram = FT_ulaziQ_periodogram/n;

    fosa = (1:Nfft)/Tchip;
    f = Vchip*(0:1/Nfft:(Nfft-1)/Nfft);
    
    figure(20)
    plot(f,FT_ulaziI_periodogram,'b--')
    hold on
    plot(f,FT_ulaziQ_periodogram,'r-')
    title('Spektar signala na ulazu - I i Q grana')
    legend('SGSS na ulazu - I grana', 'SGSS na ulazu - Q grana')
    hold off
    
    
    FT_Channel1_periodogram(1:Nfft,1) = zeros(Nfft,1);
    FT_Channel2_periodogram(1:Nfft,1) = zeros(Nfft,1);
    
    n = max(size(ChannelsOne))/Nfft;
    for i = 1:n
        FT_Channel1_pom(1:Nfft,1) = fft(ChannelsOne((i-1)*Nfft+1:i*Nfft,1),Nfft);
        FT_Channel2_pom(1:Nfft,1) = fft(ChannelsOne((i-1)*Nfft+1:i*Nfft,2),Nfft);
        FT_Channel1_periodogram = FT_Channel1_periodogram + abs(FT_Channel1_pom).^2/(Nfft);
        FT_Channel2_periodogram = FT_Channel2_periodogram + abs(FT_Channel2_pom).^2/(Nfft);
    end;
    
    FT_Channel1_periodogram = FT_Channel1_periodogram/n;
    FT_Channel2_periodogram = FT_Channel2_periodogram/n;
    
    figure(22)
    plot(f,FT_Channel1_periodogram,'b-')
    hold on
    plot(f,FT_Channel2_periodogram,'c-')
    
    title('Spektar signala nakon kompresije Channalizing sekvencama zeljenog korisnika - za oba kanala')
    legend('Kanal 1', 'Kanal 2')
   
    hold off
   
end;