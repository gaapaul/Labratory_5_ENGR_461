%lab5 461
clc
%clear all
close all

%NxN image and RGB (3 colours) 
N = 128;
RGB = 3;
k  = 1.38e-23
T = 300;
snr = 10;
err_arr = [];
err_arr2 = [];

qfun = zeros(1,15);
M=4;
No = k*T;



%N,N,R array image
image = imread('lena128.jpg');
%makes single col array with 8 bit rows (strings)
image_bin = dec2bin(image);
%single bit col array (string)
%image_bin_rehape = reshape(image_bin,[],1);
sizeofimage = size(image_bin);
c =  zeros(sizeofimage(1) * sizeofimage(2)/2,2);
for i = 1:sizeofimage(1) 
  for j = 1: sizeofimage(2)/2
     c(4*(i-1)+j,:) =  image_bin(i,((j-1)*2+1:j*2));
  end
end
%2PAM
image_bin_rehape2 = reshape(image_bin,[],1);
image_num = str2num(image_bin_rehape2);
%4PAM
image_bin_rehape = char(c); 
noisy2 = normrnd(0,sqrt(No/2),length(image_bin_rehape2),1);
noisy = noisy2(1:length(image_bin_rehape));


for z = 1:11
  
    snr = z-1
    eb=No * (10^(snr/10));
    es = eb*log2(M);
    d=sqrt((es*3)/((M^2)-1));

    snr_check = 10*log10(eb/No)
    %4pam
    scaled_image_pam = FourPAMmodem(image_bin_rehape,d);
    img_pam_err=scaled_image_pam+noisy;
    %2pam
    scaled_image_pam2 = pammod(image_num,2)*sqrt(eb);
    img_pam_err2=scaled_image_pam2'+noisy2;
    %
    %4
    img_pam_demod = FourPAMdemodem(img_pam_err,d);
    err_cnt = 0;
    for i=1:length(c)
      %for j = 1:2
          if(img_pam_demod(i,1) ~= image_bin_rehape(i,1))
          err_cnt = err_cnt + 1;
          end
         if(img_pam_demod(i,2) ~= image_bin_rehape(i,2))
          err_cnt = err_cnt + 1;
        end
      %end
    end
    %
    %2
    img_pam_demod2 = pamdemod(img_pam_err2,2);
    err_cnt2 = 0;
    for i=1:(length(img_pam_demod2))
      if(img_pam_demod2(i) ~= image_num(i))
          err_cnt2 = err_cnt2 + 1;
      end
    end
    if(z == 6 || z == 11)
      %4
      demod_image_bin = num2str(img_pam_demod);
      for i=1:(N*N*RGB)
          demod_image_reshape(i,:)='00000000';
      end
      for i = 1:sizeofimage(1) 
        for j = 1:sizeofimage(2)/2
           demod_image_reshape(i,((j-1)*2+1:j*2)) = img_pam_demod(4*(i-1)+j,:);
        end
     end

      demod_image= bin2dec(demod_image_reshape);
      demod_image_arr = zeros(N,N,RGB,'uint8');
      for k = 1:RGB 
          for j = 1:N
              for i = 1:N
                  demod_image_arr(i,j,k)=demod_image((i+((j-1)*N))+((k-1)*N*N));
              end
          end
      end
    end
    %2
    demod_image_bin2 = num2str(img_pam_demod2);
    for i=1:(N*N*RGB)
        demod_image_reshape2(i,:)='00000000';
    end
    for j = 0:(log2(N))
        for r = 1:(N*N*RGB)
            demod_image_reshape2(r,(j+1))=demod_image_bin2(r+(j*(N*N*RGB)));
        end
    end
    demod_image2= bin2dec(demod_image_reshape2);
    demod_image_arr2 = zeros(N,N,RGB,'uint8');
    for k = 1:RGB 
        for j = 1:N
            for i = 1:N
                demod_image_arr2(i,j,k)=demod_image2((i+((j-1)*N))+((k-1)*N*N));
            end
        end
    end
    if(z == 6)
        figure
        subplot(3,1,1), imshow(image);
        title('SNR 5: Transmitted  Image 4PAM');
        subplot(3,1,3), imshow(demod_image_arr);
        title('SNR 5: Received Image 4PAM');
        subplot(3,1,2), imshow(demod_image_arr2);
        title('SNR 5: Received Image 2PAM');
        %subplot(3,1,1), stem(img_pam); ylim([-2 2]); title('Modulated Signal at 5 SNR');
        %subplot(3,1,2), plot(noisy); title('Noisy Signal at 5 SNR');
        %subplot(3,1,3), stem(img_pam_demod); title('Demodulated Signal at 5 SNR');
    end
    if(z == 11)
        figure
        subplot(3,1,1), imshow(image);
        title('SNR 10: Transmitted  Image 4PAM');
        subplot(3,1,3), imshow(demod_image_arr);
        title('SNR 10: Received Image 4PAM');
        subplot(3,1,2), imshow(demod_image_arr2);
        title('SNR 10: Received Image 2PAM');
        %figure
        %subplot(3,1,1), stem(img_pam);  ylim([-2 2]); title('Modulated Signal at 10 SNR');
        %subplot(3,1,2), plot(noisy);  title('Noisy Signal at 10 SNR');
        %subplot(3,1,3), stem(img_pam_demod); title('Demodulated Signal at 10 SNR');
    end
    err_cnt;
    err_cnt = err_cnt/(128*128*3*4);
    err_arr = [err_arr err_cnt]
    
    err_cnt2;
    err_cnt2 = err_cnt2/(128*128*3*8)
    err_arr2 = [err_arr2 err_cnt2]
end

qfun = zeros(1,11);
for l = 1:11
    snr = l-1;
    eb=No * (10^(snr/10));
    es = eb*log2(M);
    d=sqrt((es*3)/((M^2)-1));
      qfun(l) =2*(1-1/M)*qfunc(d/(sqrt(No/2)));
end
for l = 1:11
    snr = l-1;
    No = k*T;
    eb= No*(10^(snr/10));
    qfun2(l) = qfunc(sqrt(2*eb/(No)));
end
snr = 0:10;
figure
plot(snr,10*log10(err_arr));
hold on 
plot(snr,10 * log10(qfun),'*k');
title('4PAM Estimated SER vs Real SER');
legend('4PAM Estimated SER','2*(1-1/M)*qfunc(d/(sqrt(No/2)))');
ylabel('BER'); xlabel('SNR')
hold off
figure
plot(snr,10*log10(err_arr2));
hold on 
plot(snr,10 * log10(qfun2),'*k');
title('2PAM Estimated BER vs Real BER');
legend('2PAM Estimated BER','Qfunc(sqrt(2Eb/N0))');
ylabel('BER'); xlabel('SNR')
hold off
figure
plot(snr,10*log10(err_arr/2),'k');
hold on 
plot(snr,10*log10(err_arr2),'--k');
plot(snr,10 * log10(qfun2),'*k');
title('2PAM Estimated and Calculated BER, and 4PAM Estimated BER');
legend('4PAM Estimated BER','2PAM Estimated BER','Qfunc(sqrt(2Eb/N0))');
ylabel('BER'); xlabel('SNR')
hold off
figure
plot(snr,abs(((10 * log10(qfun(1:11))) - 10*log10(err_arr(1:11)))./(10 * log10(qfun(1:11))))*100,'k');
title('Percent Error SER vs Real SER');
legend('Percent Error');
ylabel('Percent Error'); xlabel('SNR')
