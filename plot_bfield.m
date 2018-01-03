function plot_bfield(A,coil_coord,goodchannels)

% draw image showing strength of magnetic field
% code taken from nut_beamforming_gui.m on March 23, 2006
%
% A = (N x 1) mixing matrix
% coil_coord = optional (M x 3) matrix of coil coordinates
% goodchannels = optional (1 x N) vector of indices of coil_coord (1 <= integer <= M)

% From Kenny:
%Dear all,
%
%Here is a stand alone piece of code to plot the activation map for,
%e.g., a column of a learned mixing matrix. It takes up to three inputs
%(a single column of a mixing matrix, the (M x 3) matrix of sensor
%coordinates, and an optional input that represents the vector of indices
%of the good channels). This essentially is a stand alone piece of code
%for part of what is done in Nutmeg (nut_beamforming_gui.m).
%
%To play with it try out the following line (which uses a default set of
%275 coil coordinates and assumes that all 275 channels are good),
%
%plot_bfield(randn(275,1)) 
% 

  if nargin < 2
    coil_coord = [73.3388   12.6902  128.5524
                  62.0249   30.6806  132.7548
                  53.2659   50.4677  131.1133
                  39.8355   69.5157  126.8864
                  25.2830   84.0555  121.2299
                  7.4133   94.5140  113.9441
                  -12.2101   99.0546  101.9840
                  42.0507   33.1476  141.1546
                  33.0393   52.7714  138.8645
                  19.4620   69.0442  135.0790
                  2.4436   81.0860  131.3413
                  -16.5945   89.3050  121.1276
                  1.3744   61.1410  144.3035
                  -19.3370   75.8725  138.1702
                  15.1944   44.3772  147.9967
                  -19.2395   56.6911  149.7210
                  39.4166   11.9101  144.3720
                  21.7312   23.5430  149.9385
                  -5.4227   40.5364  154.1892
                  -28.4381   38.3261  156.5290
                  -47.5546   28.1121  155.6483
                  4.2479   11.7384  156.0542
                  -14.9377   21.4520  158.9587
                  -34.6700   10.7106  160.0417
                  119.3946   28.5455   32.4093
                  111.2351   50.1441   31.6962
                  98.6137   69.1023   28.1610
                  82.3627   83.2382   20.8422
                  122.5750   16.1550   51.8827
                  117.2440   37.4653   52.4747
                  108.9791   56.9955   51.9076
                  93.7216   74.8072   48.8261
                  76.8558   86.3516   41.9861
                  117.4395   26.3424   72.0658
                  111.0165   47.0795   71.3905
                  98.3192   67.1894   70.1010
                  82.9483   81.5154   66.0418
                  64.6988   91.2889   59.4111
                  111.9955   14.3458   90.5360
                  108.4788   35.5292   90.5323
                  97.9682   56.9294   89.8448
                  84.7030   73.5975   87.6157
                  67.5220   87.5694   80.1034
                  46.9817   97.4855   70.1813
                  100.4214   20.8975  107.6814
                  92.8445   44.0033  106.8891
                  81.1120   62.4176  105.5659
                  66.3927   78.0426  101.1841
                  48.9094   92.8344   90.9016
                  28.6873  102.9359   79.6742
                  82.8034   30.5764  120.4583
                  72.7499   49.4742  121.3237
                  59.8664   66.8386  118.4055
                  46.6527   82.4160  111.3004
                  30.9377   95.7074  101.5639
                  13.1274  102.8856   93.3444
                  -6.9885  105.3170   82.1931
                  -119.9072   25.0890   83.3041
                  -112.2254   46.6768   82.0073
                  -101.1607   66.6603   80.7719
                  -85.4185   83.1876   78.4578
                  -127.1347   15.1247   64.8642
                  -121.4257   36.5963   64.6921
                  -111.8662   57.8044   62.8612
                  -97.6483   76.0258   61.0303
                  -127.4641   26.0439   46.1641
                  -120.5174   47.8996   45.0496
                  -108.5114   67.7061   43.5051
                  -92.1929   84.1229   41.6507
                  -131.9419   15.4729   27.2984
                  -127.5510   37.8614   26.5719
                  -118.0350   58.7357   25.5904
                  -103.8126   77.0065   24.2267
                  -133.1197   27.0008    7.8298
                  -126.1457   49.2490    7.0564
                  -114.0513   69.0787    6.2902
                  -65.7999   17.9482  149.5652
                  -46.1156   49.3305  150.1910
                  -83.0539   10.8118  137.7799
                  -65.2481   40.0732  145.3159
                  -38.5655   67.6998  141.5255
                  -95.7831   20.7957  123.1261
                  -82.9247   36.6519  132.7228
                  -72.6572   57.3646  133.0034
                  -55.5358   74.1685  130.4352
                  -36.7950   83.5266  126.4124
                  -97.1391   43.2604  115.7762
                  -86.1183   60.7795  116.1416
                  -73.2138   77.7001  114.5110
                  -52.8561   88.5039  112.2972
                  -32.7777   95.1976  107.9360
                  -111.8897   11.6289  103.4182
                  -109.2569   35.0126   99.6253
                  -100.6040   54.8613   99.1551
                  -87.9151   73.2242   98.0987
                  -71.7124   87.9190   94.9937
                  -50.3421   97.7896   92.7511
                  -29.2779  102.1003   88.0445
                  47.5109   98.4321   48.5317
                  28.8695  104.3640   57.7313
                  9.4064  106.3880   67.1347
                  -22.3842  107.9319   67.6464
                  -44.6699  107.5827   72.6454
                  -66.9128   95.9259   75.4921
                  59.8219   94.5739   31.0489
                  30.8322  104.5656   36.0480
                  10.9937  110.5188   45.7824
                  -10.8927  114.3022   49.8003
                  -38.5511  115.1104   52.4909
                  -61.0961  108.8145   55.8546
                  -79.4344   90.3004   58.5638
                  64.9033   92.6535    9.9677
                  43.1784  101.0089   18.5740
                  13.4743  111.5416   23.9934
                  -8.5875  116.6058   28.3127
                  -30.3823  116.7406   32.1854
                  -53.0630  114.2365   35.9287
                  -74.6839  100.9973   39.3468
                  47.5639  100.2631   -2.5525
                  25.9315  107.1074    6.4157
                  -0.6750  117.3523    7.8225
                  -22.6747  117.8051   11.4793
                  -44.9966  117.3914   15.7563
                  -68.1220  110.4399   19.9037
                  -85.6887   91.4031   22.1519
                  30.0282  106.4541  -14.8230
                  8.4378  115.4055  -12.2947
                  -14.5190  118.8318   -8.8905
                  -37.3162  119.0778   -4.7033
                  -60.3725  115.7753   -0.1685
                  -80.7129  101.4488    2.9435
                  -97.9133   85.4934    4.9486
                  73.5514   -9.1036  127.6995
                  62.5693  -27.7319  130.6873
                  54.0926  -47.4964  127.7351
                  41.0890  -66.3542  122.1619
                  26.8075  -80.8102  115.4909
                  9.2188  -90.9604  107.4353
                  -10.4617  -95.0313   95.3967
                  42.6452  -30.9585  138.8736
                  34.0027  -50.6093  135.2862
                  20.7841  -66.9623  130.4116
                  3.8740  -78.8679  125.8230
                  -15.1429  -86.6854  115.0510
                  2.4612  -60.0799  140.0573
                  -18.6509  -74.4143  132.8523
                  15.9951  -43.2203  144.9770
                  -18.0876  -56.3861  145.9032
                  39.7223  -10.1536  143.6876
                  22.2551  -22.5325  148.4493
                  -4.5494  -40.2140  151.3928
                  -27.5919  -38.5793  153.9069
                  -46.9480  -28.5867  153.8707
                  4.4423  -11.3837  155.1025
                  -14.4666  -21.7109  157.5912
                  -34.4799  -11.4062  159.3669
                  119.7969  -17.2975   30.9365
                  111.9950  -38.8266   28.4437
                  99.9100  -57.8009   23.7628
                  83.7836  -71.7438   15.5254
                  122.8359   -6.2921   51.1632
                  117.8943  -27.7172   50.2317
                  109.9319  -47.1429   48.1102
                  95.0503  -65.1508   43.8942
                  78.3310  -76.5202   36.2783
                  117.9005  -17.9621   70.4684
                  111.8730  -38.4976   68.2476
                  99.7926  -58.8418   65.7598
                  84.5159  -73.3632   60.6866
                  66.3153  -82.8999   53.3567
                  112.2206   -7.5586   89.7322
                  109.1229  -28.7465   88.4119
                  99.0090  -50.0804   86.2031
                  85.9506  -66.7471   82.9768
                  69.0911  -80.7124   74.2346
                  48.5159  -90.2818   63.3958
                  100.8592  -15.5101  106.4276
                  93.5927  -38.5164  104.0492
                  82.2144  -56.9120  101.3663
                  67.7309  -72.5257   95.9533
                  50.5591  -87.1315   84.7231
                  30.2841  -96.8929   72.6600
                  83.2922  -26.1995  118.4965
                  73.5789  -45.3442  118.1118
                  60.8494  -62.9171  113.7715
                  48.0584  -77.8875  105.8709
                  32.5698  -91.0279   95.1800
                  15.3768  -97.6518   85.9892
                  -5.1612  -99.8751   75.2081
                  -119.3151  -22.2483   81.4916
                  -111.5288  -43.5248   79.0007
                  -100.1274  -62.9626   76.1785
                  -83.7874  -78.9021   72.8234
                  -126.7317  -11.3459   63.8750
                  -120.7585  -32.4341   62.3523
                  -110.8583  -53.3848   59.0418
                  -96.2952  -70.9378   56.0351
                  -126.9947  -20.7578   44.4355
                  -119.5318  -42.3513   41.9202
                  -107.2379  -61.6667   39.1763
                  -90.5945  -77.6261   35.9997
                  -131.5317   -9.0779   26.2683
                  -126.7636  -31.1925   24.1900
                  -116.9017  -51.6649   21.8912
                  -102.4235  -69.4736   19.1583
                  -132.6096  -19.2759    6.2865
                  -124.9424  -41.3116    3.8026
                  -112.5476  -60.6928    1.8375
                  -65.5132  -18.5250  148.2681
                  -45.1547  -49.5727  146.9287
                  -82.7611  -11.1674  136.9479
                  -64.5548  -40.1348  142.6007
                  -37.2159  -67.2626  136.8360
                  -95.3126  -20.1723  121.6084
                  -82.2995  -36.4511  130.1668
                  -71.6185  -56.8602  129.1520
                  -54.2299  -73.0447  125.2491
                  -35.3768  -81.7327  120.5743
                  -96.2379  -42.0480  113.2831
                  -84.9798  -59.3638  111.9302
                  -71.6680  -75.8772  109.1512
                  -51.2107  -86.0543  106.2527
                  -30.9860  -92.0981  101.6879
                  -111.8025  -10.0039  102.5052
                  -108.7480  -32.8862   97.4306
                  -99.6100  -52.5708   95.4955
                  -86.6496  -70.7339   93.0542
                  -70.0690  -84.7241   89.0922
                  -48.6442  -93.9278   86.1652
                  -27.5338  -97.3618   81.3323
                  49.2206  -89.6254   41.8315
                  30.6065  -96.6051   50.8373
                  11.2236  -99.5471   60.0679
                  -20.5400 -101.8621   60.4037
                  -42.8737 -102.2275   65.4236
                  -65.2205  -91.1964   69.0637
                  61.4273  -84.1162   24.8250
                  32.5830  -95.3183   29.0895
                  12.9310 -102.2913   38.3375
                  -9.0451 -106.6603   42.1579
                  -36.5776 -108.0560   44.7159
                  -59.3488 -102.5186   48.8364
                  -77.7625  -84.6777   52.4944
                  66.4846  -80.8702    3.9469
                  44.9194  -90.1262   11.9076
                  15.3447 -101.8047   16.7953
                  -6.6039 -107.4476   20.5331
                  -28.4063 -108.2211   24.5696
                  -51.0573 -106.4072   28.3147
                  -72.5838  -93.7906   32.6015
                  49.1645  -87.4590   -9.1664
                  27.6296  -95.6501   -0.6154
                  1.3295 -106.5172    0.1076
                  -20.7238 -107.6534    3.8193
                  -43.0596 -107.8582    7.9935
                  -66.1364 -101.6264   12.5254
                  -84.0566  -83.4410   16.0281
                  31.6992  -93.6042  -21.8515
                  10.4122 -102.5545  -19.9052
                  -12.5189 -107.4260  -16.7747
                  -35.2658 -108.0030  -12.5941
                  -58.3492 -105.4394   -7.8158
                  -78.7943  -91.8017   -3.7272
                  -96.0752  -76.3349   -0.6005
                  56.8737    1.3127  136.9478
                  22.0424    0.4084  150.8620
                  -15.1273   -0.1528  159.9833
                  -54.8757   -0.4392  155.4747
                  122.8440    5.6785   32.3230
                  120.2386    4.2016   71.5317
                  91.8490    2.4863  115.0426
                  -120.7390    1.4718   85.9523
                  -129.7026    2.5835   45.2531
                  -134.7757    3.7011    8.9002
                  -101.4804    0.4849  117.8778];
  end

  if nargin < 3
    goodchannels = 1:size(coil_coord,1);
  end

  if size(A,1) ~= length(goodchannels)
    disp(' ')
    disp(' ')
    disp('Size of mixing matrix does not match the number of good channels of coil_coord')
    disp('Change either A or goodchannels')
    error
  end

  % normalize coil_coord using all channels
  max_value1 = max(abs(coil_coord(:,1)));
  max_value2 = max(abs(coil_coord(:,2)));
  coil_coord(:,2) = coil_coord(:,2)*max_value1/max_value2; % make extent of first and second columns equal
  max_value2 = max(abs(coil_coord(:,2)));
  coil_coord(:,3) = coil_coord(:,3) - min(coil_coord(:,3)); % make min value 0 and max value 1
  coil_coord(:,3) = coil_coord(:,3)/max(coil_coord(:,3));

  % convert 3-d coil_coord into 2-d variable, out
  out(:,1) = coil_coord(:,1).*exp(-coil_coord(:,3));
  out(:,2) = coil_coord(:,2).*exp(-coil_coord(:,3));

  % normalize to uniformly fill a circle
  out(:,1) = out(:,1)*max_value1/max(abs(out(:,1)));
  out(:,2) = out(:,2)*max_value2/max(abs(out(:,2)));

  % adjust for size of border
  border = 4;
  resolution = 100;
  out = out*resolution/(resolution + border);

  % initialize
  L = resolution + 2*border;
  radius = 0.95*(L+1)/2; % radius of circle that represents head
  image_matrix = zeros(L); % create empty image

  % colormap
  % first row should be [0 0 0] (black) and second row should be [1 1 1], white
  c = zeros(4096+2,3); c(2,:) = ones(1,3);
  c((1:2048)+2,1) = (1:-1/2047:0)'; c((2049:4096)+2,3) = (0:1/2047:1)';
  cmap_length = size(c,1);

  % row indices where sensors lie (in 2-d image space)
  row = round((resolution-1)*(out(:,1)+max_value1)/(2*max_value1)) + 1 + border;
  col = round((resolution-1)*(out(:,2)+max_value2)/(2*max_value2)) + 1 + border;

  % ensure sensors lie in circle (that represents the head) - fix needed for BTI sensors
  sensor_radii = sqrt(((L+1)/2-row).^2 + ((L+1)/2-col).^2);
  if max(sensor_radii) > radius
    out = out*0.95*radius/max(sensor_radii);
    row = round((resolution-1)*(out(:,1)+max_value1)/(2*max_value1)) + 1 + border;
    col = round((resolution-1)*(out(:,2)+max_value2)/(2*max_value2)) + 1 + border;
  end

  % define indices
  sensor_ndx_good = (col(goodchannels)-1).*(size(image_matrix,1)) + row(goodchannels);
  sensor_ndx = (col-1).*(size(image_matrix,1)) + row;

  if size(A,1) >= 2
    Arange = max(max(A)) - min(min(A));
    maxvalue = max(max(A)) - Arange*0.1;
    minvalue = min(min(A)) + Arange*0.1;
    
    % create image (1 is black, 2 is white, the rest are defined in the colormap above)
    meandata = A(:,1); % ignore other columns
    f = find(meandata < minvalue); meandata(f) = minvalue; % introduce user-specified clipping
    f = find(meandata > maxvalue); meandata(f) = maxvalue;
    f = find(meandata < 0); f2 = find(meandata >= 0);
    meandata(f) = meandata(f)/(2*abs(minvalue)); % -0.5 to 0 are for negative A values
    meandata(f2) = meandata(f2)/(2*abs(maxvalue)); % 0 to 0.5 are for positive A values
    image_matrix(sensor_ndx_good) = meandata; % output takes values between -0.5 and 0.5
    image_matrix = fliplr(flipud(image_matrix));

    % smooth image
    maxnum = max(max(image_matrix));
    minnum = min(min(image_matrix));
    b = exp(-[-2:0.1:2].^2)'*exp(-[-2:0.1:2].^2);
    image_matrix = filter2(b,image_matrix);
    f = find(image_matrix < 0); f2 = find(image_matrix > 0);
    minnum_new = min(min(image_matrix));
    if minnum_new == 0, minnum_new = eps; end
    image_matrix(f) = image_matrix(f)*abs(minnum/minnum_new); % restore previous minimum
    image_matrix(f2) = image_matrix(f2)*maxnum/max(max(image_matrix)); % restore previous maximum
  end

  % compute x and y indices of image
  xstep = 2*max_value1/(size(image_matrix,1) - 2*border - 1);
  ystep = 2*max_value2/(size(image_matrix,2) - 2*border - 1);
  xndx = (-max_value1 - 2*xstep):xstep:(max_value1 + 2*xstep);
  yndx = (-max_value2 - 2*ystep):ystep:(max_value2 + 2*ystep);
  
  % adjust values of amplitude map
  image_matrix = image_matrix + 0.5; % output takes values between 0 and 1

  % adjust values for size of colormap
  image_matrix = round((cmap_length-3)*image_matrix) + 3; % avoid first 2 rows!

  % find indices of region lying outside of circle
  f2 = [];
  for i = 1:L
    f = find(((1:L)-(L+1)/2).^2 + (i-(L+1)/2)^2 > radius^2);
    f2 = [f2 f+L*(i-1)];
  end

  % remove portion outside of circle
  image_matrix(f2) = 2; % makes area outside of circle white

  lh = findobj('name','Plot B Field');
  if isempty(lh)
    lh = figure;
    set(lh,'name','Plot B Field')
    pos = get(lh,'position');
    set(lh,'position',[pos(1:2) 460 420])
  end
  set(lh,'Colormap',c);
  figure(lh); clf

  bf_image = image(xndx,yndx,image_matrix);
  set(gca,'XTick',[],'XTickLabel','','YTick',[],'YTickLabel','');
  a = axis; text(a(2)+0.04*(a(2)-a(1)),a(3)+0.44*(a(4)-a(3)),' Right','Rotation',270);
  title('Anterior'); ylabel('  Left'); xlabel('Posterior'); % title and axis labels must occur after text command above

  hold on
  plot(-out(goodchannels,2),-out(goodchannels,1),'ow');
  badchannels = 1:size(coil_coord,1); badchannels(goodchannels) = [];
  plot(-out(badchannels,2),-out(badchannels,1),'o','Color',[0.2 0.2 0.2]);

  return
