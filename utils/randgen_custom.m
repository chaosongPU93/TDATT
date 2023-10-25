function [meanrand,pdfmisfit,exppdf] = randgen_custom(seed,N,data,wt,pdfxloc,pdfval)
% [randnum,muHat,sigmaHat,muCI,sigmaCI] = RANDGEN_CUSTOM(data, wt, dimsize)
%
% Given a data set 'data' with non-equal weights 'wt' that has a custom PDF,
% we wish to draw random numbers with the same dimension size of the summed
% weights, by N times (realizations), from its underlying PDF. This function 
% summarizes 3 methods that could deal with the univariate case (1D).
% 
% --1.Take the underlying PDF of data as custom and call a 3rd-party function 
%     'notmymatf/randpdf'. This option is equivalent to call 'randgen_pdf'.
% --2.Use kernel smoothing function estimate for univariate and bivariate 
%     data, call a built-in function 'ksdensity'. This option is equivalent
%     to call 'randgen_ksdensity'.
% --3.Assume the underlying PDF of data is Gaussian (normal)-like in 
%     principle, call my function 'utils/randn_wt' which uses 'normfit' 
%     and 'normrnd'. This option is equivalent to call 'randgen_norm'.
%     
%
% It takes in random seed state for reproduction.
% Then it evaluate the how well the PDF of generated random numbers matches
% that of data.
% Output the median of randoms and the misfit between randoms PDF and data
% PDF in each round of generation. 
%
% By Chao Song, chaosong@princeton.edu
% First created date:   2020/10/01
% Last modified date:   2022/04/06

meanrand = zeros(N, 3);    % mean of the generated random numbers
pdfmisfit = zeros(N, 3);    % misfits of PDFs between data and generated random numbers
dimsize = [ceil(sum(wt)),1];  % size of desired randum numbers per sampling
for i = 1: N
% i = randi(1000,1);

    rdnum = [];
    
    % 1. use randpdf
    %%% NOTE:
    %%% the shifting arises bc. the obtained pdf requires certain binning, so if you want to shift
    %%% the data by amount 'shift' by directly add it to the data, it is fine; but if the bin width
    %%% and edge remains unchanged, the histogram (pdf) usually would NOT have a same shape as the
    %%% unshifted histogram (pdf), so preserve the shape, you MUST add the shift to the pdf x
    %%% location, rather than do binning again to the shifted new data
    rng(seed{i});
    [randnum, exppdf] = randpdf(pdfxloc, pdfval, dimsize);
    rdnum = [rdnum randnum];

    % 2. use ksdensity
    rng(seed{i});
    randnum = ksdensity(data, rand(dimsize), 'function', 'icdf', 'weight', wt);
    rdnum = [rdnum randnum];
    
    % 3. assume pdf is normal-like, and try to fit it with a gaussian, call 
    rng(seed{i});
    [randnum,muHat,sigmaHat,~,~] = randn_wt(data, wt, dimsize);
    rdnum = [rdnum randnum];
%     some sense of the misfit
%     [pdfmisfit1, pdfmisfit2, f66] = compare_distribution2(randnum,muHat,sigmaHat,vdist,wt,pdfxloc,pdfval);
    
%     %%% check the distribution of generated random numbers to see if it accords with the pdf of data
%     pdfmisf = zeros(size( rdnum,2), 1);
%     titlestr = {'rand gen from custom pdf';
%                 'rand gen from ksdensity';
%                 'rand gen from normal fitting';
%            };
%     for ii = 1: size(rdnum,2)
%         binw = pdfxloc(2)-pdfxloc(1);
%         [pdfmisf(ii), f55] = compare_distribution1(rdnum(:,ii),vdist,wt,pdfxloc,pdfval,exppdf,binw);
%         supertit(f55.ax, string(titlestr{ii}),12);
%         
%     end
    
    meanrand(i,:) = mean(rdnum,1);
    
    %%% output the misfits of PDFs as well but no figures
    for ii = 1: size(rdnum,2)
        binw = pdfxloc(2)-pdfxloc(1);
        pdfmisfit(i, ii) = check_randgen(rdnum(:,ii),pdfxloc,pdfval,exppdf,binw);
    end
    
end

% keyboard









