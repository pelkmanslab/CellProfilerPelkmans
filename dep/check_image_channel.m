function [intChannelNumber,intZstackNumber,intActionNumber] = check_image_channel(strImageName)

    if nargin==0
%         strImageName = 'DAPI_p53SLS - n000000.tif';
%         strImageName = '040423_frank_si_VSV_B02_25_w460.tif';
%         strImageName = '2008-08-14_HPV16_batch1_CP001-1ea_A20_9_w530.tif';

%         strImageName = '081029_Eva_FAKKO_20X_ABCA1inh_SV40_A05f00d0.png';

%         strImageName = 'Dapi - n000000.tif';
%         strImageName = 'Y:\Data\Users\Jean-Philippe\p53hitvalidation\ko_p53Pro72_plate1_triplicate1\TIFF\Well H12\Dapi_p53SLS - n000000.tif';
%         strImageName = '080611-olivia-1_A10_s1_w12E22EFEB-B167-43E0-A05F-997CCA19728A.tif'
%         strImageName = '080611-olivia-1_A10_s12E22EFEB-B167-43E0-A05F-997CCA19728A.tif'
%         strImageName = 'CellNumber_C04_s10B80DBC85-9B51-42C6-909A-72025A001DCA.png'
%         strImageName = 'ExpressionPlate2_B02_s2_w1065BDCCC-B2EB-432C-9E44-A7F1E0D9E154';
%          strImageName = '120306_Tf4881_NIKON_t0001_A02_s010_w02.png';        
         strImageName = 'blablabl_NIKON_t0001_C21_s1_z0001_w01_ledGFP.png';
    end

    
    intChannelNumber = NaN;
    intZstackNumber = NaN;
    intActionNumber = NaN;
    
    % CW
    strNomenclature1 = regexp(strImageName,'f\d\dd\d.(tif|png)','Match');
    strNomenclature1a = regexp(strImageName,'f\dd\d.(tif|png)','Match');
    strNomenclature2 = regexp(strImageName,'_w\d\d\d.(tif|png)','Match');
    strNomenclature3 = regexp(strImageName,' - n\d{2,}.(tif|png)','Match');
    
    % NIKON
    strNomenclaturePre4 = regexp(strImageName,'NIKON.*_t\d{1,}(_z\d{1,}_|_)[A-Z]\d{2,3}_s\d{1,}_w\d{1,}[^\.]*\.(tiff?|png)','Match');

    % MD MICROEXPRESS
    strNomenclature4 = regexp(strImageName,'_\w\d\d_s\d{1,}_w\d','Match');
    strNomenclature4a = regexp(strImageName,'_\w\d\d_s\d{1,}[A-Z0-9\-]{36}\.(tif|png)','Match');
    
    % CV7K
    strNomenclature5 = regexp(strImageName, ...
        '_([^_]{3})_(T\d{4})(F\d{3})(L\d{2})(A\d{2})(Z\d{2})(C\d{2})', 'Match');
    
    strChannelWavelengths = {'_w460','_w530','_w595','_w685'};
    strChannelDescriptions = {'DAPI','FITC','TRITC','CY5'};

    if not(isempty(strNomenclature1))
        %%% CELLWORX
        strChannelMatch = regexp(strImageName,'\d.(tif|png)','Match');
        intChannelNumber = str2double(strrep(strrep(char(strChannelMatch{1}),'.tif',''),'.png',''))+1;
    elseif not(isempty(strNomenclature2))
        %%% CELLWORX        
        strChannelMatch = regexp(strImageName,'_w\d\d\d.(tif|png)','Match');
        intChannelNumber = find(~cellfun('isempty',strfind(strChannelWavelengths,char(strrep(strrep(strChannelMatch,'.tif',''),'.png','')))));
   elseif not(isempty(strNomenclature1a))
        %%% CELLWORX       
        strChannelMatch = regexp(strImageName,'\d.(tif|png)','Match');
        intChannelNumber = str2double(strrep(strrep(char(strChannelMatch{1}),'.tif',''),'.png',''))+1;
   elseif not(isempty(strNomenclature3))
        %%% BD-PATHWAY
        for i = 1:length(strChannelDescriptions)
            if ~isempty(findstr(upper(strChannelDescriptions{i}),upper(strImageName)))
                intChannelNumber = i;
            end
        end
    elseif not(isempty(strNomenclaturePre4))
        %%% NIKON
        %strChannelMatch = regexp(strImageName,'NIKON.*_\w\d\d_s\d{3,}_w(\d\d)','Tokens');
        strChannelMatch = regexpi(strImageName,'NIKON.*_t\d{1,}(_z\d{1,}_|_)[A-Z]\d{2,3}_s(\d{1,})_w(\d{1,})[^\.]*\.(tiff?|png)','Tokens');
        intChannelNumber = str2double(strChannelMatch{1}(3));
        strZstack = char(strChannelMatch{1}(1));
        if ~strcmp(strZstack,'_') % does filename actually contain a z-stack info?
            intZstackNumber = str2double(strZstack(3:end-1));
        else
            intZstackNumber = NaN;
        end
    elseif not(isempty(strNomenclature4))
        %%% MD
        strChannelMatch = regexp(strImageName,'_\w\d\d_s\d{1,}_w(\d)','Tokens');
        intChannelNumber = str2double(strChannelMatch{1});
    elseif not(isempty(strNomenclature4a))
        %%% MD - if there is only one channel present
        intChannelNumber = 1;
    elseif not(isempty(strNomenclature5))    
        %%% CV7K - we have a match against "Yokogawa" filenaming tail 
        strChannelMatch = regexp(strImageName, '_([^_]{3})_(T\d{4})F(\d{3})L(\d{2})A(\d{2})Z(\d{2})C(\d{2})', 'Tokens');
        intChannelNumber = str2double(strChannelMatch{1}(7));
        intZstackNumber = str2double(strChannelMatch{1}(6));
        intActionNumber = str2double(strChannelMatch{1}(5));
    else
        warning('iBRAIN:check_image_channel','unknown file name %s',strImageName)
    end
end
