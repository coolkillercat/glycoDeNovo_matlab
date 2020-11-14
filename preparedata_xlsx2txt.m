datapath = '/Users/hong/Documents/Projects/Glycomics/data';
dataset = '20190116';
originalfiles = dir( [datapath, filesep, dataset, filesep, 'original', filesep, '*.xlsx'] );

for k = 1 : length(originalfiles)
    xlsfile = [originalfiles(k).folder, filesep, originalfiles(k).name];
    [~, sheets] = xlsfinfo( xlsfile );
    if originalfiles(k).name(1) == '~'
        continue;
    end
    for m = 1 : length(sheets)
        [num, txt, raw] = xlsread( xlsfile );
        sheetname = strrep( strrep( sheets{m}, 'Isomer ', '' ), ' ', '');
        resultfile = [datapath, filesep, dataset, filesep, strrep(originalfiles(k).name, 'xlsx', ''), sheetname, '.txt'];
        rfid = fopen( resultfile, 'w' );
        for r = 1 : size( txt, 1 )
            temp = txt{r, 1};
            if  contains( temp, '?' ) 
                temp = strrep( strrep( strrep( temp, '?', '-' ), '?', 'a' ), '?', 'b' );
                idxes = strfind( temp, ')' );
                for c = length(idxes) : -1 : 1
                    if temp(idxes(c)+1) ~= ']'
                        temp = [temp(1:idxes(c)), ' ', temp(idxes(c)+1:end)];
                    end
                end
                idxes = strfind( temp, ']' );
                for c = length(idxes) : -1 : 1
                    temp = [temp(1:idxes(c)), ' ', temp(idxes(c)+1:end)];
                end
                fprintf( rfid, '# %s\n', temp );
            else
                expidx = strfind(txt{r, 1}, '(exptl.)');
                if ~isempty(expidx)
                    fprintf( rfid, '%s\n\n', txt{r, 1}(1:expidx(1)-1) );
                elseif isempty(txt{r, 1})
                    continue;
                elseif strcmp( 'm/z', txt{r, 1} ) == 1
                    fprintf( rfid, '%s\t%s\t%s\t%s\n', txt{r, 1}, txt{r, 2}, txt{r, 3}, txt{r, 4} );
                    break;
                else
                    fprintf( rfid, '%s\n', txt{r, 1} );
                end
            end
        end
        for datar = 1 : size(num, 1)
            fprintf( rfid, '%f\t%s\t%d\t%.1f\n', num(datar,1), txt{r+datar, 2}, num(datar,3), num(datar,4) );
        end
        fclose( rfid );
        disp( resultfile );
    end
end