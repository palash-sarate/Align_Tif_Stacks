function parsed = temp(str)
    tokens = regexp(str,"\['([^']+)',\s*([0-9\.]+)\]", 'tokens', 'once');
    confidence = str2double(tokens{2});
    % confidence = confidence(1);    
    parsed = {tokens{1}, confidence};
    disp(parsed);
end

for i = 1:1%numel(stack_info.ocr_results)
    len = numel(stack_info.ocr_results{i});
    for j = 1:len
        tempdata = stack_info.ocr_results{i}{j};
        parsed = temp(tempdata);
        % disp([tempdata, parsed])
    end
end
