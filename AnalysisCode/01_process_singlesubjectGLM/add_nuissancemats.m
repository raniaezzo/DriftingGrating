function matrices = add_nuissancemats(bidsDir,sub, ses, matrices);

    fmriprep_ses = fullfile(bidsDir, 'derivatives', 'fmriprep', sub, ses, 'func');

    sessionConfoundfiles = dir(fullfile(fmriprep_ses, '*confounds*.tsv'));

    for fi=1:numel(matrices)

        disp('PLEASE VERIFY THAT RUN MATRIX AND CONFOUND FILE RUN NUMBER MATCH')
        disp(sprintf('Run number: %i', fi))

        confoundFilename = sessionConfoundfiles(fi).name;
        disp(sprintf('confoundFilename: %i', fi))
    
        % Specify the file path
        file_path = fullfile(fmriprep_ses, confoundFilename);

        copyfile(file_path, strrep(file_path, '.tsv', '.txt'));
    
        % Read the TSV file into a table
        %data_table = readtable(file_path, 'Delimiter', '\t');
        data_table = readtable(strrep(file_path, '.tsv', '.txt'), 'Delimiter', '\t');
        
        % Specify the header names you want to extract
        desired_headers = {'framewise_displacement', 'csf', 'white_matter', 'global_signal', 'trans_x', 'trans_y', 'trans_z', 'rot_x', 'rot_y', 'rot_z'};
        
        % Extract columns based on the specified header names
        selected_columns = data_table(:, desired_headers);

        % change any nans to 0
        selected_columns.framewise_displacement(isnan(selected_columns.framewise_displacement)) = 0;

        matrices{fi} = [matrices{fi} table2array(selected_columns)];

        % delete the text file that was created from the tsv file
        delete(strrep(file_path, '.tsv', '.txt'))

    end


end
