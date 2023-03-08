function vpixxShutdown(const)
    if const.vpixx == 1
        % datapixx shutdown
        Datapixx('RegWrRd');
        Datapixx('StopAllSchedules');
        Datapixx('Close');
    end
end