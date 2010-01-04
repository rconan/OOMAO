function out = deviceDriver(varargin)

driverData = varargin;

    function status = openDevice
        status = true;
    end

    function status = closeDevice
        status = true;
    end

    function out = runDevice
        out = rand(32);        
    end

out.openDevice  = @openDevice;
out.closeDevice = @closeDevice;
out.runDevice   = @runDevice;

end