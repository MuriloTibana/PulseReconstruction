classdef Laser
    properties
        frequency {mustBeNonnegative}
        FWHM {mustBeNonnegative}
        cycles {mustBeNonnegative}
        period {mustBeNonnegative}
        chirp {mustBeReal}
        amplitude {mustBeNonnegative}
        position {mustBeNumeric}
        phase {mustBeNumeric}
    end
    methods
        function obj = Laser(amplitude,frequency,FWHM,chirp,position,phase)
            % Constructor for a laser object
            %   Can pass either the amplitude, frequency, FWHM, and
            %   position or the parameters list that the prep() method
            %   generates. The phase can be included in a cartesian or
            %   polar way with the optional phase parameter.
            if nargin == 1
                obj.frequency = amplitude(1);
                obj.period = amplitude(2);
                obj.chirp = amplitude(6);
                obj.amplitude = abs(amplitude(3) + 1i*amplitude(4)); 
                obj.position = amplitude(5);
                obj.phase = angle(amplitude(3) + 1i*amplitude(4));
                obj.cycles = obj.period * obj.frequency / (2 * pi);
                obj.FWHM = obj.period * (2*sqrt(2*log(2)));
            elseif nargin == 5
                obj.amplitude = abs(amplitude);
                obj.frequency = frequency;
                obj.FWHM = FWHM;
                obj.period = FWHM / (2*sqrt(2*log(2)));
                obj.chirp = chirp;
                obj.cycles = obj.period * obj.frequency / (2 * pi);
                obj.position = position;
                obj.phase = angle(amplitude);
            elseif nargin == 6
                obj.amplitude = abs(amplitude);
                obj.phase = phase;
                obj.frequency = frequency;
                obj.FWHM = FWHM;
                obj.period = FWHM / (2*sqrt(2*log(2)));
                obj.chirp = chirp;
                obj.cycles = obj.period * obj.frequency / (2 * pi);
                obj.position = position;
            end
        end
        function prep = params(obj,single_color,chirp)
            % Conversion from laser object to list of parameters for the
            % matrix element calculations. This will also convert a list of
            % laser objects to a list of laser parameters.
            if nargin == 1
                single_color = false;
                chirp = true;
            end
            if numel(obj) == 1
                complex_amplitude = obj.amplitude .* exp(1i.*obj.phase);
                if(~single_color)
                    % Unknown frequencies
                    if(chirp)
                        prep = [obj.frequency obj.period real(complex_amplitude) imag(complex_amplitude) obj.position obj.chirp];
                    else
                        prep = [obj.frequency obj.period real(complex_amplitude) imag(complex_amplitude) obj.position];
                    end
                else
                    % Known/specified frequencies
                    if(chirp)
                        prep = [obj.period real(complex_amplitude) imag(complex_amplitude) obj.position obj.chirp];
                    else
                        prep = [obj.period real(complex_amplitude) imag(complex_amplitude) obj.position];
                    end
                end
            else
                prep = cell2mat(arrayfun(@(o) o.params(single_color,chirp), obj, 'UniformOutput', false));
            end
        end
        function value = calculate(obj,time)
            % Calculate the laser in time and return the value with
            % corresponding time.
            if numel(obj) == 1

                complex_period = 1/sqrt(1/obj.period^2 + 1i*obj.chirp);
                value = obj.amplitude .* exp(1i*obj.phase) * exp(-((time - obj.position).^2/...
                    (2)) * (1/complex_period^2) - 1i * obj.frequency .*...
                    time);
            else
                vals = arrayfun(@(o) o.calculate(time),obj,'UniformOutput',false);
                value = sum(cell2mat(vals),1);
            end
        end
        function plot(obj,time)
            % Plots a laser object in time or a list of laser objects in
            % time.
            vals = calculate(obj,time);
            figure; plot(time,real(vals),time,abs(vals)); 
            xlabel('time ')
            ylabel('real electric field')
            grid on
        end
    end
    methods(Static)
        function laser = generate(params,chirp,single_color_omega)
            laser = [];
            if(chirp)
                for i=1:size(params,1)
                    if size(params,2) == 4
                        new_laser = [single_color_omega(i) params(i,1:end-1) 0 params(i,end)];
                    elseif size(params,2) == 5
                        new_laser = [single_color_omega(i) params(i,:)];
                    else
                        new_laser = params(i,:);
                    end
                    laser = [laser; Laser(new_laser)];
                end
            else
                for i=1:size(params,1)
                    if size(params,2) == 3
                        new_laser = [single_color_omega(i) params(i,:) 0 0];
                    elseif size(params,2) == 4
                        new_laser = [single_color_omega(i) params(i,:) 0];
                    else
                        new_laser = [params(i,:) 0];
                    end
                    laser = [laser; Laser(new_laser)];
                end
            end
        end

        function wavelength = au2SI_wavelength(omega)
            % Returns wavelength in nm
            energy = omega .* 27.211;
            wavelength = 1239.8 ./ energy;
        end
        function intensity = au2SI_intensity(amplitude)
            % Returns intensity in W/cm^2
            intensity = amplitude.^2 .* 3.51e16;
        end
        function omega = SI2au_wavelength(wavelength)
            % Returns angular frequency in au
            energy = 1239.8 ./ wavelength;
            omega = energy ./ 27.211;
        end
        function amplitude = SI2au_intensity(intensity)
            % Returns field amplitude in au
            amplitude = sqrt(intensity ./ 3.51e16);
        end
        function duration_SI = SI2au_duration(duration_au)
            % Return time duration in fs
            duration_SI = duration_au / 41.32;
        end
        function duration_au = au2SI_duration(duration_SI)
            % Return time duration in fs
            duration_au = duration_SI * 41.32;
        end
    end
end
