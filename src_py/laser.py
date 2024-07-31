import numpy as np

class Laser:
    def __init__(self, amplitude,frequency,FWHM,chirp,position,phase):
        self.amplitude = np.abs(amplitude)
        self.frequency = frequency
        self.FWHM = FWHM
        self.period = FWHM / (2 * np.sqrt(2 * np.log(2)))
        self.chirp = chirp
        self.cycles = self.period * self.frequency / (2 * np.pi)
        self.position = position
        self.phase = phase

    def params(self, single_color=False, chirp=True):
        complex_amplitude = self.amplitude * np.exp(1j * self.phase)
        if single_color:
            if chirp:
                return [self.period, np.real(complex_amplitude), np.imag(complex_amplitude), self.position, self.chirp]
            else:
                return [self.period, np.real(complex_amplitude), np.imag(complex_amplitude), self.position]
        else:
            if chirp:
                return [self.frequency, self.period, np.real(complex_amplitude), np.imag(complex_amplitude), self.position, self.chirp]
            else:
                return [self.frequency, self.period, np.real(complex_amplitude), np.imag(complex_amplitude), self.position]

    def calculate(self, time):
        if isinstance(self, Laser):  
            complex_period = 1 / np.sqrt(1 / self.period ** 2 + 1j * self.chirp)
            value = self.amplitude * np.exp(1j * self.phase) * np.exp(-((time - self.position) ** 2 / (2))* (1 / complex_period ** 2) - 1j * self.frequency * time)
        else:  # Handling list of Laser instances
            print('More than one value')
            values = [obj.calculate(time) for obj in self]
            value = np.sum(values, axis=0)
        return value
    
    def plot(self, time):
        import matplotlib.pyplot as plt
        vals = self.calculate(time)
        plt.figure()
        plt.plot(time, np.real(vals), label='Real Part')
        plt.plot(time, np.abs(vals), label='Magnitude')
        plt.xlabel('Time')
        plt.ylabel('Electric Field')
        plt.grid(True)
        plt.legend()
        plt.show()

    @staticmethod
    def generate(params, chirp=True, single_color_omega=None):
        lasers = []
        if chirp:
            for i in range(len(params)):
                if len(params[i]) == 4:
                    new_laser = [single_color_omega[i], *params[i][:-1], 0, params[i][-1]]
                elif len(params[i]) == 5:
                    new_laser = [single_color_omega[i], *params[i]]
                else:
                    new_laser = params[i]
                lasers.append(Laser(new_laser))
        else:
            for i in range(len(params)):
                if len(params[i]) == 3:
                    new_laser = [single_color_omega[i], *params[i], 0, 0]
                elif len(params[i]) == 4:
                    new_laser = [single_color_omega[i], *params[i], 0]
                else:
                    new_laser = [*params[i], 0]
                lasers.append(Laser(new_laser))
        return lasers

    @staticmethod
    def au2SI_wavelength(omega):
        energy = omega * 27.211
        wavelength = 1239.8 / energy
        return wavelength

    @staticmethod
    def au2SI_intensity(amplitude):
        intensity = amplitude ** 2 * 3.51e16
        return intensity

    @staticmethod
    def SI2au_wavelength(wavelength):
        energy = 1239.8 / wavelength
        omega = energy / 27.211
        return omega

    @staticmethod
    def SI2au_intensity(intensity):
        amplitude = np.sqrt(intensity / 3.51e16)
        return amplitude

    @staticmethod
    def SI2au_duration(duration_SI):
        duration_au = duration_SI / 41.32
        return duration_au

    @staticmethod
    def au2SI_duration(duration_au):
        duration_SI = duration_au * 41.32
        return duration_SI
