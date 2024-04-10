from math import log, exp

R = 8.314     # Gas constant (J/(mol K))
F = 96.485    # Faraday's constant (C/mol)
ATM = 101325  # 1 Atmosphere (Pa)

class PEMFC_1D:
	def __init__(self,  \
		E_thermo = 1.0,  j = 0.5,  T = 343,  p_SAT = 0.307,  \
		x_O2 = 0.19,  x_H2O = 0.1,  pC = 3,  pA = 3,  \
		Deff_H2 = 0.149,  Deff_O2 = 0.0295,  D_lambda = 3.81E-6,  alpha = 0.5,  \
		j_0 = 0.0001,  tM = 125,  tA = 350,  tC = 350
	):
		self.E_thermo = E_thermo  # Thermodynamic voltage (V)
		self.j = j                # Operating current density (A/cm^2)
		self.T = T                # Temperature (K)
		self.p_SAT = p_SAT        # Vapor saturation pressure (atm)
		self.x_O2 = x_O2          # Oxygen mole fraction (unitless)
		self.x_H2O = x_H2O        # Cathode water mole fraction (unitless)
		self.pC = pC              # Cathode pressure (atm)
		self.pA = pA              # Anode pressure (atm)
		self.Deff_H2 = Deff_H2    # Effective hydrogen (or water) diffusivity (cm^2/s)
		self.Deff_O2 = Deff_O2    # Effective oxygen (or water) diffusivity (cm^2/s)
		self.D_lambda = D_lambda  # Water diffusivity in Nafion (cm^2/s)
		self.alpha = alpha        # Transfer coefficient
		self.j_0 = j_0            # Exchange current density (A/cm^2)
		self.tM = tM              # Electrolyte thickness (nano-m)
		self.tA = tA              # Anode thickness (nano-m)
		self.tC = tC              # Cathode thickness (nano-m)

		# (My) Calculate constants alpha_star and C needed for R_m
		_tM_cm = self.tM * 1e-6  # convert to cm
		_tA_m = self.tA * 1e-9   # convert to m
		_tC_m = self.tC * 1e-9
		_pA_Pa = self.pA * ATM  # convert to Pa
		_pC_Pa = self.pC * ATM
		_j_m2 = self.j * 1e-4             # convert to A/m^2
		_Deff_O2_m2 = self.Deff_O2 * 1e4  # convert to m^2/s
		_Deff_H2_m2 = self.Deff_H2 * 1e4

		_E = exp(0.000598 * self.j * _tM_cm / self.D_lambda)  # Cathode constants
		print("_E", _E)
		_KC = 1.4 * self.pC / self.p_SAT
		print("_KC", _KC)
		_DC = 2 * F * _pC_Pa * _Deff_O2_m2
		print("_DC =", _DC)
		_KA = 14 * self.pA / self.p_SAT                        # Anode constants
		print("_KA", _KA)
		_DA = 2 * F * _pA_Pa * _Deff_H2_m2
		print("_DA", _DA)

		_AC = 4.4 + _KC * _tC_m * _j_m2 * R * self.T / _DC  # Linear constants: A*alpha_star + E?*C = B
		print("_AC =", _AC)
		_AA = 4.4 + _KA * _tA_m / _DA * _j_m2 * R * self.T
		print("_AA", _AA)
		_BC = 12.6 + _KC * self.x_H2O - _KC * _tC_m * _j_m2 * R * self.T / _DC
		print("_BC =", _BC)
		_BA = _KA * self.x_H2O
		print("_BA", _BA)

		self._alpha_star = (_E * _BA - _BC) / (_E * _AA - _AC)
		self._C = _BA - _AA * self._alpha_star

	def R_m(self):
		""" (Area-specific) resistance of the membrane. (eq. 4.18, p. 120) (eq. 6.52, p. 214) """
		return 0.117  # TODO: where did this number come from??
	
	def eta_ohmic(self):
		""" Ohmic overvoltage. (eq. 6.53, p. 214) """
		return self.j * self.R_m()

	def eta_cathode(self):
		""" Cathode overpotential. (eq. 6.27, p. 209) (eq. 6.54, p. 215) """
		_pC_Pa = self.pC * ATM   # Convert atm -> Pa
		_tC_m = self.tC * 1e-9  # Convert nano-m -> m
		_Deff_O2_m2 = self.Deff_O2 * 1e-4  # Convert cm^2/s -> m^2/s

		_k = R * self.T / (4 * self.alpha * F)
		_f = (_tC_m * self.j * R * self.T) / (4 * F * _pC_Pa * _Deff_O2_m2)
		_a = self.j_0 * self.pC

		return _k * log(self.j / (_a * (self.x_O2 - _f * 1e-4))) * 1e-3  # TODO: number is off by factor of 1000

	def voltage(self):
		return self.E_thermo - self.eta_ohmic() - self.eta_cathode()


if __name__ == "__main__":
	fc = PEMFC_1D()
	print(fc.eta_cathode())
