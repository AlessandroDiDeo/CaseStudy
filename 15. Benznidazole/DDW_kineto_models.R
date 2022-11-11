#### original
ddw_kineto_mod1 <- function(t, S, parameters) {
	with(as.list(c(S, parameters)), {

		# Intracellular parasite duplication 
		N1 = N				
		N2 = N * exp(BETA)
		N3 = N * exp(2*BETA)
		N4 = N * exp(3*BETA)
		N5 = N * exp(4*BETA)

		# Amastigote cell number is the effect of intracellular duplication 
		S8 = (N1*S[2] + N2*S[4] + N3*S[5] + N4*S[6] + N5*S[7])
    
		# Total number of cells
		STOT = S[1] + S[2] + S[4] + S[5] + S[6] + S[7]
		
		
		## drug effect
	
		CONC = S[10]
		
		IKIN = t^5/(TEKIN^5 + t^5)
    IK1 = t^5/(TEK1^5 + t^5)
    IKOUT = t^5/(TEKOUT^5 + t^5)
		EKD = EMAXKD*t^5/(TEKD^5 + t^5)
		
		EKIN    = IKIN*(CONC)/(IC50KIN + CONC)
		EK1     = IK1*(CONC)/(IC50K1 + CONC)
		EKOUT   = IKOUT*(CONC)/(IC50KOUT + CONC)
		EKDEATH = EKD*CONC/(IC50KD + CONC)

		KOUT  	= KOUT0 * (1 - EKOUT)	      # Infected cells burst rate
		KPROT  	= KPROT0 				            # Cells infection rate
		KDEATH 	= KDEATH0 * (1 + EKDEATH)		# Amastigotes death rate

		## Delay rate for cells
		K1 	= K1_1 * (1 - EK1)
		K3  = 0.5 * K1

		## Cell growth
		R     = GR1
		BASE  = (S[1] + S[2] + S[4] + S[5] + S[6] + S[7])/NMAX
		KIN   = R * (1 - (BASE)^GAMMA)*(1 - EKIN)
		
		## PK
    KEL   = CL/V

		### ODE system
		dSdt=vector(len=10)

		# Healthy cells compartment
		dSdt[1] = KIN*S[1] - KPROT*S[3]*S[1]

		# Infected cells compartment
		dSdt[2] = KPROT*S[1]*S[3] - K3*S[2] - KOUT*S[2]

		# Delayed infected cells stage 1
		dSdt[4] = K3*S[2] - K1*S[4] - KOUT*S[4]

		# Delayed infected cells stage 2
		dSdt[5] = K1*S[4] - K1*S[5] - KOUT*S[5]

		# Delayed infected cells stage 3
		dSdt[6] = K1*S[5] - K1*S[6] - KOUT*S[6]

		# Delayed infected cells stage 4
		dSdt[7] = K1*S[6]

		# Trypomastigotes compartment
		dSdt[3] = KOUT*S[8] - KOUT*N5*dSdt[7] - KPROT*S[1]*S[3]

		# Amastigotes compartment
		dSdt[8] = (N1*dSdt[2] + N2*dSdt[4] + N3*dSdt[5] + N4*dSdt[6] + N5*dSdt[7]) - KDEATH*S[8]

		## PK model (set events for dosing events)
		dSdt[9]  = - KA*S[9]
		dSdt[10] = KA*S[9] - KEL*S[10]

		
		list(dSdt)
	})
}