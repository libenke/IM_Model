import numpy as np
class IM_Model_precise:
    def __init__(self, λmax = 3.34, shear_rate = 31.6, finish_time = 10, τR = 0.124, β = 0.25, δt = 0.0001, CQ = 6,\
                    Gi = [1.14E2, 2.52, 8.80E-1, 3.81E-1, 2.07E-1, 1.70E-1, 1.76E-1, 1.37E-1, 2.03E-1], \
                    τi_eq = [2.16E-3, 5.90E-2, 3.33E-1, 1.53, 6.87, 2.9E1, 1.11E2, 3.63E2, 8.92E2]):
        self.node_points = 1E4
        self.λmax = λmax
        self.τR = τR
        self.β = β
        self.CQ = CQ
        ############
        self.δt = δt
        self.finish_time = finish_time
        self._δt = None #the value passed to function _IM_Tumbling_Multimode_integral()
        self._finish_time = None #the value passed to function _IM_Tumbling_Multimode_integral()
        self._finish_time_old = None
        ############
        self.Gi = np.array(Gi)
        self.τi_eq = np.array(τi_eq)
        self.shear_rate = shear_rate
        self.κ = np.array([[0,self.shear_rate,0],[0,0,0],[0,0,0]])
        self.t_span_log = np.logspace(np.log10(self.δt),np.log10(self.finish_time),1000)
        self.Q_span = np.zeros([len(self.t_span_log),3,3]) # nonlinear strain measure
        self.S_average_span = np.zeros([len(self.t_span_log),3,3]) # orientation tensor
        self.Si_span = np.zeros([len(self.τi_eq),len(self.t_span_log),3,3]) # orientation tensor
        self.τi_span = np.zeros([len(self.τi_eq),len(self.t_span_log)]) # tau_i
        self.τd_span = np.zeros(len(self.t_span_log)) # tau_d
        self.λ_span = np.zeros(len(self.t_span_log)) # λ

    def S_integral(self, Q_span,τd_span):
        """
        the Q_span,τd_span were trancate in the time t
        usage: S_integral(Q_span[:ind+1,:,:],δt,τd_span[:ind+1])
        """
        δt = self._δt
        #\int_t'^t dt''/taud(t'')
        τd_inv = 1/τd_span
        #exp[-\int_t'^t dt''/taud(t'')]
        exp_minus_int_tp_t_τd_inv = np.exp(-(np.cumsum(τd_inv[::-1])-τd_inv[-1]/2) * δt)[::-1]
        #integral from time=0 to time t
        S00 = np.sum(τd_inv * exp_minus_int_tp_t_τd_inv * Q_span[::-1,0,0]) * δt
        S01 = np.sum(τd_inv * exp_minus_int_tp_t_τd_inv * Q_span[::-1,0,1]) * δt
        S02 = np.sum(τd_inv * exp_minus_int_tp_t_τd_inv * Q_span[::-1,0,2]) * δt
        S11 = np.sum(τd_inv * exp_minus_int_tp_t_τd_inv * Q_span[::-1,1,1]) * δt
        S12 = np.sum(τd_inv * exp_minus_int_tp_t_τd_inv * Q_span[::-1,1,2]) * δt
        S22 = np.sum(τd_inv * exp_minus_int_tp_t_τd_inv * Q_span[::-1,2,2]) * δt
        S_0_t = np.array([[S00,S01,S02],\
                        [S01,S11,S12],\
                        [S02,S12,S22]])
        #plus the integral from time = 0 to time = -infinite
        #import ipdb
        #ipdb.set_trace()
        S = S_0_t + np.exp(-np.sum(1/τd_span)*δt) * Q_span[-1,:,:]
        return S
    
    def update_finish_time(self):
        """initial the values for calculation,
        and increase the _finish_time 2 orders of magnitudes
        """
        if self._finish_time is None:
            #initialize the values if _finish_time is None.
            self._δt = self.δt
            self._finish_time = self._δt * self.node_points
            if self._finish_time >= self.finish_time*0.95:
                self._finish_time = self.finish_time
            self._finish_time_old = None
            #initialize the matrix used in the calculations
            self._t_span = np.arange(0,self._finish_time,self._δt)
            self._Q_span = np.zeros([len(self._t_span),3,3]) # nonlinear strain measure
            self._S_average_span = np.zeros([len(self._t_span),3,3]) # orientation tensor
            self._Si_span = np.zeros([len(self.τi_eq),len(self._t_span),3,3]) # orientation tensor
            self._τi_span = np.zeros([len(self.τi_eq),len(self._t_span)]) # tau_i
            self._τd_span = np.zeros(len(self._t_span)) # tau_d
            self._λ_span = np.zeros(len(self._t_span)) # λ
        else:
            #update the matrix stored in log time matrix
            self._finish_time_old = self._finish_time
            _t_span_log_mask = self.t_span_log<=self._finish_time_old
            _t_span_log = self.t_span_log[_t_span_log_mask]
            self.λ_span[_t_span_log_mask] = np.interp(_t_span_log,self._t_span,self._λ_span)
            self.τd_span[_t_span_log_mask] = np.interp(_t_span_log,self._t_span,self._τd_span)
            for i in [0, 1, 2]:
                for j in [0, 1, 2]:
                    self.Q_span[_t_span_log_mask,i,j] = np.interp(_t_span_log, self._t_span, self._Q_span[:,i,j])
                    self.S_average_span[_t_span_log_mask,i,j] = np.interp(_t_span_log, self._t_span, self._S_average_span[:,i,j])
                    for k,_ in enumerate(self.τi_eq):
                        self.Si_span[k,_t_span_log_mask,i,j] = np.interp(_t_span_log, self._t_span, self._Si_span[k,:,i,j])
            for k,_ in enumerate(self.τi_eq):
                self.τi_span[k,_t_span_log_mask] = np.interp(_t_span_log, self._t_span, self._τi_span[k,:])
            
            #increase _finish_time 1 orders of magnitudes
            if self._finish_time == self.finish_time:
                return
            self._δt = self._δt*10
            _finish_time_new = self._δt * self.node_points
            if _finish_time_new >= self.finish_time*0.95:
                self._finish_time = self.finish_time
            else:
                self._finish_time = _finish_time_new
            
            #interp back from the log time matrix to the matrix in the calculation
            self._t_span = np.arange(0,self._finish_time,self._δt)
            _t_span_mask = self._t_span<=self._finish_time_old
            _t_span_interp = self._t_span[_t_span_mask]
            self._Q_span = np.zeros([len(self._t_span),3,3]) # nonlinear strain measure
            self._S_average_span = np.zeros([len(self._t_span),3,3]) # orientation tensor
            self._Si_span = np.zeros([len(self.τi_eq),len(self._t_span),3,3]) # orientation tensor
            for i in [0,1,2]:
                for j in [0,1,2]:
                    self._Q_span[_t_span_mask,i,j] = np.interp(_t_span_interp,_t_span_log,self.Q_span[_t_span_log_mask,i,j])
                    self._S_average_span[_t_span_mask,i,j] = np.interp(_t_span_interp,_t_span_log,self.S_average_span[_t_span_log_mask,i,j])
                    for k,_ in enumerate(self.τi_eq):
                        self._Si_span[k,_t_span_mask,i,j] = np.interp(_t_span_interp,_t_span_log,self.Si_span[k,_t_span_log_mask,i,j])
            self._τi_span = np.zeros([len(self.τi_eq),len(self._t_span)]) # tau_i
            for k,_ in enumerate(self.τi_eq):
                self._τi_span[k,_t_span_mask] = np.interp(_t_span_interp,_t_span_log,self.τi_span[k,_t_span_log_mask])
            self._τd_span = np.zeros(len(self._t_span)) # tau_d
            self._τd_span[_t_span_mask] = np.interp(_t_span_interp,_t_span_log,self.τd_span[_t_span_log_mask])
            self._λ_span = np.zeros(len(self._t_span)) # λ
            self._λ_span[_t_span_mask] = np.interp(_t_span_interp,_t_span_log,self.λ_span[_t_span_log_mask])
        return
    
    def calculate_stress(self):
        self.σ_span = np.zeros([len(self.t_span_log),3,3]) # σ stress
        self.Fλ = np.zeros_like(self.t_span_log)
        λmax = self.λmax
        for ind,_ in enumerate(self.t_span_log):
            #initial the value at time = 0s
            λ=self.λ_span[ind]
            Fλ_update = (λmax**2 - λ**2/3)/(λmax**2 - λ**2)*(λmax**2 - 1)/(λmax**2 - 1/3)
            self.Fλ[ind] = Fλ_update
            σ = np.zeros([3,3])
            for i in np.arange(len(self.τi_eq)):
                σ += self.CQ * Fλ_update * λ**2 * self.Gi[i] * self.Si_span[i,ind,:,:]
            self.σ_span[ind] = σ
        return self.t_span_log, self.σ_span
    
    def IM_Tumbling_Multimode_integral(self):
        #init
        self._finish_time = None
        while (self._finish_time is None) or (self._finish_time < self.finish_time * 0.98):
            import time
            t0 = time.time()
            self.update_finish_time()
            self._IM_Tumbling_Multimode_integral()
            print(f"self._δt: {self._δt}s, self._finish_time: {self._finish_time} \n")
            print(f"spend time: {time.time()-t0:.3f}s\n")
        self.update_finish_time()
        return self.calculate_stress()
    
    def _IM_Tumbling_Multimode_integral(self):
        """
        return: t_span, τd_span, σ_span, λ_span, S_span
        """
        #
        τR = self.τR
        β = self.β
        shear_rate = self.shear_rate
        λmax = self.λmax
        κ = self.κ
        Gi = self.Gi
        τi_eq = self.τi_eq

        δt = self._δt
        t_span = self._t_span
        Q_span = self._Q_span
        S_average_span = self._S_average_span
        Si_span = self._Si_span
        τi_span = self._τi_span
        τd_span = self._τd_span
        λ_span = self._λ_span
        
        #
        τi_t = np.zeros_like(self.τi_eq) #τi at time t
        WiR = τR * shear_rate
        ω = WiR**(-0.2)/(8 * np.pi) * shear_rate
        β_angle = WiR**(-0.2)/8 * shear_rate
        φ_span = np.cos(2 * np.pi * ω * t_span) * np.exp( - β_angle * t_span)
        r = 0 # tube loss rate
        for ind,t in enumerate(t_span):
            #initial the value at time = 0s 
            if (self._finish_time_old is not None) and (t < self._finish_time_old * 0.98):
                continue

            if ind == 0:
                Q_span[0,:,:] = np.eye(3)/3
                S_average_span[0,:,:] = np.eye(3)/3
                for i in np.arange(len(τi_eq)):
                    Si_span[i,0,:,:] = np.eye(3)/3
                    τi_span[i,0] = τi_eq[i]
                τd_span[0] = np.sum(Gi * τi_span[:,0]**2)/np.sum(Gi * τi_span[:,0])
                λ_span[0] = 1
            
            #calculate the Q at time t.
            γ_history = t * shear_rate
            E = np.eye(3)
            E[0,1] = γ_history
            B = np.dot(E,E.T)
            eigs, eig_vec = np.linalg.eig(B)
            B_square_root = np.dot(eig_vec,np.dot(np.diag(np.sqrt(eigs)),eig_vec.T))
            Q_span[ind,:,:] = B_square_root/np.trace(B_square_root)
            
            #calculate the τi(t)
            for i in np.arange(len(τi_eq)):
                #Si at time t
                Si_t = Si_span[i,ind-1,:,:]
                τi_t[i] = 1/(1/τi_eq[i] + β*r)
                τi_span[i,ind] = τi_t[i]
                #Si at time t update
                Si_t_update = self.S_integral(Q_span[:ind+1,:,:],τi_span[i,:ind+1])
                Si_span[i,ind,:,:] = Si_t_update
            #τd = \sum Giτi^2 / (\sum Giτi)
            τd_span[ind] = np.sum(Gi * τi_span[:,ind]**2)/np.sum(Gi * τi_span[:,ind])
            #S_average at time t
            S_average_t_update = self.S_integral(Q_span[:ind+1,:,:],τd_span[:ind+1])
            S_average_span[ind,:,:] = S_average_t_update
            
            #calculate the dλ_dt
            λ = λ_span[ind-1]
            Fλ = (λmax**2 - λ**2/3)/(λmax**2 - λ**2)*(λmax**2 - 1)/(λmax**2 - 1/3)
            φ = φ_span[ind]
            dλ_dt = φ * λ * np.trace(np.dot(κ,S_average_t_update)) - (λ * Fλ - 1) / τR
            λ_update = λ + dλ_dt * δt
            r =  np.trace(np.dot(κ,S_average_t_update)) - 1/λ_update*dλ_dt # update r, tube loss rate
            λ_span[ind] = λ_update
        return