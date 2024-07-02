import numpy as np

def S_integral(Q_span,δt,t,τd_span):
    """
    the Q_span,τd_span were trancate in the time t
    usage: S_integral(Q_span[:ind+1,:,:],δt,τd_span[:ind+1])
    """
    #\int_t'^t dt''/taud(t'')
    τd_inv = 1/τd_span
    #exp[-\int_t'^t dt''/taud(t'')]
    exp_minus_int_tp_t_τd_inv = np.exp(-np.cumsum(τd_inv[::-1]) * δt)[::-1]
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
    S = S_0_t + np.exp(-np.sum(1/τd_span)*δt) * Q_span[-1,:,:]
    return S

def IM_Multimode_integral(λmax = 3.34, τR = 0.124, shear_rate = 31.6, δt = 0.001, finish_time = 10, \
            Gi = np.array([1.14E2, 2.52, 8.80E-1, 3.81E-1, 2.07E-1, 1.70E-1, 1.76E-1, 1.37E-1, 2.03E-1]), \
            τi_eq = np.array([2.16E-3, 5.90E-2, 3.33E-1, 1.53, 6.87, 2.9E1, 1.11E2, 3.63E2, 8.92E2]), CQ = 6):
    """
    return t_span, τd_span, σ_span, λ_span, S_average_span
    """
    κ = np.array([[0,shear_rate,0],\
              [0,0,0],\
              [0,0,0]])
    t_span = np.arange(0,finish_time,δt)
    #
    Q_span = np.zeros([len(t_span),3,3]) # nonlinear strain measure
    S_average_span = np.zeros([len(t_span),3,3]) # orientation tensor
    Si_span = np.zeros([len(τi_eq),len(t_span),3,3]) # orientation tensor
    τi_span = np.zeros([len(τi_eq),len(t_span)]) # tau_i
    τd_span = np.zeros(len(t_span)) # tau_d
    λ_span = np.zeros(len(t_span)) # λ
    σ_span = np.zeros([len(t_span),3,3]) # σ stress
    τi_t = np.zeros_like(τi_eq) #τi at time t
    for ind,t in enumerate(t_span):
        #initial the value at time = 0s 
        if ind == 0:
            Q_span[0,:,:] = np.eye(3)/3
            S_average_span[0,:,:] = np.eye(3)/3
            for i in np.arange(len(τi_eq)):
                Si_span[i,0,:,:] = np.eye(3)/3
                if τi_eq[i] > τR:
                    τi_span[i,0] = τi_eq[i]/2  + τR
                else:
                    τi_span[i,0] = τi_eq[i]/2
            τd_span[0] = np.sum(Gi * τi_span[:,0]**2)/np.sum(Gi * τi_span[:,0])
            λ_span[0] = 1
            for i in np.arange(len(τi_eq)):
                σ_span[0,:,:] += CQ * Gi[i] * Si_span[i,0,:,:]
            continue
        
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
            if τi_eq[i] > τR:
                τi_t[i] = 1/2/(1/τi_eq[i] + np.trace(np.dot(κ,Si_t))) + τR
            else:
                τi_t[i] = 1/2/(1/τi_eq[i] + np.trace(np.dot(κ,Si_t)))
            #Si at time t update
            Si_t_update = S_integral(Q_span[:ind,:,:],δt,t,τi_span[i,:ind])
            Si_span[i,ind,:,:] = Si_t_update
        τi_span[:,ind] = τi_t
        #τd = \sum Giτi^2 / (\sum Giτi)
        τd_span[ind] = np.sum(Gi * τi_span[:,ind]**2)/np.sum(Gi * τi_span[:,ind])
        #S_average at time t
        S_average_t_update = S_integral(Q_span[:ind,:,:],δt,t,τd_span[:ind])
        S_average_span[ind,:,:] = S_average_t_update
        
        #calculate the dλ_dt
        λ = λ_span[ind-1]
        Fλ = (λmax**2 - λ**2/3)/(λmax**2 - λ**2)*(λmax**2 - 1)/(λmax**2 - 1/3)
        dλ_dt = λ * np.trace(np.dot(κ,S_average_t_update)) - (λ * Fλ - 1) / τR
        λ_update = λ + dλ_dt * δt
        λ_span[ind] = λ_update
        Fλ_update = (λmax**2 - λ_update**2/3)/(λmax**2 - λ_update**2)*(λmax**2 - 1)/(λmax**2 - 1/3)
        σ_update = np.zeros([3,3])
        for i in np.arange(len(τi_eq)):
            σ_update += CQ * Fλ_update * λ_update**2 * Gi[i] * Si_span[i,ind,:,:]
        σ_span[ind] = σ_update
    return t_span, τd_span, σ_span, λ_span, S_average_span


def IM_SingleMode_integral(λmax = 3.34, τR = 0.124, shear_rate = 31.6, δt = 0.001, \
                           finish_time = 10, GN0 = 2.79E5, τd_eq = 0.8):
    """
    return: t_span, τd_span, σ_span, λ_span, S_span
    """
    #
    κ = np.array([[0,shear_rate,0],\
                  [0,0,0],\
                  [0,0,0]])
    t_span = np.arange(0,finish_time,δt)
    #
    Q_span = np.zeros([len(t_span),3,3]) # nonlinear strain measure
    S_span = np.zeros([len(t_span),3,3]) # orientation tensor
    τd_span = np.zeros(len(t_span)) # tau_d
    λ_span = np.zeros(len(t_span)) # λ
    σ_span = np.zeros([len(t_span),3,3]) # σ stress
    for ind,t in enumerate(t_span):
        # initial the value at time = 0s 
        if ind == 0:
            Q_span[0,:,:] = np.eye(3)/3
            S_span[0,:,:] = np.eye(3)/3
            τd_span[0] = τd_eq/2 + τR
            λ_span[0] = 1
            σ_span[0,:,:] = 3 * GN0 * np.eye(3)/3
            continue
        
        # calculate the Q at time t.
        γ_history = t * shear_rate
        E = np.eye(3)
        E[0,1] = γ_history
        B = np.dot(E,E.T)
        eigs, eig_vec = np.linalg.eig(B)
        B_square_root = np.dot(eig_vec,np.dot(np.diag(np.sqrt(eigs)),eig_vec.T))
        Q_span[ind,:,:] = B_square_root/np.trace(B_square_root)
        
        # calculate the τd(t)
        # S at time t
        S_t = S_span[ind-1,:,:]
        τd_t = 1/2/(1/τd_eq + np.trace(np.dot(κ,S_t))) + τR
        τd_span[ind] = τd_t
        #S at time t update
        S_update = S_integral(Q_span[:ind+1,:,:],δt,t,τd_span[:ind+1])
        S_span[ind,:,:] = S_update
        
        # calculate the dλ_dt
        λ = λ_span[ind-1]
        Fλ = (λmax**2 - λ**2/3)/(λmax**2 - λ**2)*(λmax**2 - 1)/(λmax**2 - 1/3)
        dλ_dt = λ * np.trace(np.dot(κ,S_update)) - (λ * Fλ - 1) / τR
        λ_update = λ + dλ_dt * δt
        λ_span[ind] = λ_update
        Fλ_update = (λmax**2 - λ_update**2/3)/(λmax**2 - λ_update**2) * (λmax**2 - 1)/(λmax**2 - 1/3)
        σ_update = 3 * GN0 * Fλ_update * λ_update**2 * S_update
        σ_span[ind] = σ_update
    return t_span, τd_span, σ_span, λ_span, S_span


def IM_SingleMode_differential(λmax = 3.34, τR = 0.124, shear_rate = 31.6, δt = 0.001, \
                               finish_time = 10, GN0 = 2.79E5, τd_eq = 0.8):
    """
    return: t_span, τd_span, σ_span, λ_span, S_span
    """
    #
    κ = np.array([[0,shear_rate,0],\
              [0,0,0],\
              [0,0,0]])
    t_span = np.arange(0,finish_time,δt)
    #
    S_span = np.zeros([len(t_span),3,3]) # orientation tensor S
    S2_span = np.zeros([len(t_span),3,3]) # orientation tensor S**2
    τd_span = np.zeros(len(t_span)) # tau_d
    λ_span = np.zeros(len(t_span)) # λ
    σ_span = np.zeros([len(t_span),3,3]) # σ stress
    for ind,t in enumerate(t_span):
        # initial the value at time = 0s 
        if ind == 0:
            S_span[0,:,:] = np.eye(3)/3
            S2_span[0,:,:] = np.eye(3)/9
            τd_span[0] = τd_eq/2 + τR
            λ_span[0] = 1
            σ_span[0,:,:] = 3 * GN0 * np.eye(3)/3
            continue
        # calculate the τd(t)
        S_t = S_span[ind-1,:,:] #S at time t
        S2_t = S2_span[ind-1,:,:] #S2 at time t
        τd_t = 1/2/(1/τd_eq + np.trace(np.dot(κ,S_t))) + τR
        τd_span[ind] = τd_t
        
        # dS2_dt at time t update
        dS2_dt = np.dot(κ, S2_t) + np.dot(S2_t, κ.T) - \
                2 * S2_t * np.trace(np.dot(κ, S_t)) - \
                2/τd_t * np.dot(S_t, S_t-np.eye(3)/3)
        S2_update = S2_t + dS2_dt * δt
        S2_span[ind,:,:] = S2_update
        eigs, eig_vec = np.linalg.eig(S2_update)
        S_update = np.dot(eig_vec,np.dot(np.diag(np.sqrt(eigs)),eig_vec.T))
        S_span[ind,:,:] = S_update
        
        # calculate the dλ_dt
        λ = λ_span[ind-1]
        Fλ = (λmax**2 - λ**2/3)/(λmax**2 - λ**2)*(λmax**2 - 1)/(λmax**2 - 1/3)
        dλ_dt = λ * np.trace(np.dot(κ,S_update)) - (λ * Fλ - 1) / τR
        λ_update = λ + dλ_dt * δt
        λ_span[ind] = λ_update
        Fλ_update = (λmax**2 - λ_update**2/3)/(λmax**2 - λ_update**2) * (λmax**2 - 1)/(λmax**2 - 1/3)
        σ_update = 3 * GN0 * Fλ_update * λ_update**2 * S_update
        σ_span[ind] = σ_update
    return t_span, τd_span, σ_span, λ_span, S_span


def IM_Tumbling_Multimode_integral(λmax = 3.34, τR = 0.124, shear_rate = 31.6, δt = 0.001, finish_time = 10, \
                Gi = np.array([1.14E2, 2.52, 8.80E-1, 3.81E-1, 2.07E-1, 1.70E-1, 1.76E-1, 1.37E-1, 2.03E-1]), \
                τi_eq = np.array([2.16E-3, 5.90E-2, 3.33E-1, 1.53, 6.87, 2.9E1, 1.11E2, 3.63E2, 8.92E2]), CQ = 6):
    """
    return: t_span, τd_span, σ_span, λ_span, S_span
    """
    #
    κ = np.array([[0,shear_rate,0],\
              [0,0,0],\
              [0,0,0]])
    t_span = np.arange(0,finish_time,δt)
    #
    Q_span = np.zeros([len(t_span),3,3]) # nonlinear strain measure
    S_average_span = np.zeros([len(t_span),3,3]) # orientation tensor
    Si_span = np.zeros([len(τi_eq),len(t_span),3,3]) # orientation tensor
    τi_span = np.zeros([len(τi_eq),len(t_span)]) # tau_i
    τd_span = np.zeros(len(t_span)) # tau_d
    λ_span = np.zeros(len(t_span)) # λ
    σ_span = np.zeros([len(t_span),3,3]) # σ stress
    τi_t = np.zeros_like(τi_eq) #τi at time t
    WiR = τR * shear_rate
    ω = WiR**(-0.2)/(8 * np.pi) * shear_rate
    β = WiR**(-0.2)/8 * shear_rate
    φ_span = np.cos(2 * np.pi * ω * t_span) * np.exp( - β * t_span)
    for ind,t in enumerate(t_span):
        #initial the value at time = 0s 
        if ind == 0:
            Q_span[0,:,:] = np.eye(3)/3
            S_average_span[0,:,:] = np.eye(3)/3
            for i in np.arange(len(τi_eq)):
                Si_span[i,0,:,:] = np.eye(3)/3
                if τi_eq[i] > τR:
                    τi_span[i,0] = τi_eq[i]/2  + τR
                else:
                    τi_span[i,0] = τi_eq[i]/2
            τd_span[0] = np.sum(Gi * τi_span[:,0]**2)/np.sum(Gi * τi_span[:,0])
            λ_span[0] = 1
            for i in np.arange(len(τi_eq)):
                σ_span[0,:,:] += CQ * Gi[i] * Si_span[i,0,:,:]
            continue
        
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
            if τi_eq[i] > τR:
                τi_t[i] = 1/2/(1/τi_eq[i] + np.trace(np.dot(κ,Si_t))) + τR
            else:
                τi_t[i] = 1/2/(1/τi_eq[i] + np.trace(np.dot(κ,Si_t)))
            #Si at time t update
            Si_t_update = S_integral(Q_span[:ind,:,:],δt,t,τi_span[i,:ind])
            Si_span[i,ind,:,:] = Si_t_update
        τi_span[:,ind] = τi_t
        #τd = \sum Giτi^2 / (\sum Giτi)
        τd_span[ind] = np.sum(Gi * τi_span[:,ind]**2)/np.sum(Gi * τi_span[:,ind])
        #S_average at time t
        S_average_t_update = S_integral(Q_span[:ind,:,:],δt,t,τd_span[:ind])
        S_average_span[ind,:,:] = S_average_t_update
        
        #calculate the dλ_dt
        λ = λ_span[ind-1]
        Fλ = (λmax**2 - λ**2/3)/(λmax**2 - λ**2)*(λmax**2 - 1)/(λmax**2 - 1/3)
        φ = φ_span[ind]
        dλ_dt = φ * λ * np.trace(np.dot(κ,S_average_t_update)) - (λ * Fλ - 1) / τR
        λ_update = λ + dλ_dt * δt
        λ_span[ind] = λ_update
        Fλ_update = (λmax**2 - λ_update**2/3)/(λmax**2 - λ_update**2)*(λmax**2 - 1)/(λmax**2 - 1/3)
        σ_update = np.zeros([3,3])
        for i in np.arange(len(τi_eq)):
            σ_update += CQ * Fλ_update * λ_update**2 * Gi[i] * Si_span[i,ind,:,:]
        σ_span[ind] = σ_update
    return t_span, τd_span, σ_span, λ_span, S_average_span