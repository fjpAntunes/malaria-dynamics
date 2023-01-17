

scenario_1 = dict(#Valor dos Parametros
beta_vh = 0.2, # 0 - 1 per mosquito # ok
kappa = 0.05, # 1/11 dimensionless ok
beta_hv = 0.2, # Per Mosquito ok
mu_M = 0.089, # 0.089 - 0.476 per day ok
nu = 1/15.6, # 1/15.6 +- 2.86 per day ok
alpha = 43.66,# 83 +- 48 larvae/per female mosquito ok
mu_L = 0.62, #0.62 - 0.99 per day - Larvae mortality on ponds without fish ok
mu_p = 0.31, # ok

#Parametros da Vegetação

r = 0.5/30, # 0.5 Por Mês ok
gamma = 0.65, # Proporção de vegetação retirada na limpeza ok
tau = 60, # 30-60 dias ok
H = 0.05, # Pop. realizando Limpeza - <5% ok

#Valor das Carrying Capacities Maximas
K_p_max = 0.4, #ok
K_w_max = 4,   # ok
K_0 = 0.5)     #ok

scenario_2 = dict(#Valor dos Parametros
beta_vh = 0.1, # 0 - 1 per mosquito #
kappa = 0.05, # 1/11 dimensionless
beta_hv = 0.1, # Per Mosquito
mu_M = 0.16, # 0.089 - 0.476
nu = 1/15.6, # 1/15.6 +- 2.86 per day
alpha = 43.66, # 83 +- 48 larvae/per female mosquito
mu_L = 0.62, #0.62 - 0.99 per day - Larvae mortality on ponds without fish
mu_p = 0.31,

#Parametros da Vegetação

r = 0.5/30, # 0.5 Por Mês 
gamma = 0.65, # Proporção de vegetação retirada na limpeza
tau = 45, # 30-60 dias
H = 0.05, # Pop. realizando Limpeza - <5%

#Valor das Carrying Capacities Maximas
K_p_max = 0.4,
K_w_max = 4,
K_0 = 0.5)

scenario_3 = dict(#Valor dos Parametros
beta_vh = 0.2, # 0 - 1 per mosquito #
kappa = 0.05, # 1/11 dimensionless
beta_hv = 0.2, # Per Mosquito
mu_p = 0.31, #0.99,
mu_M = 0.8, # 0.16 - 0.23 per day
mu_L = 0.99, #0.62 - 0.99 per day - Larvae mortality on ponds without fish

nu = 1/(15.6 + 2.86), # 1/15.6 +- 2.86 per day


alpha = 8.75, # 83 +- 48 larvae/per female mosquito


#Parametros da Vegetação

r = 0.5/30, # 0.5 Por Mês 
gamma = 0.65, # Proporção de vegetação retirada na limpeza
tau = 30,  # 30-60 dias
H = 0.05, # Pop. realizando Limpeza - <5%

#Valor das Carrying Capacities Maximas
K_p_max = 0.4,
K_w_max = 4,
K_0 = 0.5)