pip install pymc

import pymc as pm
from pymc.step_methods import Metropolis
from pymc import Poisson
from pymc import Normal,Uniform
from pymc import sample, Model
from pymc import find_MAP
from pymc import sample_posterior_predictive
from arviz import plot_trace, plot_autocorr, summary

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import poisson


#포아송 회귀를 위한 synthetic data 생성

from scipy.stats import poisson

np.random.seed(321)
X1 = np.random.randn(50).reshape(-1,1)
X2 = np.random.randn(50).reshape(-1,1)
eps = np.random.normal(0,0.0001,50)
X = np.concatenate((X1,X2),axis=1)
Beta = np.array([0.6,-0.3])
lamb = np.exp(np.matmul(X,Beta)+eps)
Y = np.zeros(50)
for i in range(len(Y)):
  Y[i] = poisson.rvs(lamb[i],size=1)[0]
data = pd.DataFrame({'X1':X1.reshape(-1),'X2':X2.reshape(-1),'Y':Y})

data.head(5)



#PyMC 기반 포아송 회귀문제 정의 및 MCMC 샘플링 설정

with Model() as python_ex:
  beta1 = Normal('beta1',0,2)
  beta2 = Normal('beta2',0,2)

  Lk = Poisson('Y',mu=np.exp(beta1*data['X1']+beta2*data['X2']),observed=data['Y'])
with python_ex:
  results = sample(step=pm.HamiltonianMC(),draws=1500,chains=2, tune=500, random_seed=123)

#arviz 패키지로 HMC 기반 MCMC 샘플링 결과 시각화
plot_trace(results,compact=True,var_names=('beta1','beta2'),chain_prop={'color': ['darkcyan', 'red']})

#arviz 패키지로 r_hat, HDI 확인
summary(results)
