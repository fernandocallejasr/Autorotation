import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
sns.set()

class AutoGiro():
    def __init__(self, AG_Data, hdens, list_Voo_Vh = [0], alpha = [0]):
        '''Clase que calcula la Potencia requerida por la hélice del motor  
           para un AutoGiro en vuelo horizontal a una misma altitud
           Inputs: 
               AG_Data (diccionario): Datos del AutoGiro
               hdens (float): Altitud densidad
               list_Voo_Vh (list): Lista que contiene las relaciones de velocidad
                                   Voo/Vh
               alpha (list): Lista que corresponde al ángulo de las relaciones de 
                             velocidad de la lista Voo/Vh
        '''
        
        self.AG_Data = AG_Data
        
        self.W  = AG_Data['W']    # [N]
        self.Pd = AG_Data['Pd']   # [hp]
        self.D  = AG_Data['D']    # [m]
        self.n  = AG_Data['n']    # [1]
        self.f  = AG_Data['f']    # [m^2]
        self.A  = np.pi*(1/4)*(self.D)**2  # [m^2]
        
        self.hdens = hdens        # [feet]
        self.dens   = 1.225*(1-self.hdens*(0.00357/518.4))**(1/0.235)   # [kg/m^3]
        self.dens_ratio = self.dens/1.225   # [1]
        
        self.list_Voo_Vh = list_Voo_Vh
        self. alpha = alpha
        
        self.list_Dp   = [0 for i in Voo_Vh]   # [N]
        self.list_Tp   = [0 for i in Voo_Vh]   # [N]
        self.list_Pshp = [0 for i in Voo_Vh]   # [hp]
        self.list_Voo  = [0 for i in Voo_Vh]   # [m/s]
        self.list_T    = [0 for i in Voo_Vh]   # [N]
        self.list_Pc   = [0 for i in Voo_Vh]   # [hp]
        self.list_Vc   = [0 for i in Voo_Vh]   # [m/s]
        self.DF   = pd.DataFrame({'Voo/Vh': self.list_Voo_Vh, 'alpha':self.alpha})
        
    def calc_T(self):
        """Método que calcula el empuje a cada Voo/Vh
        """
        for i in range(len(self.list_Voo_Vh)):
            self.list_T[i] = self.W/np.cos(np.deg2rad(self.alpha[i]))
        
    def calc_Voo(self):
        """Método que calcula el valor de la velocidad de avance, Voo, a la hdens dada
        """
        i = 0
        for V in self.list_Voo_Vh:
            self.list_Voo[i] = V*(self.list_T[i]/(2*self.dens*self.A))**0.5
            i += 1
        
    def calc_Dp(self):
        """Método que calcula el valor de Dp a las Voo y hdens
        """
        i = 0
        for V in self.list_Voo:
            self.list_Dp[i] = 0.5*self.dens*(V**2)*self.f
            i += 1
        
    def calc_Tp(self):
        """Método que calcula el valor de Tp a las Voo y hdens
        """
        i = 0
        for V in self.list_Voo:
            self.list_Tp[i] = self.list_T[i]*np.sin(np.deg2rad(alpha[i])) + self.list_Dp[i]
            i += 1
            
    def calc_Pshp(self):
        """Método que calcula el valor de Pshp (en realidad es P_req) a las Voo y hdens
        """
        i = 0
        for V in self.list_Voo:
            self.list_Pshp[i] = (self.list_Tp[i]*self.list_Voo[i])/745.7
            i += 1
            
    def calc_Pc(self):
        """Método que calcula el valor de Pc a las Voo y hdens
        """
        i = 0
        for V in self.list_Voo:
            self.list_Pc[i] = self.dens_ratio*self.Pd*self.n - self.list_Pshp[i]
            i += 1
            
    def calc_Vc(self):
        """Método que calcula el valor de Vc a las Voo y hdens
        """
        i = 0
        for V in self.list_Voo:
            self.list_Vc[i] = 745.7*self.list_Pc[i]/self.list_T[i]
            i += 1
            
    def calc_All(self):
        """Método que calcula T, Voo, Dp, Tp, Pshp, Pc, Vc para una hdens dada
        """
        self.__init__(self.AG_Data, self.hdens, self.list_Voo_Vh, self.alpha)
        self.calc_T()
        self.calc_Voo()
        self.calc_Dp()
        self.calc_Tp()
        self.calc_Pshp()
        self.calc_Pc()
        self.calc_Vc()
            
    def data_frame(self):
        """Método que calcula las potencias y muestra en un Data Frame a la hdens inicializada
        """
        self.calc_All()
        
        DF_P = pd.DataFrame({'Voo [m/s]': np.round(self.list_Voo, 2), 'Pshp [hp]': np.round(self.list_Pshp, 2), 'Pc [hp]': np.round(self.list_Pc, 2), 'Vc [m/s]': np.round(self.list_Vc, 2)})
        self.DF = self.DF.join(DF_P)

        return self.DF
    
    def calc_techo(self, alt_inicial):
        """Método que varía la altitud densidad hasta encontrar el techo de servicio
           sólo lo hace para una Voo/Vh = 1.8 & alpha = 30°
           Inputs:
               alt_inicial (float): Altitud para comenzar a iterar 
           Outputs:
               Techo de servicio
        """
        Error = 10
        alt_dens = alt_inicial
        
        while Error > 0.01:
            self.hdens = alt_dens
            self.calc_All()
            P_R    = self.list_Pshp[3]
            P_disp = self.dens_ratio*self.Pd*self.n
            Error  = P_disp - P_R
            alt_dens += 0.1

        return f'Techo de servicio: {alt_dens} ft'
        
    #def vel_sigma(self, sigma):
    #    for j in range(len(sigma)):
    #       i = 0
    #        for x in self.list_Voo_Vh:
    #            self.list_T[i] = self.W/np.cos(np.deg2rad(self.alpha[i]))
    #            self.list_Voo[i] = x*(self.list_T[i]/(2*1.225*sigma[j]*self.A))**0.5
    #            i += 1
        
    def calc_potencia(self, sigma = [1]):
        """Método que calcula Pshp & Pc para distintas velocidades Voo & alpha a la hdens dada
           También muestra la información en un Data Frame que incluye:
               Voo/Vh, alpha, Voo, Pshp, Pc & Vc
        """
        self.__init__(self.AG_Data, self.hdens, self.list_Voo_Vh, self.alpha)
        for j in range(len(sigma)):
            i = 0
            for x in self.list_Voo_Vh:
                self.list_T[i] = self.W/np.cos(np.deg2rad(self.alpha[i]))
                self.list_Voo[i] = x*(self.list_T[i]/(2*1.225*sigma[j]*self.A))**0.5
                self.list_Dp[i] = 0.5*1.225*sigma[j]*(self.list_Voo[i]**2)*self.f
                self.list_Tp[i] = self.list_T[i]*np.sin(np.deg2rad(alpha[i])) + self.list_Dp[i] 
                self.list_Pshp[i] = (self.list_Tp[i]*self.list_Voo[i])/745.7
                self.list_Pc[i] = sigma[j]*self.Pd*self.n - self.list_Pshp[i]
                self.list_Vc[i] = 745.7*self.list_Pc[i]/self.list_T[i]
                i += 1
            DF_P = pd.DataFrame({f'Voo_@s = {sigma[j]} [m/s]': np.round(self.list_Voo, 2), f'Pshp_@s = {sigma[j]} [hp]': np.round(self.list_Pshp, 2), f'Pc_@s = {sigma[j]} [hp]': np.round(self.list_Pc, 2), f'Vc_@s = {sigma[j]} [m/s]': np.round(self.list_Vc, 2)})
            self.DF = self.DF.join(DF_P)
            
        return self.DF


def main():
    Voo_Vh = [1.8,  1.8, 1.8, 1.8,  3,  4, 4.25, 4.52,  6,  6.4,    8,   9]
    alpha =   [80,   60,  40,  30, 10,  7,    6,    5,  4,    3,  3.5,   2]

    AG1_Data = {
               'W': 450*9.81,   #[N]
               'Pd': 100,       #[hp]
               'D' : 8.4,       #[m]
               'n' : 0.7,       #[1]
               'f' : 0.557418   #[m^2]
                }

    Ejemplo = AutoGiro(AG1_Data, 0, Voo_Vh, alpha)
    Ejemplo.data_frame()


import unittest

AG1_Data = {
           'W': 450*9.81,   #[N]
           'Pd': 100,       #[hp]
           'D' : 8.4,       #[m]
           'n' : 0.7,       #[1]
           'f' : 0.557418   #[m^2]
            }

Voo_Vh = [4]; alpha = [7]

class Test_bemClass(unittest.TestCase):
    def setUp(self):
        self.AGH = AutoGiro(AG1_Data, 0, Voo_Vh, alpha)
        self.AGH.calc_All()
        
    def test_initializaion(self):
        #self.assertEqual(self.clothing.color, 'orange', 'color should be orange')
        self.assertEqual(round(self.AGH.list_T[0], 2), 4447.65, 'T = 4447.65 N')
        self.assertEqual(round(self.AGH.A, 2) , 55.42, 'Área = 55.41 m^2')
        self.assertEqual(round(self.AGH.list_Voo[0], 2), 22.89, 'Voo = 22.890 m/s')
        self.assertEqual(round(self.AGH.list_Dp[0], 2), 178.95, 'Dp = 178.95 N')
        self.assertEqual(round(self.AGH.list_Tp[0], 2), 720.98, 'Tp = 720.98 N')
        self.assertEqual(round(self.AGH.list_Pshp[0], 2), 22.13, 'Pshp = 22.13')
        self.assertEqual(round(self.AGH.list_Pc[0], 2), 47.87, 'Pc = 47.87 hp')
        self.assertEqual(round(self.AGH.list_Vc[0], 2), 8.03, 'Vc = 8.03 m/s')
        
tests = Test_bemClass()
tests_loaded = unittest.TestLoader().loadTestsFromModule(tests)
unittest.TextTestRunner().run(tests_loaded)

if __name__ == '__main__':
    main()

